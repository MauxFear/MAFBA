#!/usr/bin/env python3
"""
Kcat Monte Carlo sampling (all reactions).

Quantifies how uncertainty in enzyme turnover numbers (kcat) affects
acetate threshold (λ_ac), maximum growth rate (μ_max), and flux allocations.
Outputs are written to the directory passed as output_dir (e.g. data/outputs/Figure_S5/all_reactions).
"""

import cobra
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging
from datetime import datetime
from scipy import stats
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 300

class KcatMonteCarloAnalysis:
    """Monte Carlo analysis of kcat uncertainty on MAFBA model predictions."""
    
    def __init__(self, model_path, kcat_data_path, output_dir, solver='glpk'):
        """
        Initialize analysis.
        
        Parameters:
        -----------
        model_path : str
            Path to SBML model file
        kcat_data_path : str
            Path to Excel file with kcat data
        output_dir : str
            Directory for outputs
        solver : str
            Solver to use ('gurobi', 'glpk', 'cplex', etc.)
        """
        self.model_path = Path(model_path)
        self.kcat_data_path = Path(kcat_data_path)
        self.output_dir = Path(output_dir)
        self.solver = solver
        
        # Create a date-labeled subfolder for the output
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.output_dir = self.output_dir / timestamp
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Create output directories
        self.figures_dir = self.output_dir / 'figures'
        self.logs_dir = self.output_dir / 'logs'
        for d in [self.figures_dir, self.logs_dir]:
            d.mkdir(parents=True, exist_ok=True)
        
        # Initialize data containers
        self.model = None
        self.kcat_df = None
        self.protein_constraint_id = 'Prot_Const'
        self.membrane_constraint_id = 'Memb_Const'
        self.baseline_membrane_bound = None
        self.baseline_costs = {}
        self.factors = {}
        self.results = []
        
        logger.info(f"Initialized analysis with output dir: {self.output_dir}")
        logger.info(f"Solver: {self.solver}")
    
    def load_model(self):
        """Load SBML model and identify protein constraint."""
        logger.info(f"Loading model from {self.model_path}")
        self.model = cobra.io.read_sbml_model(str(self.model_path))
        
        # Set solver
        try:
            self.model.solver = self.solver
            logger.info(f"Set solver to: {self.solver}")
        except Exception as e:
            logger.warning(f"Could not set solver to {self.solver}: {e}")
            logger.warning(f"Using default solver: {self.model.solver}")
        
        logger.info(f"Model loaded: {len(self.model.reactions)} reactions, "
                   f"{len(self.model.metabolites)} metabolites")
        
        # For MAFBA model, the protein constraint is 'Protein_Allocation'
        if 'Protein_Allocation' in self.model.metabolites:
            self.protein_constraint_id = 'Protein_Allocation'
            protein_met = self.model.metabolites.get_by_id(self.protein_constraint_id)
            logger.info(f"Identified protein constraint: {self.protein_constraint_id}")
            logger.info(f"  Protein constraint appears in {len(protein_met.reactions)} reactions")
        else:
            # Fallback: search for protein-related metabolites
            protein_candidates = [m for m in self.model.metabolites 
                                if 'protein' in m.id.lower() or 'protein' in m.name.lower()]
            
            if protein_candidates:
                self.protein_constraint_id = protein_candidates[0].id
                logger.info(f"Identified protein constraint (fallback): {self.protein_constraint_id}")
            else:
                logger.error("Could not identify protein constraint metabolite")
                raise ValueError("Protein constraint metabolite not found")
        
        # Store baseline membrane constraint bound
        if self.membrane_constraint_id in self.model.reactions:
            membrane_rxn = self.model.reactions.get_by_id(self.membrane_constraint_id)
            self.baseline_membrane_bound = membrane_rxn.upper_bound
            logger.info(f"Identified membrane constraint: {self.membrane_constraint_id}")
            logger.info(f"  Baseline upper bound: {self.baseline_membrane_bound}")
        else:
            logger.warning(f"Membrane constraint {self.membrane_constraint_id} not found")
            self.baseline_membrane_bound = None
    
    def load_kcat_data(self, kcat_threshold=1e7):
        """
        Load and filter kcat data from Excel file.
        
        Parameters:
        -----------
        kcat_threshold : float
            Filter out kcat values above this threshold (default 1e7)
        """
        logger.info(f"Loading kcat data from {self.kcat_data_path}")
        self.kcat_df = pd.read_excel(self.kcat_data_path, sheet_name='Reactions')
        
        logger.info(f"Loaded {len(self.kcat_df)} reactions")
        logger.info(f"kcat statistics before filtering:\n{self.kcat_df['kcat'].describe()}")
        
        # Filter out extreme values
        original_count = len(self.kcat_df)
        self.kcat_df = self.kcat_df[self.kcat_df['kcat'] < kcat_threshold].copy()
        filtered_count = len(self.kcat_df)
        
        logger.info(f"Filtered kcat values >= {kcat_threshold}: "
                   f"{original_count - filtered_count} reactions removed")
        logger.info(f"Remaining reactions: {filtered_count}")
        logger.info(f"kcat statistics after filtering:\n{self.kcat_df['kcat'].describe()}")
    
    def load_sector_data(self):
        """Load proteome sector data."""
        logger.info("Loading proteome sector data")
        self.sector_df = pd.read_excel(self.kcat_data_path, sheet_name='Reactions')
        logger.info(f"Loaded {len(self.sector_df)} reactions from sector data")
    
    def extract_baseline_costs(self):
        """
        Extract baseline proteome costs from model S matrix.
        
        Costs are the coefficients in the protein constraint equation.
        For MAFBA models: cost = factor / kcat
        """
        logger.info("Extracting baseline costs from protein constraint")
        
        if self.protein_constraint_id is None:
            logger.error("Protein constraint not identified")
            return
        
        protein_met = self.model.metabolites.get_by_id(self.protein_constraint_id)
        
        # Extract coefficients from S matrix for protein constraint
        # Store both the cost and the sign for restoration
        for rxn in protein_met.reactions:
            coeff = rxn.metabolites[protein_met]
            if coeff != 0:
                self.baseline_costs[rxn.id] = coeff  # Store with sign
        
        logger.info(f"Extracted baseline costs for {len(self.baseline_costs)} reactions")
        
        # Calculate factors (factor = cost * kcat)
        for rxn_id, cost in self.baseline_costs.items():
            # Find kcat for this reaction
            kcat_row = self.kcat_df[self.kcat_df['rxns'] == rxn_id]
            if not kcat_row.empty:
                kcat = kcat_row.iloc[0]['kcat']
                self.factors[rxn_id] = abs(cost) * kcat  # Use absolute value for factor
        
        logger.info(f"Calculated factors for {len(self.factors)} reactions")
        
        # Log statistics
        if self.factors:
            factors_array = np.array(list(self.factors.values()))
            logger.info(f"Factor statistics:\n{pd.Series(factors_array).describe()}")
    
    def restore_baseline_costs(self):
        """
        Restore baseline costs to the model.
        
        This is called before each Monte Carlo iteration to ensure
        each iteration starts from the same baseline state.
        """
        protein_met = self.model.metabolites.get_by_id(self.protein_constraint_id)
        
        for rxn_id, baseline_coeff in self.baseline_costs.items():
            rxn = self.model.reactions.get_by_id(rxn_id)
            current_coeff = rxn.metabolites.get(protein_met, 0)
            if current_coeff != 0:
                # Restore to baseline
                rxn.add_metabolites({protein_met: -current_coeff + baseline_coeff})
    
    def set_gamma(self, gamma):
        """
        Set membrane constraint multiplier (gamma).
        
        Parameters:
        -----------
        gamma : float
            Multiplier for membrane constraint upper bound
            gamma = 1.0 means baseline constraint
            gamma < 1.0 means tighter membrane constraint
        """
        if self.baseline_membrane_bound is None:
            logger.warning("Baseline membrane bound not available, skipping gamma adjustment")
            return
        
        membrane_rxn = self.model.reactions.get_by_id(self.membrane_constraint_id)
        new_bound = gamma * self.baseline_membrane_bound
        membrane_rxn.upper_bound = new_bound
        logger.info(f"Set gamma = {gamma:.2f}, membrane bound = {new_bound:.6f}")
    
    def update_model_costs(self, kcat_sample):
        """
        Update model costs based on sampled kcat values.
        
        Note: Updates costs in-place on the original model to avoid
        Python 3.13 deepcopy issues with COBRApy.
        
        Parameters:
        -----------
        kcat_sample : dict
            Dictionary mapping reaction IDs to sampled kcat values
        """
        protein_met = self.model.metabolites.get_by_id(self.protein_constraint_id)
        
        for rxn_id, new_kcat in kcat_sample.items():
            if rxn_id in self.factors:
                # Calculate new cost: cost = factor / kcat
                new_cost = self.factors[rxn_id] / new_kcat
                
                # Update coefficient in S matrix
                rxn = self.model.reactions.get_by_id(rxn_id)
                old_coeff = rxn.metabolites.get(protein_met, 0)
                if old_coeff != 0:
                    # Preserve sign of coefficient
                    sign = 1 if self.baseline_costs[rxn_id] > 0 else -1
                    rxn.add_metabolites({protein_met: -old_coeff + sign * new_cost})
    
    def run_glucose_scan(self, model, glucose_uptake_values):
        """
        Run growth optimization for a range of glucose uptake values.
        
        Parameters:
        -----------
        model : cobra.Model
            Model to optimize
        glucose_uptake_values : list
            List of glucose uptake rates to scan
        
        Returns:
        --------
        list
            List of dictionaries with optimization results for each uptake value
        """
        scan_results = []
        
        glc_exchange_rxn_id = 'EX_glc__D_e_b'
        if glc_exchange_rxn_id not in model.reactions:
            logger.error(f"Glucose exchange reaction '{glc_exchange_rxn_id}' not found.")
            return scan_results

        for uptake_val in glucose_uptake_values:
            with model:
                model.reactions.get_by_id(glc_exchange_rxn_id).upper_bound = uptake_val
                
                results = {
                    'glucose_uptake_bound': uptake_val,
                    'solver_status': None,
                    'objective_value': np.nan,
                    'growth_rate': np.nan,
                    'glucose_uptake': np.nan,
                    'acetate_secretion': np.nan,
                    'oxygen_uptake': np.nan,
                    'co2_production': np.nan,
                    'atp_production': np.nan,
                }
                
                try:
                    # Optimize
                    solution = model.optimize()
                    
                    results['solver_status'] = solution.status
                    results['objective_value'] = solution.objective_value
                    results['growth_rate'] = solution.objective_value
                    
                    # Extract key fluxes
                    flux_map = {
                        'glucose_uptake': 'EX_glc__D_e_b',
                        'acetate_secretion': 'EX_ac_e_f',
                        'oxygen_uptake': 'EX_o2_e_b',
                        'co2_production': 'EX_co2_e_f',
                        'atp_production': 'ATPS4rpp_f',
                    }
                    
                    for key, rxn_id in flux_map.items():
                        if rxn_id in model.reactions:
                            results[key] = solution.fluxes.get(rxn_id, 0.0)
                    
                    # Store all fluxes
                    results['fluxes'] = solution.fluxes.to_dict()

                except Exception as e:
                    logger.warning(f"Optimization failed for glucose uptake {uptake_val}: {e}")
                    results['solver_status'] = 'failed'
                
                scan_results.append(results)
                
        return scan_results
    
    def monte_carlo_sampling(self, n_iterations=500, cv=0.1, gamma_values=[1.0]):
        """
        Run Monte Carlo sampling over kcat values for multiple gamma values.
        
        Parameters:
        -----------
        n_iterations : int
            Number of Monte Carlo iterations per gamma value
        cv : float
            Coefficient of variation for log-normal sampling (e.g., 0.1 = 10% uncertainty)
            This represents the measurement uncertainty in kcat values
        gamma_values : list
            List of gamma values (membrane constraint multipliers) to test
        """
        logger.info(f"Starting Monte Carlo sampling: {n_iterations} iterations per gamma")
        logger.info(f"Coefficient of Variation (CV): {cv:.2%} ({cv:.3f})")
        logger.info(f"Gamma values: {gamma_values}")
        
        # Define glucose uptake values to scan
        points_below_10 = np.linspace(0.2, 9.8, 15)
        points_above_10 = 10 + 25 * (np.linspace(0, 1, 25)**2)
        glucose_uptake_values = np.concatenate([points_below_10, points_above_10])

        # Get reactions with both kcat data and factors
        valid_rxns = [rxn for rxn in self.kcat_df['rxns'] if rxn in self.factors]
        logger.info(f"Sampling {len(valid_rxns)} reactions with valid kcat and factors")
        
        # Initialize result dataframes
        self.flux_df = pd.DataFrame()
        self.aceT_df = pd.DataFrame()
        self.sectors_df = pd.DataFrame()
        self.coeffs_df = pd.DataFrame()
        self.kcats_df = pd.DataFrame()

        # Loop over gamma values
        for gamma in tqdm(gamma_values, desc="Gamma values"):
            logger.info(f"\n{'='*50}")
            logger.info(f"Processing gamma = {gamma:.2f}")
            logger.info(f"{'='*50}")
            
            # Set gamma value for membrane constraint
            self.set_gamma(gamma)
            
            # Monte Carlo iterations for this gamma
            for i in tqdm(range(1, n_iterations + 1), desc=f"γ={gamma:.2f} iterations", leave=False):
                # Restore baseline costs before each iteration
                self.restore_baseline_costs()
                
                # Sample kcat values using FIXED CV (removes Layer 1 randomness)
                # Only lognormal sampling provides Monte Carlo variability
                sigma = np.sqrt(np.log(cv**2 + 1))
                kcat_sample = {}
                for rxn_id in valid_rxns:
                    orig_kcat = self.kcat_df[self.kcat_df['rxns'] == rxn_id].iloc[0]['kcat']
                    mu = np.log(orig_kcat) - 0.5 * sigma**2
                    new_kcat = np.random.lognormal(mu, sigma)
                    kcat_sample[rxn_id] = new_kcat
                
                # Store kcat sample with gamma label
                kcat_sample_df = pd.DataFrame([kcat_sample])
                kcat_sample_df['iteration'] = i
                kcat_sample_df['gamma'] = gamma
                self.kcats_df = pd.concat([self.kcats_df, kcat_sample_df], ignore_index=True)

                # Update model with sampled kcats (in-place)
                self.update_model_costs(kcat_sample)
                
                # Store coefficients
                coeffs_sample = {}
                protein_met = self.model.metabolites.get_by_id(self.protein_constraint_id)
                for rxn in self.model.reactions:
                    if protein_met in rxn.metabolites:
                        coeffs_sample[rxn.id] = rxn.metabolites[protein_met]
                coeffs_sample_df = pd.DataFrame([coeffs_sample])
                coeffs_sample_df['iteration'] = i
                coeffs_sample_df['gamma'] = gamma
                self.coeffs_df = pd.concat([self.coeffs_df, coeffs_sample_df], ignore_index=True)

                # Run glucose scan
                scan_results = self.run_glucose_scan(self.model, glucose_uptake_values)
                
                # Process scan results
                iter_flux_df = pd.DataFrame([res['fluxes'] for res in scan_results])
                iter_flux_df['iteration'] = i
                iter_flux_df['gamma'] = gamma
                iter_flux_df['glucose_uptake_bound'] = [res['glucose_uptake_bound'] for res in scan_results]
                self.flux_df = pd.concat([self.flux_df, iter_flux_df], ignore_index=True)
                
                # Calculate acetate threshold
                # Sanitize growth rates: set to 0 if None or negative
                growth_rates = [max(0, res['growth_rate']) if res['growth_rate'] is not None else 0 
                               for res in scan_results]
                # Sanitize acetate secretion: set to 0 if None or negative
                acetate_secretions = [max(0, res['acetate_secretion']) if res['acetate_secretion'] is not None else 0 
                                     for res in scan_results]
                
                # Simple threshold detection: first point with acetate secretion > 1e-6
                acetate_threshold = 0
                max_growth = max(growth_rates) if growth_rates else 0
                for j in range(len(acetate_secretions)):
                    if acetate_secretions[j] > 1e-6:
                        acetate_threshold = growth_rates[j]
                        break
                
                aceT_row = pd.DataFrame([{
                    'iteration': i, 
                    'gamma': gamma,
                    'acetate_threshold': acetate_threshold,
                    'max_growth_rate': max_growth
                }])
                self.aceT_df = pd.concat([self.aceT_df, aceT_row], ignore_index=True)

                # Process sectors data
                sector_fluxes = {}
                for sector in self.sector_df['reactionProtGroups'].unique():
                    sector_rxns = self.sector_df[self.sector_df['reactionProtGroups'] == sector]['rxns']
                    sector_fluxes[sector] = iter_flux_df[sector_rxns].sum(axis=1)
                
                sector_fluxes_df = pd.DataFrame(sector_fluxes)
                sector_fluxes_df['iteration'] = i
                sector_fluxes_df['gamma'] = gamma
                sector_fluxes_df['glucose_uptake_bound'] = [res['glucose_uptake_bound'] for res in scan_results]
                self.sectors_df = pd.concat([self.sectors_df, sector_fluxes_df], ignore_index=True)
                
                if i % 50 == 0:
                    logger.info(f"Gamma {gamma:.2f}: Completed {i}/{n_iterations} iterations")
            
            logger.info(f"Completed gamma = {gamma:.2f}")
        
        logger.info("Monte Carlo sampling completed for all gamma values")
        
        # Save iteration results
        self.flux_df.to_csv(self.output_dir / 'flux_data.csv', index=False)
        self.aceT_df.to_csv(self.output_dir / 'aceT_data.csv', index=False)
        self.coeffs_df.to_csv(self.output_dir / 'coeffs_data.csv', index=False)
        self.kcats_df.to_csv(self.output_dir / 'kcats_data.csv', index=False)
        self.sectors_df.to_csv(self.output_dir / 'sectors_data.csv', index=False)
        logger.info("Saved iteration results to CSV files.")
    
    def compute_summary_statistics(self):
        """Compute and save summary statistics grouped by gamma."""
        logger.info("Computing summary statistics")

        # Acetate threshold summary by gamma
        aceT_summary = self.aceT_df.groupby('gamma')[['acetate_threshold', 'max_growth_rate']].describe()
        aceT_summary.to_csv(self.output_dir / 'aceT_summary.csv')
        logger.info("Saved acetate threshold summary.")

        # Summary for relevant reactions
        relevant_reactions = [
            'ATPS4rpp_f', 'PGI_f', 'FBA_f', 'TKT1_f', 'GND', 'CS', 
            'ACONTa_f', 'ICDHyr_f', 'AKGDH', 'SUCDi', 'FUM_f', 'MDH_f'
        ]
        
        # Calculate mean flux for each iteration and gamma
        mean_fluxes = self.flux_df.groupby(['gamma', 'iteration'])[relevant_reactions].mean().reset_index()
        
        # Summary of mean fluxes by gamma
        flux_summary = mean_fluxes.groupby('gamma')[relevant_reactions].describe()
        flux_summary.to_csv(self.output_dir / 'flux_summary.csv')
        logger.info("Saved flux summary for relevant reactions.")

        # Summary of k-cat values for relevant reactions by gamma
        kcats_summary = self.kcats_df.groupby('gamma')[relevant_reactions].describe()
        kcats_summary.to_csv(self.output_dir / 'kcats_summary.csv')
        logger.info("Saved k-cat summary for relevant reactions.")
    
    def run_full_analysis(self, n_iterations=500, kcat_threshold=1e7, cv=0.1, gamma_values=[1.0]):
        """
        Run complete analysis pipeline with gamma sweep.
        
        Parameters:
        -----------
        n_iterations : int
            Number of Monte Carlo iterations per gamma value (default 500)
        kcat_threshold : float
            Filter kcat values above this threshold (default 1e7)
        cv : float
            Coefficient of variation for kcat sampling (default 0.1 = 10% uncertainty)
            Can be specified as decimal (0.1) or will interpret percentage (10 → 0.1)
        gamma_values : list
            List of gamma values (membrane constraint multipliers) to test (default [1.0])
        """
        # Convert percentage to decimal if needed (e.g., 10 → 0.1)
        if cv > 1.0:
            logger.info(f"Converting CV from percentage: {cv}% → {cv/100:.3f}")
            cv = cv / 100.0
        
        logger.info("="*60)
        logger.info("KCAT MONTE CARLO SAMPLING (ALL REACTIONS)")
        logger.info("="*60)
        logger.info(f"Coefficient of Variation: {cv:.2%} (uncertainty in kcat measurements)")
        logger.info(f"Gamma values: {gamma_values}")
        logger.info(f"Iterations per gamma: {n_iterations}")
        logger.info(f"Total optimizations: ~{len(gamma_values) * n_iterations * 40}")
        
        start_time = datetime.now()
        
        # Step 1: Load model
        self.load_model()
        
        # Step 2: Load kcat data
        self.load_kcat_data(kcat_threshold=kcat_threshold)
        
        # Step 2a: Load sector data
        self.load_sector_data()
        
        # Step 3: Extract baseline costs
        self.extract_baseline_costs()
        
        # Step 4: Run Monte Carlo sampling with gamma sweep
        self.monte_carlo_sampling(
            n_iterations=n_iterations, 
            cv=cv,
            gamma_values=gamma_values
        )
        
        # Step 5: Compute statistics
        self.compute_summary_statistics()
        
        # Save metadata
        metadata = {
            'analysis_name': 'kcat_sampling_all_reactions',
            'date': datetime.now().isoformat(),
            'solver': self.solver,
            'n_iterations_per_gamma': n_iterations,
            'gamma_values': gamma_values,
            'n_gamma_values': len(gamma_values),
            'total_iterations': len(gamma_values) * n_iterations,
            'kcat_threshold': kcat_threshold,
            'cv': cv,
            'cv_percentage': f"{cv*100:.1f}%",
            'model_path': str(self.model_path),
            'kcat_data_path': str(self.kcat_data_path),
            'n_reactions_sampled': len([rxn for rxn in self.kcat_df['rxns'] if rxn in self.factors]),
            'baseline_membrane_bound': self.baseline_membrane_bound,
            'runtime_seconds': (datetime.now() - start_time).total_seconds(),
        }
        
        with open(self.output_dir / 'analysis_metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)
        
        end_time = datetime.now()
        runtime = (end_time - start_time).total_seconds()
        
        logger.info("="*60)
        logger.info(f"ANALYSIS COMPLETED in {runtime:.1f} seconds ({runtime/60:.1f} minutes)")
        logger.info(f"Total iterations: {len(gamma_values) * n_iterations}")
        logger.info(f"Results saved to: {self.output_dir}")
        logger.info("="*60)

