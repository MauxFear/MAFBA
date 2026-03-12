#!/usr/bin/env python3
"""
Protein export gamma sweep analysis.

Simulates the effect of dedicating a fraction of the ribosomal proteome budget
to protein translocation machinery (Sec/Tat) on acetate secretion and proteome
sector allocation across a range of membrane constraint multipliers (γ).

Outputs (written to output_dir/):
    flux_data.csv        – per-point flux solutions
    sectors_data.csv     – per-point proteome sector totals
    aceL_data.csv        – acetate-threshold (λac) detection rows
    coeffs_data.csv      – membrane allocation coefficients per export fraction
    kcats_data.csv       – scenario-level summary (μmax, λac, bounds)
    analysis_metadata.json
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
from tqdm import tqdm
import warnings
import os

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 300

BIOMASS_RXN = 'BIOMASS_Ec_iML1515_WT_75p37M'
GLUCOSE_UPTAKE = 'EX_glc__D_e_b'
ACETATE_SECRETION = 'EX_ac_e_f'

PROTEIN_CONSTRAINT = 'Prot_Const'
MEMBRANE_CONSTRAINT = 'Memb_Const'

PROTEIN_ALLOCATION = 'Protein_Allocation'
MEMBRANE_ALLOCATION = 'Membrane_Allocation'

TRANSLOCATION_KEYWORDS = [
    'secy', 'sece', 'secf', 'secg', 'seca', 'secb', 'secd',
    'tata', 'tatb', 'tatc', 'tate', 'tat', 'yidc', 'ffh', 'ftsy',
    'transloc', 'translocase', 'signal recognition particle'
]


class ProteinExportGammaSweep:
    """
    Gamma-sweep analysis with translocation machinery export cost.

    Applies a fraction of the ribosomal proteome cost as an additional
    membrane burden and sweeps over γ values to compute acetate thresholds
    and proteome sector allocation.
    """

    def __init__(self, model_path, attributes_path, output_dir, solver='gurobi'):
        self.model_path = Path(model_path)
        self.attributes_path = Path(attributes_path)
        self.solver = solver

        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.output_dir = Path(output_dir) / timestamp
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.figures_dir = self.output_dir / 'figures'
        self.logs_dir = self.output_dir / 'logs'
        for d in [self.figures_dir, self.logs_dir]:
            d.mkdir(parents=True, exist_ok=True)

        self.base_model = None
        self.ribosome_cost_raw = None
        self.translocation_cost = 0.0
        self.ribosome_cost_adjusted = None
        self.translocation_rxns = []
        self.baseline_membrane_coeff = None
        self.sector_mapping = None
        self.phi_max = None

        self.flux_data_path = None
        self.sectors_data_path = None
        self.aceL_data_path = None
        self.coeffs_data_path = None
        self.kcats_data_path = None

        logger.info(f"Output directory: {self.output_dir}")

    def load_model(self):
        logger.info(f"Loading model from {self.model_path}")
        self.base_model = cobra.io.read_sbml_model(str(self.model_path))

        try:
            self.base_model.solver = self.solver
        except Exception as e:
            logger.warning(f"Could not set solver to {self.solver}: {e}")

        self._verify_model()
        self._compute_ribosome_cost()

        biomass = self.base_model.reactions.get_by_id(BIOMASS_RXN)
        membrane_met = self.base_model.metabolites.get_by_id(MEMBRANE_ALLOCATION)
        self.baseline_membrane_coeff = biomass.metabolites.get(membrane_met, 0.0)
        logger.info(f"Baseline membrane coefficient: {self.baseline_membrane_coeff:.6f}")

        self.phi_max = self.base_model.reactions.get_by_id(PROTEIN_CONSTRAINT).upper_bound
        logger.info(f"φ_max: {self.phi_max:.4f}")

        self._load_sector_mapping()
        return self

    def _verify_model(self):
        required_metabolites = [PROTEIN_ALLOCATION, MEMBRANE_ALLOCATION]
        required_reactions = [BIOMASS_RXN, GLUCOSE_UPTAKE, PROTEIN_CONSTRAINT, MEMBRANE_CONSTRAINT]

        for met_id in required_metabolites:
            assert met_id in [m.id for m in self.base_model.metabolites], \
                f"Required metabolite {met_id} not found"

        for rxn_id in required_reactions:
            assert rxn_id in [r.id for r in self.base_model.reactions], \
                f"Required reaction {rxn_id} not found"

        solution = self.base_model.optimize()
        assert solution.status == 'optimal', f"Baseline optimization failed: {solution.status}"
        logger.info(f"Baseline growth rate: {solution.objective_value:.4f} hr⁻¹")

    def _compute_ribosome_cost(self):
        df = pd.read_excel(self.attributes_path, sheet_name='Reactions')
        mapping = dict(zip(df['rxns'], df['reactionProtGroups']))

        protein_met = self.base_model.metabolites.get_by_id(PROTEIN_ALLOCATION)
        ribosome_cost = 0.0
        for rxn in protein_met.reactions:
            if mapping.get(rxn.id) == 'phiR':
                ribosome_cost += abs(rxn.metabolites[protein_met])

        if ribosome_cost == 0:
            raise ValueError("Ribosome cost is zero – phiR reactions not found")

        self.ribosome_cost_raw = ribosome_cost
        logger.info(f"Ribosome cost: {ribosome_cost:.6f}")

        self._estimate_translocation_cost()
        self.ribosome_cost_adjusted = max(self.ribosome_cost_raw - self.translocation_cost, 0.0)
        logger.info(
            "Ribosome cost adjusted: raw=%.6f, translocation=%.6f, adjusted=%.6f",
            self.ribosome_cost_raw, self.translocation_cost, self.ribosome_cost_adjusted
        )

    def _estimate_translocation_cost(self):
        override = os.getenv("TRANSLOCATION_COST_OVERRIDE")
        protein_met = self.base_model.metabolites.get_by_id(PROTEIN_ALLOCATION)

        if override:
            try:
                self.translocation_cost = float(override)
                logger.info("Using environment override for translocation cost: %.6f", self.translocation_cost)
                return
            except ValueError:
                logger.error("Invalid TRANSLOCATION_COST_OVERRIDE; using heuristic")

        matches = [r for r in self.base_model.reactions
                   if any(k in f"{r.id} {r.name or ''}".lower() for k in TRANSLOCATION_KEYWORDS)]

        if not matches:
            logger.warning("No translocation reactions detected; translocation cost set to 0")
            self.translocation_cost = 0.0
            self.translocation_rxns = []
            return

        trans_cost = 0.0
        tracked = []
        for rxn in matches:
            coeff = rxn.metabolites.get(protein_met, 0.0)
            if abs(coeff) > 0:
                trans_cost += abs(coeff)
                tracked.append(rxn.id)

        self.translocation_cost = trans_cost
        self.translocation_rxns = tracked
        logger.info("Translocation cost: %.6f from %d reactions", trans_cost, len(tracked))

    def _load_sector_mapping(self):
        df = pd.read_excel(self.attributes_path, sheet_name='Reactions')
        self.sector_mapping = dict(zip(df['rxns'], df['reactionProtGroups'].fillna('unassigned')))

    def _prepare_output_files(self):
        reaction_ids = [r.id for r in self.base_model.reactions]

        self.flux_data_path = self.output_dir / 'flux_data.csv'
        pd.DataFrame(columns=['export_fraction', 'gamma', 'glc_uptake', 'mu',
                               'objective_status', 'is_lambda_ac'] + reaction_ids
                     ).to_csv(self.flux_data_path, index=False)

        self.sectors_data_path = self.output_dir / 'sectors_data.csv'
        pd.DataFrame(columns=['export_fraction', 'gamma', 'glc_uptake', 'mu',
                               'phiEc', 'phiEm', 'phiEr', 'phiT', 'phiCm', 'phiR', 'phiCc', 'unassigned']
                     ).to_csv(self.sectors_data_path, index=False)

        self.aceL_data_path = self.output_dir / 'aceL_data.csv'
        pd.DataFrame(columns=['export_fraction', 'gamma', 'lambda_ac_mu', 'lambda_ac_glucose',
                               'lambda_ac_acetate', 'phi_M_at_lambda', 'notes']
                     ).to_csv(self.aceL_data_path, index=False)

        self.coeffs_data_path = self.output_dir / 'coeffs_data.csv'
        pd.DataFrame(columns=['export_fraction'] + reaction_ids).to_csv(self.coeffs_data_path, index=False)

        self.kcats_data_path = self.output_dir / 'kcats_data.csv'
        pd.DataFrame(columns=['export_fraction', 'gamma', 'mu_max', 'phi_Mmax_bound',
                               'phi_Pmax_bound', 'lambda_ac_mu', 'lambda_ac_acetate', 'notes']
                     ).to_csv(self.kcats_data_path, index=False)

    def run_analysis(self, export_fractions, gamma_values, acetate_threshold=0.01):
        logger.info("STARTING PROTEIN EXPORT GAMMA SWEEP")
        logger.info(f"Export fractions: {export_fractions}")
        logger.info(f"Gamma values: {gamma_values}")

        self._prepare_output_files()
        glucose_grid = self._generate_glucose_grid()
        lambda_ac_recorded = set()

        for export_fraction in tqdm(export_fractions, desc="Export fractions"):
            export_cost = export_fraction * self.ribosome_cost_adjusted
            logger.info(f"Export fraction: {export_fraction:.2f}, cost: {export_cost:.6f}")

            model = self.base_model.copy()
            self._apply_export_cost(model, export_cost)
            self._save_membrane_coefficients(model, export_fraction)

            for gamma in tqdm(gamma_values, desc=f"  γ sweep (f={export_fraction:.2f})", leave=False):
                model.reactions.get_by_id(MEMBRANE_CONSTRAINT).upper_bound = gamma * self.phi_max
                lambda_ac_data = self._run_glucose_sweep(
                    model, export_fraction, gamma, glucose_grid,
                    acetate_threshold, lambda_ac_recorded
                )
                self._save_scenario_metadata(export_fraction, gamma, lambda_ac_data)

        self._save_analysis_metadata(export_fractions, gamma_values, acetate_threshold, glucose_grid)
        logger.info(f"Analysis complete. Results in: {self.output_dir}")

    def _generate_glucose_grid(self):
        coarse = np.linspace(0.2, 9.8, 15)
        dense = 10 + 25 * (np.linspace(0, 1, 25) ** 2)
        return np.concatenate([coarse, dense])

    def _apply_export_cost(self, model, export_cost):
        biomass = model.reactions.get_by_id(BIOMASS_RXN)
        membrane_met = model.metabolites.get_by_id(MEMBRANE_ALLOCATION)
        desired = self.baseline_membrane_coeff + export_cost
        current = biomass.metabolites.get(membrane_met, 0.0)
        delta = desired - current
        if abs(delta) > 1e-12:
            biomass.add_metabolites({membrane_met: delta})

    def _save_membrane_coefficients(self, model, export_fraction):
        membrane_met = model.metabolites.get_by_id(MEMBRANE_ALLOCATION)
        row = {'export_fraction': export_fraction}
        row.update({r.id: r.metabolites.get(membrane_met, 0.0) for r in model.reactions})
        pd.DataFrame([row]).to_csv(self.coeffs_data_path, mode='a', header=False, index=False)

    def _run_glucose_sweep(self, model, export_fraction, gamma, glucose_grid,
                           acetate_threshold, lambda_ac_recorded):
        lambda_ac_data = None
        key = (export_fraction, gamma)

        glc_rxn = model.reactions.get_by_id(GLUCOSE_UPTAKE)
        biomass_rxn = model.reactions.get_by_id(BIOMASS_RXN)
        prot_rxn = model.reactions.get_by_id(PROTEIN_CONSTRAINT)

        for glc_uptake in glucose_grid:
            glc_rxn.lower_bound = 0.0
            glc_rxn.upper_bound = glc_uptake

            try:
                model.objective = biomass_rxn
                model.objective_direction = "max"
                solution = model.optimize()

                if solution.status == 'optimal':
                    mu_opt = solution.objective_value
                    orig_lb = biomass_rxn.lower_bound
                    biomass_rxn.lower_bound = mu_opt
                    model.objective = prot_rxn
                    model.objective_direction = "min"
                    pfba_sol = model.optimize()
                    if pfba_sol.status == 'optimal':
                        solution = pfba_sol
                    biomass_rxn.lower_bound = orig_lb
                    model.objective = biomass_rxn
                    model.objective_direction = "max"

                if solution.status == 'optimal':
                    fluxes = solution.fluxes
                    mu = fluxes.get(BIOMASS_RXN, np.nan)
                    acetate_flux = fluxes.get(ACETATE_SECRETION, 0)
                    is_lambda_ac = False

                    if key not in lambda_ac_recorded and acetate_flux >= acetate_threshold:
                        is_lambda_ac = True
                        lambda_ac_recorded.add(key)
                        phi_M = fluxes.get(MEMBRANE_CONSTRAINT, np.nan)
                        lambda_ac_data = {
                            'export_fraction': export_fraction, 'gamma': gamma,
                            'lambda_ac_mu': mu, 'lambda_ac_glucose': glc_uptake,
                            'lambda_ac_acetate': acetate_flux, 'phi_M_at_lambda': phi_M,
                            'notes': 'Acetate onset detected'
                        }
                        pd.DataFrame([lambda_ac_data]).to_csv(
                            self.aceL_data_path, mode='a', header=False, index=False)
                else:
                    mu = np.nan
                    fluxes = pd.Series(np.nan, index=model.reactions.list_attr('id'))
                    is_lambda_ac = False

            except Exception as e:
                logger.error(f"Optimization failed at glc={glc_uptake}: {e}")
                mu = np.nan
                fluxes = pd.Series(np.nan, index=model.reactions.list_attr('id'))
                solution = None
                is_lambda_ac = False

            self._save_flux_data(export_fraction, gamma, glc_uptake, mu,
                                 solution.status if solution else 'failed', is_lambda_ac, fluxes)
            self._save_sector_data(export_fraction, gamma, glc_uptake, mu, fluxes)

        return lambda_ac_data

    def _save_flux_data(self, export_fraction, gamma, glc_uptake, mu, status, is_lambda_ac, fluxes):
        row = {'export_fraction': export_fraction, 'gamma': gamma,
               'glc_uptake': glc_uptake, 'mu': mu,
               'objective_status': status, 'is_lambda_ac': is_lambda_ac}
        row.update(fluxes.to_dict())
        pd.DataFrame([row]).to_csv(self.flux_data_path, mode='a', header=False, index=False)

    def _save_sector_data(self, export_fraction, gamma, glc_uptake, mu, fluxes):
        sector_totals = {
            'phiEc': 0.0, 'phiEm': 0.0, 'phiEr': 0.0, 'phiT': 0.0,
            'phiCm': 0.0, 'phiR': 0.0, 'phiCc': 0.0, 'unassigned': 0.0
        }
        protein_met = self.base_model.metabolites.get_by_id(PROTEIN_ALLOCATION)

        for rxn_id, flux in fluxes.items():
            if pd.notna(flux) and abs(flux) > 1e-9:
                if rxn_id in self.base_model.reactions:
                    rxn = self.base_model.reactions.get_by_id(rxn_id)
                    cost = rxn.metabolites.get(protein_met, 0.0)
                    frac = abs(flux * cost)
                sector = self.sector_mapping.get(rxn_id, 'unassigned')
                if sector in sector_totals:
                    sector_totals[sector] += frac

        phiR_raw = sector_totals.get('phiR', 0.0)
        if self.translocation_cost > 0 and pd.notna(mu):
            phiT_val = min(max(mu * self.translocation_cost, 0.0), phiR_raw)
            sector_totals['phiT'] = phiT_val
            sector_totals['phiR'] = max(phiR_raw - phiT_val, 0.0)

        row = {'export_fraction': export_fraction, 'gamma': gamma,
               'glc_uptake': glc_uptake, 'mu': mu, **sector_totals}
        pd.DataFrame([row]).to_csv(self.sectors_data_path, mode='a', header=False, index=False)

    def _save_scenario_metadata(self, export_fraction, gamma, lambda_ac_data):
        flux_df = pd.read_csv(self.flux_data_path)
        mu_max = flux_df[
            (flux_df['export_fraction'] == export_fraction) &
            (flux_df['gamma'] == gamma)
        ]['mu'].max()

        if lambda_ac_data is None:
            meta_row = {
                'export_fraction': export_fraction, 'gamma': gamma, 'mu_max': mu_max,
                'phi_Mmax_bound': gamma * self.phi_max, 'phi_Pmax_bound': self.phi_max,
                'lambda_ac_mu': np.nan, 'lambda_ac_acetate': np.nan,
                'notes': 'Acetate threshold not reached'
            }
        else:
            meta_row = {
                'export_fraction': export_fraction, 'gamma': gamma, 'mu_max': mu_max,
                'phi_Mmax_bound': gamma * self.phi_max, 'phi_Pmax_bound': self.phi_max,
                'lambda_ac_mu': lambda_ac_data['lambda_ac_mu'],
                'lambda_ac_acetate': lambda_ac_data['lambda_ac_acetate'],
                'notes': lambda_ac_data['notes']
            }
        pd.DataFrame([meta_row]).to_csv(self.kcats_data_path, mode='a', header=False, index=False)

    def _save_analysis_metadata(self, export_fractions, gamma_values, acetate_threshold, glucose_grid):
        metadata = {
            'timestamp': datetime.now().isoformat(),
            'model_path': str(self.model_path),
            'attributes_path': str(self.attributes_path),
            'solver': self.solver,
            'parameters': {
                'export_fractions': export_fractions,
                'gamma_values': gamma_values,
                'acetate_threshold': acetate_threshold,
                'n_glucose_points': len(glucose_grid),
                'glucose_range': [float(glucose_grid.min()), float(glucose_grid.max())]
            },
            'model_info': {
                'ribosome_cost': float(self.ribosome_cost_adjusted),
                'ribosome_cost_raw': float(self.ribosome_cost_raw),
                'translocation_cost': float(self.translocation_cost),
                'baseline_membrane_coeff': float(self.baseline_membrane_coeff),
                'phi_max': float(self.phi_max),
                'n_reactions': len(self.base_model.reactions),
                'n_metabolites': len(self.base_model.metabolites)
            },
            'output_directory': str(self.output_dir)
        }
        with open(self.output_dir / 'analysis_metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)
