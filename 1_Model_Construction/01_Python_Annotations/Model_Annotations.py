#%%

import pandas as pd
import numpy as np
import cobra
import cameo
import json
import os
from six import iteritems, string_types
import pprint
from pathlib import Path
from cobra.core import Gene, Reaction
from cobra.util.solver import set_objective
from six import iteritems
from IPython.display import display

# Set configurations
np.set_printoptions(precision=4)

def convert_to_irreversible(cobra_model, suffixes = ('forward', 'reverse')):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    cobra_model: A Model object which will be modified in place.
    suffixes: A list of the two suffixes to use to rename the splitted reactions. Forward reaction suffix first followed by Reverse reaction suffix. 
    """
 
    reactions_to_add = []
    coefficients = {}
    for reaction in cobra_model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0:
            r_name = str(reaction.id)
            reaction.id = r_name + "_{}".format(suffixes[0])
            reverse_reaction = Reaction(r_name + "_{}".format(suffixes[1]))
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                            for k, v in iteritems(reaction._metabolites)}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    cobra_model.add_reactions(reactions_to_add)
    set_objective(cobra_model, coefficients, additive=True)
    return cobra_model

def revert_to_reversible(cobra_model):
    """This function will convert an irreversible model made by
    convert_to_irreversible into a reversible model.

    cobra_model : cobra.Model
        A model which will be modified in place.
    
    """
    
    reverse_reactions = [x for x in cobra_model.reactions
                         if "reflection" in x.notes and
                         x.id.endswith('_reverse')]

    # If there are no reverse reactions, then there is nothing to do
    if len(reverse_reactions) == 0:
        return

    for reverse in reverse_reactions:
        forward_id = reverse.notes.pop("reflection")
        forward = cobra_model.reactions.get_by_id(forward_id)
        forward.lower_bound = -reverse.upper_bound
        if forward.upper_bound == 0:
            forward.upper_bound = -reverse.lower_bound

        if "reflection" in forward.notes:
            forward.notes.pop("reflection")

    # Since the metabolites and genes are all still in
    # use we can do this faster removal step.  We can
    # probably speed things up here.
    cobra_model.remove_reactions(reverse_reactions)
    return cobra_model


class Annotate_Model(object):
    """docstring for Annotate_Model."""
    def __init__(self, dict_paths, suffixes = ['f', 'b'], modelV = 0,
    dateVersion= 'AlterOriginal', optionModel = False ,  optionIrreversible = False, additionalKcats = False, optionDefaultKcat = 'mean', threshold= 250000  ):
        
        super(Annotate_Model, self).__init__()
        self.dict_paths = dict_paths
        self.suffixes, self.suffixes_dashed, self.suffixes_PAM = self.get_suffixes(suffixes)
        self.modelV = modelV
        self.dateVersion = dateVersion
        self.optionDefaultKcat = optionDefaultKcat
        self.threshold = threshold
        self.rxn_no_record = []
        self.model_biggiML, self.irreversible_Model = self.loadModels(optionIrreversible, optionModel)
        print(f'Annotating Model {self.irreversible_Model.id}')
        #Preparing the data for the annotations
        self.dfmap_iMLsMOMENT, self.df_alterData, self.df_kcatData  = self.get_kcat_dataframes()
        self.df_kcatData = self.choose_kcat(additionalKcats)
        self.update_reaction_bounds_for_kcat()
        #Adding the annotations
        self.irreversible_Model = self.update_reaction_annotations()
        self.rxnsDatacsv = self.get_excelFile()
        self.confirm_model()


    def get_suffixes(self, suffixes):
        suffixes_dashed = ['_'+suffix for suffix in suffixes]
        suffixes_PAM = ('_f', '_b')
        return suffixes, suffixes_dashed, suffixes_PAM

    def check_model_info(self, model):
        print('Model: {}'.format(model.id))
        print('%i reactions' % len(model.reactions))
        print('%i metabolites' % len(model.metabolites))
        print('%i genes' % len(model.genes))
        
    def check_reaction(self, reaction):
        print('Reaction ID: {}'.format(reaction.id))
        print('Gene rule: {}'.format(reaction.gene_reaction_rule))
        print('Reaction Subsystem: {}'.format(reaction.subsystem))
        print('Reaction Annotations: {}'.format(reaction.annotation))

    def loadModels(self,optionIrreversible, optionModel ):
        # Load the BiGG model
        #   - Have all the information that we will include.
        model_biggiML = cameo.models.bigg.iML1515
        self.check_model_info(model_biggiML)
        # Convert the Model to irreversible (Optional if needed)
        if optionModel != False:
            irreversible_Model = optionModel
            irreversible_Model.id = 'iML1515_irreversible'
            self.check_model_info(irreversible_Model)
        elif optionIrreversible == True:
            irreversible_Model = convert_to_irreversible(model_biggiML)
            print('Model was converted to irreversible')
            self.check_model_info(model_biggiML)
        elif optionIrreversible == 'alter': 
            # Load the irreversible model from Alter Matlab format
            irreversible_Model = cobra.io.load_matlab_model(self.dict_paths['model_path'])
            irreversible_Model.id = 'iML1515_irreversible'
            self.check_model_info(irreversible_Model)
        else:
            # Load the irreversible model from a SBML format
            irreversible_Model = cobra.io.read_sbml_model(self.dict_paths['model_path'])
            irreversible_Model.id = 'iML1515_irreversible'
            self.check_model_info(irreversible_Model)
        return model_biggiML, irreversible_Model

    def get_kcat_dataframes(self):
        """
        Retrieves kcat dataframes from sMOMENT mapping and Alter et al. data.

        Parameters:
        - path_sMOMENT_reactions_kcat_map (str): Path to the sMOMENT kcat reactions mapping.
        - alter_data_path (str): Path to the Alter et al. data.

        Returns:
        - dfmap_iML1515s (pd.DataFrame): DataFrame containing sMOMENT kcat reactions mapping.
        - df_alterData (pd.DataFrame): DataFrame containing Alter et al. kcat data.
        - df_kcatData (pd.DataFrame): DataFrame for PAM data with forward and reverse columns.
        """

        # Load sMOMENT kcat reactions mapping
        dfmap_iMLsMOMENT = pd.read_json(self.dict_paths['sMOMENT_reactions_kcat_map'], orient='index')
        # display(dfmap_iMLsMOMENT.head(5))
        # print('Length dfmap: {}'.format(len(dfmap_iMLsMOMENT)))

        # Load Alter et al. data
        df_alterData = pd.read_csv(self.dict_paths['alter_data'], index_col='Reaction ID')
        names = {'Enzyme molar mass [g/mol]': 'E_mw', 'Turnover number (kcat) [1/s]': 'kcat', 'EC number': 'ec-code'}
        df_alterData.rename(columns=names, inplace=True)
        display(df_alterData.head(5))

        # Extract reaction IDs from PAM data
        reacs_ids_PAMdata_nosuffix = [id for id in df_alterData.index if id[len(self.suffixes_PAM[0])*-1:] not in self.suffixes_PAM]
        reacs_ids_PAMdata_wtsuffix = [id for id in df_alterData.index if id[len(self.suffixes_PAM[0])*-1:] in self.suffixes_PAM]
        reacs_ids_PAMdata_wosuffix = [id[:len(self.suffixes_PAM[0])*-1] for id in df_alterData.index if id[len(self.suffixes_PAM[0])*-1:] in self.suffixes_PAM]
        # Check if the lengths are consistent
        if len(df_alterData.index) == len(reacs_ids_PAMdata_nosuffix) + len(reacs_ids_PAMdata_wtsuffix):
            print('Same length: {} == {} + {}'.format(len(df_alterData.index), len(reacs_ids_PAMdata_nosuffix),
            len(reacs_ids_PAMdata_wtsuffix)))
        else:
            print('Not the same length: {} == {} + {}'.format(len(df_alterData.index), len(reacs_ids_PAMdata_nosuffix),
            len(reacs_ids_PAMdata_wtsuffix)))
        # Create a list of reactions ids with no suffix to make a consolidated list
        consolidated_list_reactionids_PAM = list(set(reacs_ids_PAMdata_nosuffix + reacs_ids_PAMdata_wosuffix))
        # print('Length consolidated list of ids in PAM data: {}'.format(len(consolidated_list_reactionids_PAM)))

        # Extract reaction IDs from model
        reacs_ids_model_nosuffix = [reaction.id for reaction in self.irreversible_Model.reactions if reaction.id[len(self.suffixes_dashed[0])*-1:] not in self.suffixes_dashed]
        reacs_ids_model_wtsuffix = [reaction.id for reaction in self.irreversible_Model.reactions if reaction.id[len(self.suffixes_dashed[0])*-1:] in self.suffixes_dashed]
        reacs_ids_model_wosuffix = [reaction.id[:len(self.suffixes_dashed[0])*-1] for reaction in self.irreversible_Model.reactions if reaction.id[len(self.suffixes_dashed[0])*-1:] in self.suffixes] 
        # Check if the lengths are consistent
        if len(self.irreversible_Model.reactions) == len(reacs_ids_model_nosuffix) + len(reacs_ids_model_wtsuffix):
            print('Same length: {} == {} + {}'.format(len(self.irreversible_Model.reactions), len(reacs_ids_model_nosuffix),
            len(reacs_ids_model_wtsuffix)))
        else:
            print('Not the same length: {} == {} + {}'.format(len(self.irreversible_Model.reactions), len(reacs_ids_model_nosuffix),
            len(reacs_ids_model_wtsuffix)))
        # Create a list of reactions ids with no suffix to make a consolidated list
        consolidated_list_reactionids_model = list(set(reacs_ids_model_nosuffix + reacs_ids_model_wosuffix))
        # print('Length consolidated list of reaction ids in irreversible model: {}'.format(len(consolidated_list_reactionids_model)))
        
        # Mergue both consolidated lists to make a final consolidated list. 
        final_consolidated_list = list(set(consolidated_list_reactionids_model)|set(consolidated_list_reactionids_PAM))

# CHECKED: logic seems correct based on the output of test
        # Create a new dataframe with a consolidated list of ids
        df_kcatData = pd.DataFrame(index=final_consolidated_list, columns=["forward", "reverse"])
        # Fill df_kcatData with alter kcats as default, can be updated later. 
        reacs_ids = df_alterData.index
        for id in reacs_ids:
            if id[len(self.suffixes_PAM[0])*-1:] in self.suffixes_PAM:
                suffix = id[len(self.suffixes_PAM[0])*-1:]
                if suffix == self.suffixes_PAM[0]:
                    df_kcatData.loc[id[:len(self.suffixes_dashed[0])*-1], 'forward'] = df_alterData.loc[id, 'kcat']

                else:
                    df_kcatData.loc[id[:len(self.suffixes_dashed[1])*-1], 'reverse'] = df_alterData.loc[id, 'kcat']
            else:
                df_kcatData.loc[id, 'forward'] = df_alterData.loc[id, 'kcat']
        test = 0
        if test==1:
            # List of reaction IDs to check
            reactions_list = ['2DGULRx', 'SPODM', 'SUCTARTtpp_f', 'GALM2pp']

            # Check if the reaction IDs are available in df_kcatData
            available_reactions = [reaction_id for reaction_id in reactions_list if reaction_id in df_kcatData.index]

            # Print out the corresponding rows
            if available_reactions:
                display(df_kcatData.loc[available_reactions, :])
            else:
                print("No matching reactions found in df_kcatData.")

        return dfmap_iMLsMOMENT, df_alterData, df_kcatData

    def choose_kcat(self, additionalKcats = False):
        """
        Choose kcat values based on the specified model version.

        Parameters:
        - df_kcatData (pd.DataFrame): DataFrame to be updated.
        - dfmap_iMLsMOMENT (pd.DataFrame): DataFrame containing sMOMENT kcat data.
        - modelV (int): Model version (default is 3).

        Returns:
        - df_kcatData (pd.DataFrame): Updated DataFrame.
        """

        reacs_ids = self.dfmap_iMLsMOMENT.dropna(subset=['forward']).index

        if self.modelV == 3:
            print(f'Model with combined kcats selecting the higher kcat')
            for id in reacs_ids:
                prev_forwardValue, prev_reverseValue = self.df_kcatData.loc[id, ['forward', 'reverse']]
                new_forwardValue, new_reverseValue = self.dmap_iMLsMOMENT.loc[id, ['forward', 'reverse']]
                if new_forwardValue > prev_forwardValue or np.isnan(prev_forwardValue):
                    self.df_kcatData.loc[id, 'forward'] = new_forwardValue
                if new_reverseValue > prev_reverseValue or np.isnan(prev_reverseValue):
                    self.df_kcatData.loc[id, 'reverse'] = new_reverseValue
        elif self.modelV == 2:
            print(f'Model with combined kcats priority sMOMENT kcats')
            for id in reacs_ids:
                new_forwardValue, new_reverseValue = self.dfmap_iMLsMOMENT.loc[id, ['forward', 'reverse']]
                self.df_kcatData.loc[id, 'forward'] = new_forwardValue
                self.df_kcatData.loc[id, 'reverse'] = new_reverseValue
        elif self.modelV == 1:
            print(f'Model with combined kcats getting the mean of both datasets')
            for id in reacs_ids:
                new_forwardValue, new_reverseValue = self.dfmap_iMLsMOMENT.loc[id, ['forward', 'reverse']]
                old_forwardValue, old_reverseValue = self.df_kcatData.loc[id, ['forward', 'reverse']]
                self.df_kcatData.loc[id, 'forward'] = (new_forwardValue + old_forwardValue)/2
                self.df_kcatData.loc[id, 'reverse'] = (new_reverseValue + old_reverseValue )/2
        else:
            print(f'Model only using the Alter kcats')

        if additionalKcats != False:
            reacs_ids = additionalKcats.index
            for id in reacs_ids:
                if id[len(self.suffixes[0])*-1:] in self.suffixes:
                    suffix = id[len(self.suffixes[0])*-1:]
                    if suffix == self.suffixes[0]:
                        self.df_kcatData.loc[id[:len(self.suffixes[0])*-1], 'forward'] = additionalKcats.loc[id, 'kcat']
                    else:
                        self.df_kcatData.loc[id[:len(self.suffixes[1])*-1], 'reverse'] = additionalKcats.loc[id, 'kcat']
                else:
                    self.df_kcatData.loc[id, 'forward'] = additionalKcats.loc[id, 'kcat']
        # display(self.df_kcatData.loc[reacs_ids,:])
# CHECKED: Aad a test to check output for plobematic reactions
        test = 0
        if test :
            print('\nTest for choose_kcat : \n')
            
            # List of reaction IDs to check
            reactions_list = ['2DGULRx', 'SPODM', 'SUCTARTtpp_f', 'GALM2pp']

            # Check if the reaction IDs are available in df_kcatData
            available_reactions = [reaction_id for reaction_id in reactions_list if reaction_id in self.df_kcatData.index]

            # Print out the corresponding rows
            if available_reactions:
                display(self.df_kcatData.loc[available_reactions, :])
            else:
                print("No matching reactions found in df_kcatData.")

            raise ValueError("Testrun finnished. Execution halted.")

        return self.df_kcatData

    def update_reaction_bounds_for_kcat(self):
        """
        Update the model reaction bounds based on kcat data.

        Parameters:
        - model (cobra.Model): Constraint-Based Reconstruction and Analysis (COBRA) model.
        - df_kcatData (pd.DataFrame): DataFrame containing kcat data.

        Returns:
        - rxnsNoReverse_bound (dict): Dictionary of reaction IDs with updated lower bounds.
        - rxnsNoForward_bound (dict): Dictionary of reaction IDs with updated upper bounds.
        """

        # Getting a list of reactions with a reverse kcat
        list_reverse_reactions_ids = self.df_kcatData.loc[~self.df_kcatData.isnull()['reverse']].index.values
        # print(list_reverse_reactions_ids)

        # Update the reactions boundaries
        rxnsNoReverse_bound = {}
        rxnsNoForward_bound = {}

        for reaction_id in list_reverse_reactions_ids:
            if reaction_id in self.irreversible_Model.reactions:
                reaction = self.irreversible_Model.reactions.get_by_id(reaction_id)
            else: 
                reaction_id_suffix = reaction_id + self.suffixes_dashed[1]
                if reaction_id_suffix in self.irreversible_Model.reactions:
                    reaction = self.irreversible_Model.reactions.get_by_id(reaction_id_suffix)
                else:
                    print(f'Reaction {reaction_id} not included in the Model {self.irreversible_Model.id}')
                    continue
            # Update lower bound for reversible reactions
            if reaction.lower_bound > -1000 and reaction.lower_bound < 0 :
                print(f'The reaction {reaction.id} has a lower bound of {reaction.lower_bound}')
                rxnsNoReverse_bound[reaction.id] = reaction.lower_bound
                reaction.lower_bound = 0

            # Update upper bound for reversible reactions
            if reaction.upper_bound < 1000 :
                print(f'The reaction {reaction.id} has an upper bound of {reaction.upper_bound}')
                rxnsNoForward_bound[reaction.id] = reaction.upper_bound
                reaction.upper_bound = 1000
        # Display the reaction bounds that were updated
        print("Reactions with updated lower bounds:", rxnsNoReverse_bound)
        print("Reactions with updated upper bounds:", rxnsNoForward_bound)
        
        return rxnsNoReverse_bound, rxnsNoForward_bound

    def add_stepdb_annotations(self):
            """
            Add STEPdb annotations to the genes in the given CobraPy model.

            Parameters:
            - model (cobra.Model): The CobraPy metabolic model.
            - stepdb_path (str): Path to the STEPdb file.

            Returns:
            - model (cobra.Model): The modified CobraPy metabolic model.
            """

            # Load STEPdb
            step_database = pd.read_csv(self.dict_paths['STEPdb_data'], index_col='Gene_ID')
            display(step_database)
            step_database = step_database.fillna('nan')

            geneLocalizations = [x for x in step_database['Gene_Location']]

            # Adding GENE annotations to the model from STEPdb
            for ind, gene in enumerate(self.irreversible_Model.genes):
                gene_id = gene.id

                # Extract information from STEPdb and add annotations
                gene.annotation['stepdb.subsystem'] = step_database.loc[gene_id, 'Gene_Subsys']
                gene.annotation['stepdb.location'] = step_database.loc[gene_id, 'Gene_Location']
                gene.annotation['stepdb.weight'] = str(step_database.loc[gene_id, 'Molecular_Weight'])
                gene.annotation['stepdb.name'] = step_database.loc[gene_id, 'Gene_Names_x']
                gene.annotation['stepdb.tmLength'] = step_database.loc[gene_id, 'Transmembrane_Length']
                gene.annotation['stepdb.tmCount'] = str(step_database.loc[gene_id, 'Count_Transmembrane'])
                gene.annotation['stepdb.chainLength'] = step_database.loc[gene_id, 'Chain_Length']
                gene.annotation['stepdb.tmPerc'] = str(step_database.loc[gene_id, 'Transmembrane_Percentage'])
                gene.annotation['stepdb.tmWeight'] = str(step_database.loc[gene_id, 'Transmembrane_MolWeight'])
                # gene.annotation['stepdb.tmSA'] = str(step_database.loc[gene_id, 'Surface_Area_Calculated'])

                # Uncomment the line below if you want to replace the gene in the model
                # model.genes[ind] = gene

            # Display information about a specific gene in the model
            example_gene = self.irreversible_Model.genes[45]
            print(example_gene.annotation)
            return [self.irreversible_Model , step_database, geneLocalizations]

    def parse_geneRule(self, reaction):
        ''' Function to parse the gene rules of a metabolic model to calculate the w parameter'''
        rule = reaction.gene_reaction_rule
        or_rule = list(map(lambda x: x.strip(' ()'), rule.split('or')))
        parsed_rule = []
        for group in or_rule:
            and_rule = list(map(lambda x: x.strip(' ()'), group.split('and')))
            parsed_rule.append(and_rule)
        return parsed_rule

    def get_geneInfoLocMw(self, geneid):
        if geneid in self.step_database.index :
            gene_location = self.step_database.loc[geneid, 'Gene_Location']
            gene_mW = float(self.step_database.loc[geneid, 'Molecular_Weight'])
            if gene_mW < 0 or gene_mW == np.nan:
                gene_mW = 0
            gene_mWMb = float(self.step_database.loc[geneid, 'Transmembrane_MolWeight'])
            if gene_mWMb < 0 or gene_mWMb == np.nan:
                gene_mWMb = 0
            counts_tmb_gene = self.step_database.loc[geneid, 'Count_Transmembrane']    
        else:
            gene_location = ''
            gene_mW = 0
            gene_mWMb = 0
            counts_tmb_gene = 0
        if gene_mWMb != 0 and gene_location == 'Cytoplasmic':
            print(f'check gene {geneid}')
        return gene_location, gene_mW, gene_mWMb,counts_tmb_gene
        
    def get_reactionLocMw(self, reaction):
        # Get the gene rules for the reaction
        gene_rules = self.parse_geneRule(reaction)
        
        # Initialize lists to store location and molecular weight information
        locations = []
        mw_totals = []
        mw_Mb_totals = []
        counts_tmb_totals = []

        # Loop through each set of OR conditions in the gene rules
        for or_conditions in gene_rules:
            # Initialize lists to store information for AND conditions
            and_locations = []
            and_mw = []
            and_mw_Mb = []
            and_counts_tmb= []

            # Loop through each AND condition
            for gene_id in or_conditions:
                # Get the location and molecular weight for the gene
                gene_location, gene_mW, gene_mWMb,counts_tmb_gene = self.get_geneInfoLocMw(gene_id)  
                # Append location and molecular weight for the AND condition
                and_locations.append(gene_location)
                and_mw.append(gene_mW)
                # Append the molecular weight for transmembrane
                and_mw_Mb.append(gene_mWMb)
                and_counts_tmb.append(counts_tmb_gene)

            # For 'Inner Membrane', sum the molecular weights for transmembrane
            mw_Mb_totals.append(sum(and_mw_Mb))
            # Append total molecular weight for the OR condition
            mw_totals.append(sum(and_mw))
            # Append total counts of transmembrane domains
            counts_tmb_totals.append(sum(and_counts_tmb))
            # Determine the location for the OR condition
            check_membraneLoc = [location for location in and_locations if 'Inner Membrane' in location ]
            if check_membraneLoc: 
                locations.append('Inner Membrane')
            elif reaction.subsystem == 'Extracellular exchange':
                locations.append('Exchange')
            else:
                locations.append('Cytoplasmic')

        # Find the minimum total molecular weight
        mw_min = self.find_minimum_value(mw_totals,gene_rules)
        index_mw_min = mw_totals.index(mw_min)
        # Find the minimum total molecular weight for transmembrane 
        mw_Mb_min = self.find_minimum_value(mw_Mb_totals,gene_rules)
        index_mw_Mb_min = mw_Mb_totals.index(mw_Mb_min)
        counts_tmb_min = counts_tmb_totals[index_mw_Mb_min] 

        #Get crossed min for testing
        crossed_mw = mw_totals[index_mw_Mb_min]
        crossed_mw_Mb = mw_Mb_totals[index_mw_min]
 
        if index_mw_min != index_mw_Mb_min :
            # print(f'For reaction {reaction.id} the mW and mWMb minimums are not the same.')
            # print(f'Crossed minimums mW -> {crossed_mw}, mWMb -> {crossed_mw_Mb}')
            # Selecting the minimum of the total mW to select its Membrane weigth
            mw_Mb_min = crossed_mw_Mb #
            # Selecting the minimum of the total counts for transmembrane domains
            counts_tmb_min = counts_tmb_totals[index_mw_min] #

        # Get the corresponding location
        location_min = locations[index_mw_min]
        if mw_Mb_min == 0 and location_min == 'Inner Membrane':
            print(f'check reaction {reaction.id}')
        return location_min, mw_min, mw_Mb_min, counts_tmb_min, mw_totals, mw_Mb_totals, locations, counts_tmb_totals
    
    def find_minimum_value(self, weight_list,gene_rules):
        # Check if the list contains only zeros
        flattened_grules = [item for sublist in gene_rules for item in sublist]
        if 's0001' in flattened_grules:
            return 0
        elif all(weight == 0 for weight in weight_list):
            return 0
        else:
            # Filter out zero values and find the minimum of the remaining values
            non_zero_values = [weight for weight in weight_list if weight != 0]
            # If the list has only zeros, return 0; otherwise, return the minimum non-zero value
            return min(non_zero_values, default=0)

    def get_Defaultkcat(self):
        # Obtain the lists with the locations and kcats for each reaction
        reactions_locations = [x.annotation['location'] for x in self.irreversible_Model.reactions]  
        reactions_kcats = [float(x.annotation['kcat']) for x in self.irreversible_Model.reactions]
        # Obtain the unique locations to filter the data
        unique_locations = set(reactions_locations)
        pre_df = pd.DataFrame({'Locations':reactions_locations, 'Kcats': reactions_kcats }) 
        # Initialize a dictionary to store default values per location
        dict_defaultKcat = {}
        # Obtain the statistic option to be used
        statOption = self.optionDefaultKcat
        for location in unique_locations:
            # Filter the dataframe based on location and kcat below a threshold
            estimation_df = pre_df.loc[(pre_df['Locations'] == location) & (pre_df['Kcats'] <= self.threshold) & (pre_df['Kcats'] > 0)].copy()
            dict_defaultKcat[location] = estimation_df['Kcats'].apply(statOption)
            
        return dict_defaultKcat

    def correct_kcat(self):
        for ind, reaction in enumerate(self.irreversible_Model.reactions):
            kcat = float(reaction.annotation['kcat']) 
            # CHECKED: The issue is the threshold used since its value is less than some elevated kcats. 
            if kcat == None or kcat <= 0 or kcat > 250000 or kcat == np.nan: 
                self.rxn_no_record.append(reaction.id)
                location = reaction.annotation['location']
                if location == 'Exchange' or 'diffusion' in reaction.name:
                    kcat = 100000000
                else:
                    # kcat = self.dict_defaultKcat[location]
                    kcat = 100000000
                print(f'Reaction {reaction.id} kcat was set to {kcat:.4f}')
            reaction.annotation['kcat'] = str(kcat)
            self.irreversible_Model. reactions[ind] = reaction
        return self.irreversible_Model

    def annotate_reactionLocMw(self):
        for ind , reaction in enumerate(self.irreversible_Model.reactions):
            location_min, mw_min, mw_Mb_min, counts_tmb_min, mw_totals, mw_Mb_totals, locations, counts_tmb_totals = self.get_reactionLocMw(reaction)
            reaction.annotation['location'] = location_min
            reaction.annotation['weight_mw_c'] = str(mw_min)
            reaction.annotation['weight_mw_m'] = str(mw_Mb_min)
            reaction.annotation['number_tMb'] = str(counts_tmb_min)
            reaction.annotation['mw_totals'] = str(mw_totals)
            reaction.annotation['mw_Mb_totals'] = str(mw_Mb_totals)
            reaction.annotation['locations_All'] = str(locations)
            reaction.annotation['counts_tmb_totals'] = str(counts_tmb_totals)
            reaction.annotation['surface_area_tMb'] = str(counts_tmb_min*1.5)
            self.irreversible_Model.reactions[ind] = reaction 
        print(self.irreversible_Model.reactions[45].annotation)
        return self.irreversible_Model


    def manual_update_model(self):
        # Load Update data
        update_path = self.dict_paths.get('update_data')
        if not update_path:
            print("No update_data sheet provided. Skipping optional manual updates.")
            return
        update_data_dict = pd.read_excel(update_path, sheet_name=None)
        
        # Mapping between Excel column names and model annotation field names
        match_fields_ids = {
            'Location': 'location',
            'Molecular_weight': 'weight_mw_c',
            'Molecular_weight_Tmb': 'weight_mw_m',
            'MSA_parameter': 'surface_area_tMb',
            'Kcat': 'kcat'
        }
        
        # Iterate over each sheet in the Excel file
        for key_id, update_data_df in update_data_dict.items():
            # Skip empty sheets
            if update_data_df.empty:
                continue
            
            # Check if the sheet contains the required columns
            if 'reaction_id' not in update_data_df.columns or 'value' not in update_data_df.columns:
                print(f"Error: Sheet '{key_id}' is missing 'reaction_id' or 'value' columns.")
                continue
            
            update_data_df.set_index('reaction_id', inplace=True)
            field_id = match_fields_ids.get(key_id)
            
            if not field_id:
                print(f"Error: No matching field ID found for sheet '{key_id}'.")
                continue

            # Iterate over each reaction in the sheet
            for reaction_id, value in update_data_df['value'].items():
                # Get the reaction from the model
                reaction = self.irreversible_Model.reactions.get_by_id(reaction_id)
                if not reaction:
                    print(f"Warning: Reaction '{reaction_id}' not found in the model.")
                    continue
                
                # Update the annotation field
                reaction.annotation[field_id] = value
        
        return self.irreversible_Model


    def update_reaction_annotations(self):
        """
        Update reaction annotations in the given CobraPy model using data from the kcat database.

        Parameters:
        - model (cobra.Model): The CobraPy metabolic model.
        - kcat_path (str): Path to the kcat database JSON file.
        - df_alterData (pd.DataFrame): DataFrame containing Alter model data.
        - df_kcatData (pd.DataFrame): DataFrame containing kcat data.
        - suffixes (list): List of suffixes.

        Returns:
        - model (cobra.Model): The modified CobraPy metabolic model.
        """

        # Load kcat database
        combined_kcat_database = json.load(open(self.dict_paths['combined_kcatdb']))

        self.irreversible_Model,self.step_database, geneLocalizations = self.add_stepdb_annotations()
        self.irreversible_Model = self.annotate_reactionLocMw()

        # Iterate over reactions in the model
        for ind, reaction in enumerate(self.irreversible_Model.reactions):
            reaction_id = reaction.id
            reaction_suffix = reaction.id[-2:]

            # Update 'step.subsystem' annotation
            reaction.annotation['step.subsystem'] = reaction.subsystem

            # Determine reaction ID without suffix
            if reaction_suffix in self.suffixes_dashed :
                reaction_id_without_suffix = reaction_id[:-2]
            else:
                reaction_id_without_suffix = reaction_id

            # Update 'ec-code' annotation
            if ("ec-code" not in reaction.annotation.keys() or reaction.annotation['ec-code'] is None) and \
                    reaction_id in self.df_alterData.index.values:
                reaction.annotation['ec-code'] = str(self.df_alterData.loc[reaction_id]['ec-code']).split(', ')

            if "ec-code" in reaction.annotation.keys() and reaction.annotation['ec-code'] is None:
                reaction.annotation['ec-code'] = ['nan']

            if "ec-code" in reaction.annotation.keys():
                ec_code_list = reaction.annotation['ec-code']

                # Update 'ec-code' with information from kcat database
                for reaction_id_entry in ec_code_list:
                    if reaction_id_entry in combined_kcat_database.keys():
                        entry = combined_kcat_database[reaction_id_entry]
                        if 'TRANSFER' in entry.keys():
                            reaction.annotation['ec-code'].append(entry['TRANSFER'])

                # Remove duplicates
                reaction.annotation['ec-code'] = list(dict.fromkeys(reaction.annotation['ec-code']))

            # Update 'kcat_alt' annotation
            if reaction_id in self.df_alterData.index.values:
                reaction.annotation['kcat_alt'] = str(self.df_alterData.loc[reaction_id]['kcat'])

            # Adding kcat values from the kcatDB
            if reaction_id_without_suffix in self.df_kcatData.index.values:
                if reaction_suffix == '_b':
                    reaction.annotation['kcat'] = str(self.df_kcatData.loc[reaction_id_without_suffix]['reverse'])
                else:
                    reaction.annotation['kcat'] = str(self.df_kcatData.loc[reaction_id_without_suffix]['forward'])
            else:
                reaction.annotation['kcat'] = str(np.nan)
                    
            # Update the reaction in the model
            self.irreversible_Model.reactions[ind] = reaction
        self.dict_defaultKcat = self.get_Defaultkcat()
        # self.irreversible_Model = self.manual_update_model()
        # TODO: UPDATE REACTIONS BASE ON MANUAL ADJUSTS
        # self.irreversible_Model = self.correct_kcat()
# REVIEW: Check logic to see if match other problematic reactions
        test = 0
        if test:
            print('\nTest for choose_kcat : \n')

            # Display information about a specific reaction in the model
            if 'ATPS4rpp_b' in self.irreversible_Model.reactions:
                example_reaction = self.irreversible_Model.reactions.ATPS4rpp_b
                print(example_reaction.annotation)

            # List of reaction IDs to check
            reactions_list = ['2DGULRx', 'SPODM', 'SUCTARTtpp_f', 'GALM2pp']

            # Check if the reaction IDs are available in df_kcatData
            modelRxnsIds = [r.id for r in self.irreversible_Model.reactions]
            available_reactions = [reaction_id for reaction_id in reactions_list if reaction_id in modelRxnsIds]

            # Print out the corresponding rows
            if available_reactions:
                for reactionId in available_reactions:
                    reaction = self.irreversible_Model.reactions.get_by_id(reactionId)
                    self.check_reaction(reaction)
            else:
                print("No matching reactions found in the model.")

            raise ValueError("Testrun finnished. Execution halted.")
        return self.irreversible_Model

    def confirm_model(self):
        """
        Confirm model based on specific annotations for reactions ATPS4rpp_b and G3PD2_b.
        Returns True if the model meets the conditions, False otherwise.
        """
        for reaction_id in ['ATPS4rpp_b', 'G3PD2_b']:
            reaction = self.irreversible_Model.reactions.get_by_id(reaction_id)

            # Check annotations for location, weight_mw_c, weight_mw_m, and kcat
            location = reaction.annotation.get('location', '')
            weight_mw_c = reaction.annotation.get('weight_mw_c', 0)
            weight_mw_m = reaction.annotation.get('weight_mw_m', 0)
            kcat = reaction.annotation.get('kcat', 0)

            # Conditions for membrane (you may need to adjust these based on your data)
            membrane_conditions = location == 'Inner Membrane' and float(weight_mw_c) > 0 and float(weight_mw_m) > 0 and float(kcat) > 0

            # Conditions for cytoplasmic (you may need to adjust these based on your data)
            cytoplasmic_conditions = location == 'Cytoplasmic' and float(weight_mw_c) > 0 and float(weight_mw_m) == 0 and float(kcat) > 0

            # Check if the conditions are met for either membrane or cytoplasmic
            if not (membrane_conditions or cytoplasmic_conditions):
                print(f"Model confirmation failed for reaction {reaction_id}.")
                return False

        print("Model confirmed.")
        return True

    def get_excelFile(self,output_path=False):
        """
        Save model-related data, including kcat data, reaction attributes, and specific reaction information.

        Parameters:
        - model (cobra.Model): The CobraPy metabolic model.
        - df_alterData (pd.DataFrame): DataFrame containing Alter model data.
        - df_kcatData (pd.DataFrame): DataFrame containing kcat data.
        - model_id (str): Identifier for the model.
        - model_version (int): Version of the model.
        - date_version (str): Date version.

        Returns:
        - None
        """
        if output_path == False:
            model_id = self.irreversible_Model.id
            if 'output_folder' in self.dict_paths.keys():
                output_folder = self.dict_paths['output_folder']
            else:
                output_folder = ''
            if output_folder:
                os.makedirs(os.path.join(output_folder, 'models_InputData'), exist_ok=True)
            output_path = f'{output_folder}models_InputData/{model_id}_kcatDataBaseV{self.modelV}_{self.dateVersion}.csv'

        # Save kcat data as CSV
        self.df_kcatData.to_csv(output_path)

        # Create a dataframe for reaction attributes
        pre_df = {
            'Rxn_ID': [],
            'Rxn_name': [],
            'Rxn_Loc': [],
            'Rxn_Subsystem': [],
            'Rxn_kcat': [],
            'Rxn_genes_rules': [],
            'Rxn_mW_c' : [],
            'Rxn_mW_m' : [],
            'Rxn_numTMb' : [], 
            'Rxn_mW_c_Total': [],
            'Rxn_mW_m_Total': [],
            'Rxn_locations_Total': [],
            'Rxn_count_Tmb_Total': []
        }
        labels_match = { 'Rxn_Loc':'location', 'Rxn_kcat':'kcat', 'Rxn_mW_c':'weight_mw_c', 'Rxn_mW_m':'weight_mw_m', 'Rxn_numTMb':'number_tMb', 'Rxn_mW_c_Total': 'mw_totals',             'Rxn_mW_m_Total': 'mw_Mb_totals', 'Rxn_locations_Total': 'locations_All', 'Rxn_count_Tmb_Total': 'counts_tmb_totals' }
    

        for reaction in self.irreversible_Model.reactions:
            pre_df['Rxn_name'].append(reaction.name)
            pre_df['Rxn_ID'].append(reaction.id)
            pre_df['Rxn_genes_rules'].append(reaction.gene_reaction_rule)
            pre_df['Rxn_Subsystem'].append(reaction.subsystem)

            for key, label in labels_match.items():
                if label in reaction.annotation.keys():
                    pre_df[key].append(reaction.annotation[label])
                else:
                    pre_df[key].append('')

        df_rxns = pd.DataFrame(pre_df)

        # Save reaction attribute data as CSV
        if output_folder:
            os.makedirs(os.path.join(output_folder, 'models_InputData'), exist_ok=True)
        df_rxns.to_csv(f'{output_folder}models_InputData/{model_id}_reactionDataV{self.modelV}_{self.dateVersion}.csv')

        return df_rxns
    
    def save_model_as_xml(self, output_path=False):
        """
        Save the annotated model as an XML file.

        Parameters:
        - output_path (str): Path to save the XML file.
        """
        if hasattr(self, 'irreversible_Model'):

            if output_path == False:
                model_id = self.irreversible_Model.id
                if 'output_folder' in self.dict_paths.keys():
                    output_folder = self.dict_paths['output_folder']
                else:
                    output_folder = ''
                if output_folder:
                    os.makedirs(os.path.join(output_folder, 'model_OutputData'), exist_ok=True)
                output_path = f'{output_folder}model_OutputData/annotated_{model_id}_V{self.modelV}_{self.dateVersion}.xml'

            cobra.io.write_sbml_model(self.irreversible_Model, output_path)
            print(f"Annotated model saved as XML file: {output_path}")
        else:
            print("No annotated model available. Please run the annotation process first.")





#%%


if __name__ == '__main__':
    base_dir = Path(__file__).resolve().parent
    repo_root = base_dir.parents[2]
    data_input = repo_root / "data" / "input"
    stepdb_dir = data_input / "STEPdb"

    suffixes_PAM = ['f', 'b']
    # Model and version information
    modelV = 1
    dateVersion = 'MAFBA'
    model_id = 'iML1515'
    model_code = f"{model_id}_{dateVersion}"
    output_root = repo_root / "1_Model_Construction" / "03_MAFBA_Models_Outputs" / model_code
    (output_root / "models_InputData").mkdir(parents=True, exist_ok=True)
    (output_root / "model_OutputData").mkdir(parents=True, exist_ok=True)
    model_path = data_input / f"mafba_{model_id}_ecV{modelV}_g1_{dateVersion}.xml"
    update_data = None  # optional; set this to an update Excel path if needed

    dict_paths = {
        'model_path': str(model_path),
        'sMOMENT_reactions_kcat_map': str(data_input / "iML1515_star_reactions_kcat_mapping_combined.json"),
        'alter_data': str(data_input / "model_construction" / "PAM_alter_data.csv"),
        'STEPdb_data': str(stepdb_dir / "STEPdb_iML1515_prepared_corrected.csv"),
        'combined_kcatdb': str(data_input / "combined_brenda_sabio_rk_iML_irr.json"),
        'input_folder': str(data_input) + '/',
        'output_folder': str(output_root) + '/',
        'update_data': update_data
    }
    optionIrreversible = 'alter'
    additionalKcats = False
    optionDefaultKcat = 'mean'

    # Load the irreversible model from Alter's paper
    model_iML1515_alter = cobra.io.load_matlab_model(dict_paths['model_path'])
    model_iML1515_alter.id = 'iML1515'

    annotatedModel = Annotate_Model(dict_paths, suffixes_PAM, modelV, dateVersion,  model_iML1515_alter,  optionIrreversible, additionalKcats, optionDefaultKcat)

    print(annotatedModel.irreversible_Model.reactions[117].annotation)


# %%

# annotatedModel.save_model_as_xml()
    

# %%
