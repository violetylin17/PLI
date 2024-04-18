import pandas as pd
import numpy as np
def cluster_list(component_df):
    # Input: component df
    # Output: List, drug-cluster distribution
    cluster = list(component_df.idxmax()) # Extract maximum of every col
    
    return cluster

def find_profile(profile_matrix):
    # Input: profile matrix
    # Output: find non-zero term of every residue
    feature_number = len(profile_matrix.columns)
    output = {}
    
    for i in range(feature_number):
        output[i] = {j:profile_matrix[i].to_dict()[j] for j in profile_matrix[i].to_dict()\
                     if profile_matrix[i].to_dict()[j] > 0}
    return output

def matrix_orth_check(dataframe):
    # Input: dataframe
    # Output: np.matrix, matrix * matrix.T
    matrix = np.matrix(dataframe)
    
    output = matrix * matrix.T
    return output


def ScoreAnalysis(store):
    #Input : Score.h5
    #Output: list of dataframe, st hb sthb hy total
    
    ratom = store["receptor/atoms"]
    # Drug-Residue interaction
    
    # Create receptor atom id to residue name mapping
    ratom["resname"] = ratom["segid"].apply(str) + "/" + ratom["resid"].apply(str)
    t1 = ratom.set_index("id") # pd of receptor atom but use id as index of the dataframe 
    radict = {index:atom_info['resname'] for (index,atom_info) in t1.to_dict('index').items()}
    
    # Create residue name list
    t2 = list(radict.values()) 
    t3 = np.unique(t2, return_index=True)[1] # atom id
    ralist = [t2[index] for index in sorted(t3)]
    

    
    # Create receptor atom id to residue name mapping
    ratom["resname"] = ratom["segid"].apply(str) + "/" + ratom["resid"].apply(str) + "/" + ratom["resn"].apply(str)
    t1 = ratom.set_index("id") # pd of receptor atom but use id as index of the dataframe 
    radict = {index:atom_info['resname'] for (index,atom_info) in t1.to_dict('index').items()}
    
    # Create residue name list
    t2 = list(radict.values()) 
    t3 = np.unique(t2, return_index=True)[1] # atom id
    ralist = [t2[index] for index in sorted(t3)]
    
    
    DRinteraction = pd.DataFrame()
    temp_dri = []
    temp_dri_st = []
    temp_dri_hb = []
    temp_dri_sthb=[]
    temp_dri_hy=[]
    temp_dri_total=[]
    for ligand in store.root.ligands:
        ligand_name_mode = str(ligand).split(' ')[0].split('/')[2]
    
        ligand_name = ligand_name_mode.split(',')[0]
        ligand_mode = ligand_name_mode.split(',')[1]
    
        ### Only remain first 1
        if ligand_mode != '1':
            continue
        else:
            pass
        ###
    
        interaction_df = store[str(ligand.interaction).split(" ")[0]]
    
        # Interaction per receptor atom
        ratom_interaction = interaction_df.groupby('rid', as_index=False).sum()[['rid','st','hy','hb']]
        ratom_interaction['resname'] = ratom_interaction['rid'].map(radict)
    
        # Sum the interaction if map to same resname
        residue_interaction = ratom_interaction.groupby('resname', as_index=False).sum()[['resname','st','hy','hb']]
        residue_interaction.set_index('resname', inplace=True) # Inplace True: dont create new object
    
        # Fill residue-interaction table with zero
        residue_interaction = residue_interaction.reindex(ralist, fill_value=0.0) # Fill zero to no-value residue
    
        # Define sthy = st+hy 
        residue_interaction['sthb'] = residue_interaction['st'] + residue_interaction['hb']
        residue_interaction['total'] = residue_interaction['sthb'] + residue_interaction['hy']
        residue_interaction = residue_interaction[['st','hb','sthb','hy','total']]
        
        residue_interaction_st = residue_interaction[['st']]
        residue_interaction_hb = residue_interaction[['hb']]
        residue_interaction_sthb = residue_interaction[['sthb']]   #sthy dataframe 
        residue_interaction_hy = residue_interaction[['hy']]       #hb dataframe
        residue_interaction_total = residue_interaction[['total']] #total dataframe
        # Create multi-level index
        stack_table = residue_interaction.stack()
        stack_table.index.set_names(['resname','type'], inplace=True) 
    
    
        residue_interaction = pd.DataFrame(stack_table)
        residue_interaction.columns = [f'{ligand_name},{ligand_mode}']  
        residue_interaction_st.columns = [f'{ligand_name},{ligand_mode}'] 
        residue_interaction_hb.columns = [f'{ligand_name},{ligand_mode}'] 
        residue_interaction_sthb.columns = [f'{ligand_name},{ligand_mode}'] 
        residue_interaction_hy.columns = [f'{ligand_name},{ligand_mode}'] 
        residue_interaction_total.columns = [f'{ligand_name},{ligand_mode}']
    
    
        temp_dri.append(residue_interaction)
        temp_dri_st.append(residue_interaction_st)
        temp_dri_hb.append(residue_interaction_hb)
        temp_dri_sthb.append(residue_interaction_sthb)
        temp_dri_hy.append(residue_interaction_hy)
        temp_dri_total.append(residue_interaction_total)
    
    DRinteraction = pd.concat(temp_dri, axis=1)
    DRinteraction_st = pd.concat(temp_dri_st,axis=1)
    DRinteraction_hb = pd.concat(temp_dri_hb,axis=1)
    DRinteraction_sthb = pd.concat(temp_dri_sthb,axis=1)
    DRinteraction_hy = pd.concat(temp_dri_hy,axis=1)
    DRinteraction_total = pd.concat(temp_dri_total, axis=1)
    return [DRinteraction, DRinteraction_st, DRinteraction_hb, DRinteraction_sthb, \
            DRinteraction_hy, DRinteraction_total]

def Union_set(interaction_pd):
    # Input: drug-residue score dataframe 
    # output: residue union set 
    residue_list = [i for i in interaction_pd.index]
    non_zero_interaction = []
    for i in interaction_pd.columns :
        single_drug = interaction_pd[i]
        temp = {residue : float(single_drug[residue]) for residue in residue_list if single_drug[residue] != 0}
        non_zero_interaction.append(temp)
    
    Union_set = set()
    for drug in non_zero_interaction:
        Union_set = Union_set | set(drug.keys())
    
    return Union_set

def significant_residue_set(interaction_pd):
    # Input : drug-residue score dataframe (dtype -> pd.df )
    # Output: which residue get min energy through all drug (dtype->set)
    return set(pd.DataFrame.idxmin(interaction_pd).values)

def neg_extractor(interaction_pd):
    return interaction_pd.where(interaction_pd<0.0, 0)

def pos_extractor(interaction_pd):
    return interaction_pd.where(interaction_pd>0.0, 0)

def ratio_df(interaction_pd):
    #Input :drug-residue score dataframe
    #Output:the ratio energy of every col (ex: -1;-2;-3 -> 1/6 1/3 1/2)
    # -3/-6 = 1/2 do not need to change sign
    for i in interaction_pd.columns:
        interaction_pd[i] = interaction_pd[i] / interaction_pd.sum(0)[i]
    
    return interaction_pd


def torsion_num(mol_REMARK):
    # Input  :the object of rdkit.Chem.PandasTools.LoadSDF, single mol.REMARK
    # Output :the number of active torsion of mol
    
    active_tor_line = mol_REMARK.split('\n')[2]
    active_tor = int(active_tor_line.split()[0])
    
    return active_tor