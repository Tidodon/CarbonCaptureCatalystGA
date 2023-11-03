from rdkit import Chem
#from rdkit.Chem import AllChem


#### Fix to amine to zwitterion TS (TS1).

patt_prim_amine, patt_seco_amine =Chem.MolFromSmarts("[ND1]") ,Chem.MolFromSmarts("[ND2]") 
##Replacement chunks and recognitions for zwitterion formation.
repl_zwitterion_prim = Chem.MolFromSmiles("[NH2+]C(=O)[O-]")  # NCCO -> [O-](=O)C-[NH2+]CCO
repl_zwitterion_seco= Chem.MolFromSmiles("[NH+]C(=O)[O-]")# CNCC ->

##Replacement chunks and recognitions for Carbamate formation.
repl_carbamate_prim= Chem.MolFromSmiles("[NH]C(=O)[O-]")# NCCO -> [O-](=O)C-[NH2+]CCO
repl_carbamate_seco = Chem.MolFromSmiles("NC(=O)[O-]") # CNCC ->

##Replacement chunks and recognitions for Carbamate pair product formation.
repl_carba_pair_prim = Chem.MolFromSmiles("[H][N+]([H])([H])C(=O)[O-]")
repl_carba_pair_seco = Chem.MolFromSmiles("[H][N+]([H])C(=O)[O-]")




def compute_dG(mol_smile, options=[]):
    pass

def prepare_smiles(mol_smile) -> tuple:
    """
    outputs tuple, with elements(by order): zwitterion, carbamate, and the ionized
    """

    mol = Chem.MolFromSmiles(mol_smile)
    zwitterions_p = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt_prim_amine, replacement=repl_zwitterion_prim) ## primary amine substitutions
    zwitterions_s = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt_seco_amine, replacement=repl_zwitterion_seco) ## seconadary amine substitutions
    zwitterions = zwitterions_p+zwitterions_s ### Combine products
    zwitterions = set([Chem.MolToSmiles(p) for p in zwitterions]) ###Remove duplicates


    ##Similar for carbamates.
    carbamates_p = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt_prim_amine, replacement=repl_carbamate_prim) ## primary amine substitutions
    carbamates_s = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt_seco_amine, replacement=repl_carbamate_seco) ## seconadary amine substitutions
    carbamates = carbamates_p+carbamates_s ### Combine products
    carbamates = set([Chem.MolToSmiles(p) for p in carbamates]) ###Remove duplicates

    ## Carbamate pairs;

    carba_pairs_p = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt_prim_amine, replacement=repl_carba_pair_prim) ## primary amine substitutions
    carba_pairs_s = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt_seco_amine, replacement=repl_carba_pair_seco) ## seconadary amine substitutions
    carba_pairs   = carba_pairs_p+carba_pairs_s ### Combine products
    carba_pairs   = set([Chem.MolToSmiles(p) for p in carba_pairs]) ###Remove duplicates

    carbamates.remove(mol_smile)
    carba_pairs.remove(mol_smile)
    zwitterions.remove(mol_smile)

    outs = []

    carbamates_mols  = [v for v in map(Chem.MolFromSmiles, carbamates )]  
    zwitterions_mols = [v for v in map(Chem.MolFromSmiles, zwitterions)]  
    carba_pairs_mols = [clean_over_valence_mol(v) for v in carba_pairs] 

    for carbamates_mol in carbamates_mols:
        zwitt = [x for x in zwitterions_mols if x.HasSubstructMatch(carbamates_mol)][0]
        carba_pair = [x for x in carba_pairs_mols if x.HasSubstructMatch(carbamates_mol)][0]

        outs.append(tuple([zwitt, carbamates_mol, carba_pair]))

    return outs

def clean_over_valence_mol(smile):
    """
    Generates a mol object for a molecule the exceeds explicit valence limitations. Useful for ionized compounds.
    """
    mol = Chem.MolFromSmiles(smile,sanitize=False)    
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
    return mol

def prepare_input(mols):

    #reactants : embed -> opt -> coords
    # products : embed -> opt -> coords

    #rea = smile_mol + zwitterion
    # prods = carbamate + charged smile_mol


    pass

if __name__ == "__main__":
    test_smile = "NCCO"
    print(prepare_smiles(test_smile))