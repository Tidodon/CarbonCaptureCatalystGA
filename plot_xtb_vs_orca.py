import sqlite3
import os
import pandas as pd
import sys
from rdkit import Chem

current_path = os.getcwd() # outputs a string

if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
    sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/")
    sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/catalystGA")

elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
    sys.path.append("/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/catalystGA")
    sys.path.append("/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/")
else:
    print("Path is different than testing or running environemnt")
    print("The path is: ", current_path)

from sql_utils import check_if_in_db

from dH_utils import compute_dH_list





if __name__ == "__main__":
   
    current_path = os.getcwd()

    database_path = ""
    amines_csv_path  = ""

    database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/molecules_data.db'
    if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
        database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/molecules_data.db'
        amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines_with_ions.csv"

    elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
        database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/molecules_data.db"
        amines_csv_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines_with_ions.csv"

    else:
        print("Outside predefined working directories.")


    molecules = pd.read_csv(amines_csv_path)

    patts                = [Chem.MolFromSmarts("[ND1]"),Chem.MolFromSmarts("[ND2]"),Chem.MolFromSmarts("[ND3]")]
    repls                = [Chem.MolFromSmarts("[NH3+]"),Chem.MolFromSmarts("[NH2+]"),Chem.MolFromSmarts("[NH+]")]#
    misc_smiles          = ("O=C=O", "[H]O[H]", "[H]OC(=O)[O-]")

    misc_xtb_data = {}
    xtb_dH_data = []

    #### xtb
    for misc in misc_smiles:
        o = check_if_in_db(database_path=database_path, smile=misc, method="gfn_2", solvation="alpb", table="misc")
        if o:
            misc_xtb_data.append(o)


    for smile in molecules["SMILES"]:
        rea = check_if_in_db(database_path=database_path, smile=smile, method="gfn_2", solvation="alpb", table="reactants")
        pro = []
        mol = Chem.MolFromSmiles(smile)
        products = []
        for patt, repl in zip(patts, repls):
            prods = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt, replacement=repl)
            products += list(set([Chem.MolToSmiles(p) for p in prods]))

        for p_smile in products:
            pro = check_if_in_db(database_path=database_path, smile=p_smile, method="gfn_2", solvation="alpb", table="products")

