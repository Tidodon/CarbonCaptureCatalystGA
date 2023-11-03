
from energy_utils import energy_utils
import sql_utils
from rdkit import Chem
import sqlite3



def compute_dH(smile, list_of_options, cursor):


    results = []
    
    for options in list_of_options:
        
        mol_energy_object = energy_utils(smile = smile, options=options, cursor=cursor)

        reactant_energy = mol_energy_object.compute_and_weight_energy()

        product_energies = mol_energy_object.compute_amine_products()

        miscs = mol_energy_object.get_miscs()



        pass

def multi_step():

    pass
    

if __name__ == "__main__":


    database_path = "."
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()


    dH_list = compute_dH(smile, list_of_options, cursor)


    conn.commit()
    conn.close()