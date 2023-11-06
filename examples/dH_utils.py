
from energy_utils import energy_utils
import sql_utils
from rdkit import Chem
import sqlite3



def compute_dH(smile, list_of_options, cursor):


    results = []

    for options in list_of_options:
        
        mol_energy_object = energy_utils(smile = smile, options=options, cursor=cursor)
        if results:
            mol_energy_object.misc_results = results[-1][0]
            mol_energy_object.reac_results = results[-1][1]
            mol_energy_object.prod_results = results[-1][2]

        reactant_energy = mol_energy_object.compute_and_weight_energy()

        product_energies = mol_energy_object.compute_amine_products()

        miscs = mol_energy_object.get_miscs()

        results.append([miscs, reactant_energy, product_energies])


    dHs_list = compute_dH_list(reactant_energy, product_energies, miscs)
    return dHs_list


def compute_dH_list(reactant_energy, product_energies, miscs):
    """
    Generates dHs values given 
    """
    print("Reactants: ", reactant_energy[0], miscs["[H]O[H]"][0], miscs["O=C=O"][0])
    print("products: ", miscs["[H]OC(=O)[O-]"][0], [product[1][0] for product in product_energies])
    reactants = reactant_energy[0] + miscs["[H]O[H]"][0] + miscs["O=C=O"][0]
    products_possibilites = [miscs["[H]OC(=O)[O-]"][0]+ product[1][0] for product in product_energies]

    return [ abs(pos_product - reactants) for pos_product in products_possibilites]

def multi_step():
    """
    
    """

    pass
    

if __name__ == "__main__":


    database_path = "examples/molecules_data.db"
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    
    list_of_options = [{"program":"xtb","method":"gfn_2", "opt":"tight", "solvation":"gbsa", "solvent":"water"}]
    smile = "NCCN"

    dH_list = compute_dH(smile, list_of_options, cursor)
    print("outs ", list(map( energy_utils.hartree_to_kjmol ,dH_list)))

    conn.commit()
    conn.close()