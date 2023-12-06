
from energy_utils import energy_utils
import sql_utils
from rdkit import Chem
import sqlite3



def compute_dH_data(smile, list_of_options, database_path):

    results = []

    for options in list_of_options:
        print("step: ", options)
        
        mol_energy_object = energy_utils(smile = smile, options=options, database_path=database_path)

        reac_confs, reac_atoms = [], []
        ## Assigns the computed data from previous step. 
        print("results: ", results)
        if results:
            
            mol_energy_object.reac_results = results[-1][0]
            mol_energy_object.prod_results = results[-1][1]
            mol_energy_object.misc_results = results[-1][2]

        if mol_energy_object.reac_results:
            print("mol_energy_object.reac_results: ",  mol_energy_object.reac_results)
            _, reac_confs, reac_atoms = mol_energy_object.reac_results[smile]
            for i,val in  enumerate(mol_energy_object.reac_results[smile]):
                print(i, val)
            print("reac_atoms :", reac_atoms)

        ## Outputs a dictionary of reactant data.
        reactant_energy = mol_energy_object.compute_and_weight_energy(mol=smile, precomputed_confs=reac_confs, precomputed_atoms=reac_atoms, table="reactants")
        
        ## Outputs a list of dictionaries one for each product.
        product_energies = mol_energy_object.compute_amine_products() ## Ou

        ## outputs a list of dictionaries one for each misc:
        miscs = mol_energy_object.get_miscs() 

        results.append([ reactant_energy, product_energies, miscs])

    #dHs_list = compute_dH_list(smile, reactant_energy, product_energies, miscs)

    return results


def compute_dH_list(smile, reactant_energy, product_energies, miscs_energies):
    
    """
    Generates dHs values given dictionaries holding reactant, product and misc molecules information.
    """
    print(reactant_energy[smile][0], miscs_energies["[H]O[H]"][0], miscs_energies["O=C=O"][0], miscs_energies["[H]OC(=O)[O-]"][0], [product[1] for product in product_energies.items() ])
    reactants = reactant_energy[smile][0] + miscs_energies["[H]O[H]"][0] + miscs_energies["O=C=O"][0]

    products_possibilites = [(name, miscs_energies["[H]OC(=O)[O-]"][0] + product[0]) for name, product in product_energies.items()]

    return [ (pos_product[0], abs(pos_product[1] - reactants)) for pos_product in products_possibilites]

    

if __name__ == "__main__":


    database_path = "examples/molecules_data.db"

    list_of_options = [{"program":"xtb","method":"gfn_1", "opt":"tight", "solvation":"gbsa", "solvent":"water"},
                       {"program":"xtb","method":"gfn_2", "opt":"tight", "solvation":"alpb", "solvent":"water"}]
    smile = "[K+].[O-]C(=O)[C@@H]1CCCN1"

    results = compute_dH_data(smile, list_of_options, database_path=database_path)
    # dHs_raw = [dh[1] for dh in dH_list]
    # print("outs ", list(map( energy_utils.hartree_to_kjmol ,dHs_raw)))
    
    
    sql_utils.insert_result_to_db( results, list_of_options, database_path=database_path)

    sql_utils.print_table_contents("miscs", "reactants", "products", database_path=database_path)

