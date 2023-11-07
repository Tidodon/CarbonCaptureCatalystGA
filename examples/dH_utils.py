
from energy_utils import energy_utils
import sql_utils
from rdkit import Chem
import sqlite3



def compute_dH(smile, list_of_options, cursor):


    results = []


    for options in list_of_options:
        
        mol_energy_object = energy_utils(smile = smile, options=options, cursor=cursor)

        reac_confs, reac_atoms = [],[]
        ## Assigns the computed data from previous step. 
        if results:
            mol_energy_object.misc_results = results[-1][0]
            mol_energy_object.reac_results = results[-1][1]
            mol_energy_object.prod_results = results[-1][2]

        if mol_energy_object.reac_results:
            _, reac_confs, reac_atoms = mol_energy_object.reac_results

        reactant_energy = mol_energy_object.compute_and_weight_energy(precomputed_confs=reac_confs, precomputed_atoms=reac_atoms)

        product_energies = mol_energy_object.compute_amine_products()

        miscs = mol_energy_object.get_miscs()

        results.append([miscs, reactant_energy, product_energies])
    #print("RESULTS")
    # for res in results:
    #     for r in res:
    #         print(r)
    dHs_list = compute_dH_list(reactant_energy, product_energies, miscs)
    print("dHs list: ", dHs_list)
    return dHs_list


def compute_dH_list(reactant_energy, product_energies, miscs):
    """
    Generates dHs values given 
    """

    reactants = reactant_energy[0] + miscs["[H]O[H]"][0] + miscs["O=C=O"][0]

    products_possibilites = [(name, miscs["[H]OC(=O)[O-]"][0] + product[0]) for name, product in product_energies.items()]

    return [ (pos_product[0], abs(pos_product[1] - reactants)) for pos_product in products_possibilites]

    

if __name__ == "__main__":


    database_path = "examples/molecules_data.db"
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    
    list_of_options = [{"program":"xtb","method":"gfn_1", "opt":"tight", "solvation":"gbsa", "solvent":"water"},
                       {"program":"xtb","method":"gfn_2", "opt":"tight", "solvation":"alpb", "solvent":"water"}]
    smile = "NCCO"

    dH_list = compute_dH(smile, list_of_options, cursor)
    dHs_raw = [dh[1] for dh in dH_list]
    # print("outs ", list(map( energy_utils.hartree_to_kjmol ,dHs_raw)))

    conn.commit()
    conn.close()
