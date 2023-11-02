#Module for recovering and and computing energies.
import os
import sys
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

from catalystGA import GA
from catalystGA.reproduction_utils import graph_crossover, graph_mutate
from catalystGA.utils import MoleculeOptions
import sql_utils
from rdkit import Chem
from rdkit.Chem import AllChem
from xtb import xtb_calculate
from orca import orca_calculate
import math
import sqlite3

class energy_utils:
    
    def __init__(self, smile: str, options: dict, cursor: sqlite3.Connection) -> None:
        self.smile       = smile
        self.options     = options
        self.cursor      = cursor
        self.misc_smiles = ["O=C=O", "[H]O[H]", "[H]OC(=O)[O-]"]
        self.K_B         = 3.166811563 * math.pow(10,-6)         # E_h/K
        self.T_K         = 313                                   # Kelvin
        self.n_cores     = 1                                     # Cores to be used in the computation.

    def get_miscs(self):
        """
        Overview function to recover miscalleneous molecules data,carbon dioxide water hydrogen-carbonate, from 
        database if precomputed or 
        """

        method, solvation = self.options['method'], self.options['solvation']
        miscs_data = []
        misc_check = sql_utils.check_if_miscs_in_db(self.cursor, method=method, solvation=solvation)
        if misc_check:
            miscs_data = sql_utils.recover_miscs(misc_check)
        else:
            miscs_data = self.compute_miscs()

        return miscs_data

    def compute_miscs(self):

        miscs_dict = {}
        
        for misc_smile in self.misc_smiles:
            name_energy, name_coords = self.compute_and_weight_energy()
            miscs_dict[misc_smile] = name_energy
            miscs_dict["opt_conformer_coords_of_"+misc_smile] = name_coords
        return miscs_dict

    def compute_and_weight_energy(self, separate=True):
        """
        Overhead to deal with potential ionic compounds and compute their energies and coordinates.
        
        separate : condition for whether ions should be optimized together or separates
        """

        if "." in self.smile:
            sub_smiles = self.smile.split(".")
        else:
            sub_smiles = [self.smile]

        tot_e = 0
        maxsub = max(sub_smiles, key=len) #Picks the largest ionic part of the molecule or simply the molecule itself if it is not ionic.
        conf_coords = []

        for sub_smile in sub_smiles:

            confs = self.calculate_energy(sub_smile)
            if sub_smile == maxsub:
                 conf_coords = confs
            tot_e += self.weight_energy(conf_coords)
            
            return tot_e, [conf[1] for conf in conf_coords]

    def weight_energy(self, confs):
        """
        Computes the weighted energy of a molecule. Weights are dependent on the temperature and conformers included.
        adjustment parameter makes sure that Boltzmann exponent do not cause overflow.
        E : the energy weigthed by Boltzmann parameters.
        Z : Normalization factor from the sum of Boltzmann factors.
        """

        adjustment = confs[0][2] 

        boltz_exponents = [((val[2]-adjustment))/(self.K_B * self.T_K) for val in confs ]

        boltzmann_pop_reactants = [math.exp(-boltz_exp) for boltz_exp in boltz_exponents]

        Z = sum(boltzmann_pop_reactants)
        E = sum([reactant_pop*conf_e[2] for reactant_pop, conf_e in zip(boltzmann_pop_reactants, confs)])

        return E/Z
        
    def calculate_energy(self, mol_smile):
        ###Computes an energy for a mol object defined by its SMILES/SMARTS string. 
        # The energy is weighted by the contribution of individual conformers.

        mol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(mol_smile)))

        threshold = 0.5
        if  Chem.RemoveHs(self.mol).GetNumAtoms() > 9:
            threshold = 1.0

        # if Chem.MolFromSmiles(mol).GetNumAtoms() == 1 :
        #     print("SUBOPTIONS: ", mol, self.options)
        #     _ = self.options.pop("opt")

        _ = Chem.rdDistGeom.EmbedMultipleConfs(
                    mol,
                    clearConfs=True,
                     # maxAttempts = 10,
                    numConfs=500,
                    useRandomCoords=True,
                    pruneRmsThresh=threshold,
                    #randomSeed=5
                    )
        
        optimized_confs = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant='MMFF94')

        ##### Conformer pruning based on whether optimization converged and energy magnitude.
        num_confs = self.mol.GetNumConformers()

        min_e_conf = min([optimized_confs[i][1] for i in range(num_confs)])
        for i in range(num_confs-1, -1, -1):
            #print(optimized_confs[i])
            if  bool(optimized_confs[i][0]):
                self.mol.RemoveConformer(i)
            elif (optimized_confs[i][1] - min_e_conf) > 10:
                #print((optimized_confs[i][1] - min_e_conf))
            self.mol.RemoveConformer(i)

        num_confs = self.mol.GetNumConformers()



        # Cluster conformers
        # diffmat = AllChem.GetConformerRMSMatrix(self.mol, prealigned=False)
        # clt = Butina.ClusterData(diffmat, num_confs, distThresh=0.05,
        #                      isDistData=True, reordering=True)
        # # Get unique conformers
        # best_conformers =[ conformer for conformer in (self.mol).GetConformers()]
        # centroid_idx = [c[0] for c in clt] # centroid indexes
        # print(f"Centroids for clustering are {centroid_idx}")
        # unique_best_conformers = [best_conformers[i] for i in centroid_idx]
        # print( f"Number of conformers of {Chem.MolToSmiles(self.mol)} post cluster pruning is {len(unique_best_conformers)}.")
        atoms = [atom.GetSymbol() for atom in self.mol.GetAtoms()]
        
        if self.program =="xtb":
            xtb_options = self.prepare_xtb_options()
            print("self options: ", self.options)
            try:
                xtb_options["charge"] =  self.options.pop("charge")
            except:
                xtb_options["charge"] = Chem.rdmolops.GetFormalCharge(self.mol)
            print("XTB options: ", xtb_options)
            try:
                res = [xtb_calculate(atoms=atoms, coords=conformer.GetPositions(), options=xtb_options, n_cores=n_cores) for conformer in (self.mol).GetConformers()]
                return res
            except:
                print("Incorrect termination of XTB.")
                print(self.smiles, self.options, xtb_options)
            return [[atoms, (self.mol).GetConformers()[0].GetPositions(), 1000000]]

        elif self.program == "orca":
            # charge=0
            # try: 
            #     charge = self.options.pop("charge") #Removes Charge key/value and returns the value to be used in orca_calculate
            # except:
            #     charge = Chem.rdmolops.GetFormalCharge(self.mol)#0
            charge = Chem.rdmolops.GetFormalCharge(self.mol)
            orca_options = self.prepare_orca_options()
            print("orca options: ", orca_options)
            #### Prepare orca output to same format as xtb output:
            try:
                res = [orca_calculate(atoms=atoms, coords=conformer.GetPositions(), options=orca_options, n_cores=n_cores, charge=charge) for conformer in (self.mol).GetConformers()]
                return [[v['atoms'], v['opt_coords'], v['electronic_energy']] for v in res]
            except:
                print("Incorrect termination of Orca. -> atoms/opt_coords/electronic_energy dict keys don't respond")
                print(self.smiles, self.options, orca_options)
                return [[atoms, (self.mol).GetConformers()[0].GetPositions(), 1000000]]
        else:
            raise "Incorrect specification of the QM program."
    
if __name__ == "__main__":
    pass

