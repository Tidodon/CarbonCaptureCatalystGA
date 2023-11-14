


#Module for recovering and and computing energies.
import os
import sys
import numpy as np
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

"""
TODO:
- implement combined ion computation functionality
"""
class energy_utils:

    patts                = [Chem.MolFromSmarts("[ND1]"),Chem.MolFromSmarts("[ND2]"),Chem.MolFromSmarts("[ND3]")]
    repls                = [Chem.MolFromSmarts("[NH3+]"),Chem.MolFromSmarts("[NH2+]"),Chem.MolFromSmarts("[NH+]")]#
    misc_smiles          = ("O=C=O", "[H]O[H]", "[H]OC(=O)[O-]")
    K_B                  = 3.166811563 * math.pow(10,-6)         # E_h/K
    
    def __init__(self, smile: str, options: dict, database_path: str) -> None:
        self.smile         = smile
        self.options       = options 
        self.database_path = database_path
        self.T_K           = 313                                   # Kelvin
        self.n_cores       = 1                                     # Cores to be used in the computation.
        self.amine_type    = tuple(True if Chem.MolFromSmiles(self.smile).HasSubstructMatch(patt) else False for patt in self.patts)#Respectively primary/secondary/tertiary amine WITHOUT explicit hydrogens.
        self.misc_results  = {}
        self.reac_results  = {}
        self.prod_results  = {}

    @staticmethod
    def hartree_to_kjmol(hartree) -> float:
        ### Hartree to Joule value from NIST, and Avogadro's number taken from wiki
        joule = 4.3597447222071 * 10**(-18) #Joule/Hartree
        Na = 6.02214076 * 10**23 # 1/mol
        return hartree * joule * Na * 0.001  
    
    def get_miscs(self):
        """
        Overview function to recover miscalleneous molecules data,carbon dioxide water hydrogen-carbonate, from 
        database if precomputed or 
        """

        for misc_smile in self.misc_smiles:
            self.misc_results    =  self.misc_results | self.compute_and_weight_energy(mol=misc_smile, precomputed_confs=[], precomputed_atoms=[], table="miscs")
        
        return self.misc_results

    def compute_and_weight_energy(self, mol = "", separate=True, precomputed_confs=[], precomputed_atoms=[], table=""):
        """
        Overhead to deal with potential ionic compounds and compute their energies and coordinates.
        
        separate : condition for whether ions should be optimized together or separates NOT implemented yet.

        Outputs: 
        - total energy boltzmann energy weighted by Boltzmann factors
        - list of conformer coordinates, if the compound was ionic only the largest part is retained here.
        """

        mol_check = sql_utils.check_if_in_db(database_path=self.database_path, smile=mol, method=self.options["method"], solvation=self.options["solvation"], table=table)
        if mol_check:
            print("Retrieveing energy!")
            return mol_check


        if "." in mol:
            sub_smiles = mol.split(".")
        else:
            sub_smiles = [mol]

        tot_e = 0
        maxsub = max(sub_smiles, key=len) #Picks the largest ionic part of the molecule or simply the molecule itself if it is not ionic.
        conf_coords = []

        for sub_smile in sub_smiles:

            confs = self.calculate_energy(sub_smile, precomputed_confs=precomputed_confs, precomputed_atoms=precomputed_atoms)
            if sub_smile == maxsub:
                 conf_coords = confs
            tot_e += self.weight_energy(conf_coords)

            return { mol : (tot_e, [conf[1] for conf in conf_coords], [conf[0] for conf in conf_coords])}#[atom.GetSymbol() for atom in Chem.MolFromSmiles(maxsub).GetAtoms()]

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
    
    @staticmethod
    def get_conformers(mol_smile, thresholds=[(9, 1.0),(0, 0.4)], energy_cutoff=10) -> list:
        """
        thresholds: a list of tuples defining cutoffs for conformer generation. I assume they are ordered from largest to smallest.
        First  index: length of the molecules over which threshold applies
        Second index: cutoff threshold

        Important for larger molecules that easily can have 1000's of conformers.
        """
        mol     = Chem.MolFromSmiles(mol_smile)
        molsize = Chem.RemoveHs(mol).GetNumAtoms()
        mol     = Chem.AddHs(mol)

        for threshold in thresholds:
            if threshold[0] < molsize:
                threshold_set = threshold[1]

        _ = Chem.rdDistGeom.EmbedMultipleConfs(
                    mol,
                    # clearConfs=True,
                    # maxAttempts = 10,
                    numConfs=500,
                    useRandomCoords=True,
                    pruneRmsThresh=threshold_set,
                    #randomSeed=5
                    )
        
        optimized_confs = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant='MMFF94')

        num_confs = mol.GetNumConformers()

        min_e_conf = min([optimized_confs[i][1] for i in range(num_confs)])
        for i in range(num_confs-1, -1, -1):
            if  bool(optimized_confs[i][0]):
                mol.RemoveConformer(i)
            elif (optimized_confs[i][1] - min_e_conf) > energy_cutoff:
                mol.RemoveConformer(i)

        num_confs = mol.GetNumConformers()

        return [conformer.GetPositions() for conformer in mol.GetConformers()], [atom.GetSymbol() for atom in mol.GetAtoms()]
        
    def calculate_energy(self, submol_smile, precomputed_confs=[], precomputed_atoms=[]):
        ###Computes an energy for a mol object defined by its SMILES/SMARTS string. 
        # The energy is weighted by the contribution of individual conformers.

        mol         = Chem.AddHs(Chem.MolFromSmiles(submol_smile))

        submol_size = mol.GetNumAtoms()
        
        conformers  = []
        atoms       = []

        if precomputed_confs:
            conformers, atoms = precomputed_confs, precomputed_atoms[0] ##just to keep the ordering
            conformers = [ np.array(arr) for arr in conformers]
            print("Precomputed conformers: ", conformers, atoms)
        else:
            conformers, atoms = self.get_conformers(submol_smile)
            print("First confoermers: ", conformers, atoms)
        
        if self.options["program"] =="xtb":
            xtb_options = self.prepare_xtb_options(submol_size)
            # print("self options: ", self.options)
            xtb_options["charge"] = Chem.rdmolops.GetFormalCharge(mol)
            # print("XTB options: ", xtb_options)
            try:
                res = [xtb_calculate(atoms=atoms, coords=conformer, options=xtb_options, n_cores=self.n_cores) for conformer in conformers]
                return res
            except:
                print("Incorrect termination of XTB.")
                print(submol_smile, self.options, xtb_options)
            return [[atoms, mol.GetConformerss()[0].GetPositions(), 1000000]]

        elif self.program == "orca":
            # charge=0
            # try: 
            #     charge = self.options.pop("charge") #Removes Charge key/value and returns the value to be used in orca_calculate
            # except:
            #     charge = Chem.rdmolops.GetFormalCharge(self.mol)#0

            charge = Chem.rdmolops.GetFormalCharge(mol)
            orca_options = self.prepare_orca_options(submol_size)
            print("orca options: ", orca_options)
            #### Prepare orca output to same format as xtb output:
            try:
                res = [orca_calculate(atoms=atoms, coords=conformer, options=orca_options, n_cores=self.n_cores, charge=charge) for conformer in conformers]
                #if 'opt_coords' in res[0]:
                return [[v['atoms'], v['opt_coords'], v['electronic_energy']] for v in res]
                    
            except:
                print("Incorrect recovery of Orca output. -> atoms/opt_coords/electronic_energy dict keys don't respond")
                print(self.smiles, self.options, orca_options)
                return [[atoms, (self.mol).GetConformers()[0].GetPositions(), 1000000]]
        else:
            raise "Incorrect specification of the QM program."
        
    def prepare_xtb_options(self, mol_size) -> dict:
        """
        Build an input dicionary for the xtb_calculate function to perfrom xtb calculations.
        """
        xtb_options = {}
        
        try:
            mtd, tp = self.options["method"].split("_")
            xtb_options[mtd]=int(tp)
        except:
            print("Unspecified QM method")

        try:
            xtb_options[self.options["solvation"]] = self.options["solvent"]
        except:
            print("Unspecified solvation")

        if bool(self.options.get("opt")) and (mol_size > 1):
            xtb_options["opt"] = self.options["opt"]

        return xtb_options
    

    def prepare_orca_options(self, mol_size) -> dict:
        orca_options = {}
        try:
            #mtd, tp = self.options["method"].split("_") ## Probably unnecessary for orca calculations.
            orca_options[(f'{self.options["method"]}').upper()] = ""
        except:
            print("Unspecified QM method")

        try:
            orca_options[self.options["solvation"].upper()]=self.options["solvent"].lower().capitalize()
            #options_string += f' {self.options["solvation"]}({self.options["solvent"]})'.lower()
        except:
            print("Unspecified solvation")

        if self.options["opt"] and (mol_size > 1):
            #options_string += ' OPT'
            orca_options["OPT"] = ""

        else:
            print("Unspecified optimization")
        return orca_options
    
    def cat_products(self, patt, repl):
        """
        A generator method that gives smiles representation 
        of the possible product molecules given pattern(patt) 
        and replacement(repl). It gives energy values for each 
        of the products. 

        Arguments:
        patt: recognization pattern given by a mol object 
        repl: replacement of the pattern, given by a mol object. 
        """

        # Sanitization step. DO NOT REMOVE. Otherwise Conformer embedding in self.calculate_energy() breaks:
        mol = Chem.MolFromSmiles(self.smile)

        # Replacement step. It is dependenent on the whether the molecule was sanitized or not.

        products = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt, replacement=repl)

        products = [Chem.MolFromSmiles(m) for m in set([Chem.MolToSmiles(p) for p in products])]

        precomputed_confs = []
        precomputed_atoms = []

        for prod in products:
            print("Check products ",Chem.MolToSmiles(prod))
            prod_smile = Chem.MolToSmiles(prod)
            
            if self.prod_results:
                _, precomputed_confs, precomputed_atoms = self.prod_results[prod_smile]

            yield self.compute_and_weight_energy(prod_smile, precomputed_confs=precomputed_confs, precomputed_atoms=precomputed_atoms, table="products")

    def compute_amine_products(self):

        cats = {}
        if self.amine_type[0]: # If mol has primary amine
            for prod in self.cat_products(patt=self.patts[0], repl=self.repls[0]):
                cats = cats | prod

        if self.amine_type[1]: # If mol has secondary amine
            for prod in self.cat_products(patt=self.patts[1], repl=self.repls[1]):
                cats = cats | prod

        if self.amine_type[2]: # If mol has tertiary amine.
            for prod in self.cat_products(patt=self.patts[2], repl=self.repls[2]):
                cats = cats | prod

        return cats  
    
if __name__ == "__main__":

    current_path = os.getcwd()
    if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
        amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/molecules_data.db'
    elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
        amines_csv_path = "/lustre/hpc/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = "/lustre/hpc/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/molecules_data.db"
    else:
        print("Path is different than testing or running environemnt")
        print("The path is: ", current_path)



    options =  {"program":"xtb", "opt":True, "method":"gfn_2", "solvation":"alpb", "solvent":"water"}
    test = energy_utils(smile="NCCO", options=options, database_path=database_path)
    out = test.compute_and_weight_energy(mol="NCCO", table="reactants")
    amineprot = test.compute_amine_products()

    print(amineprot)

    lst = [[[2,3.5],[45,2]]]
    stringed_list = sql_utils.opt_coords_to_csv_string(lst)
    atoms= ["A", "B"]
    atoms_ord = sql_utils.atoms_ord_to_csv_string(atoms)

    # params=("[NH3+]CCO","-14.9303","gfn_2", "alpb", stringed_list, atoms_ord)
    # cursor.execute("INSERT INTO products (smiles, energy, method, solvation, opt_coords, atoms_ord) VALUES(?,?,?,?,?,?)", params)

    # amine_products_all = test.compute_amine_products()
    # print("amine produxt LL", amine_products_all)
    # print("Computed energy of NCCO: ", test.compute_and_weight_energy())