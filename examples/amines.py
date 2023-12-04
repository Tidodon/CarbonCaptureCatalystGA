import math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.ML.Cluster import Butina
import copy
import sys
import os


#import ts_utils

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
### Modules
from dH_utils import compute_dH_data, compute_dH_list
from sql_utils import insert_result_to_db
import sqlite3

#### TODOS:
#### - Intermediate naming for prim/seco/tert amines in the scoring function
#### - dH scoring for A.B ionic compounds. Retrieve reactant data to compute.
#### - Deal with possible duplicates in DB: non-exact but similar namings, two rows with identical methods used but different energies.
#### - Pack Misc values recovery into a single function.
#### - ga.run outputs -> mol, score, dH, k, "which k?"
#### - amine_products_all -> order by energies to only recover relevant protonation sites.
#### - CO2_energy, H2O_energy, OCOO_energy = self.compute_misc_energies() -> Implement as a separate function.


class AmineCatalyst:
    save_attributes = {}  # any other attributes to save to the database


    ### Boltzmann constant
    K_B = 3.166811563 * math.pow(10,-6) # E_h/K

    #The checks are to be performed with molecules with implicit hydrogens. (-> NH+ presence?)
    patts = [Chem.MolFromSmarts("[ND1]"),Chem.MolFromSmarts("[ND2]"),Chem.MolFromSmarts("[ND3]")]
    #
    repls =  [Chem.MolFromSmarts("[NH3+]"),Chem.MolFromSmarts("[NH2+]"),Chem.MolFromSmarts("[NH+]")]#

    def __init__(self, mol: Chem.Mol) -> None:
        self.mol       = mol
        self.list_of_options   = [ ] # Specify program, method , solvation, opt, solvent. One dict of those for each computation step.
        self.database_path = ""
        self.score     = math.nan
        self.fitness   = math.nan
        self.timing    = math.nan
        self.T_K       = 313 # Kelvin
        self.error     = ""
        self.idx       = (-1, -1)
        self.amine_type = tuple(True if mol.HasSubstructMatch(patt) else False for patt in self.patts)#Respectively primary/secondary/tertiary amine WITHOUT explicit hydrogens.
        self.dHabs = math.nan #Heat of absorbtion
        self.kabs = math.nan #k of reaction limiting step. amine->bicarbonate for tertiary amines
        self.results = [] #A dictionary with possible keys: "miscs", "reactant", "products", each will be only generated
        #if the corresponding values are not found in the databse.

    @property
    def smiles(self) -> str:
        """Yields SMILES string of molecule, needed for Database.

        Returns:
            str: SMILES string
        """
        return Chem.MolToSmiles(Chem.RemoveHs(self.mol))  

    @staticmethod
    def hartree_to_kcalmol(hartree) -> float:
        ### Conversion value take from wiki
        return hartree * 627.509474063

    @staticmethod
    def hartree_to_ev(hartree) -> float:
        ### Conversion value take from wiki
        return hartree * 27.211386245988 
    
    @staticmethod
    def kcalmol_to_kjmol(kcalmol) -> float:
        ### Conversion value taken from wiki.
        return kcalmol * 4.184
    
    @staticmethod
    def hartree_to_kjmol(hartree) -> float:
        ### Hartree to Joule value from NIST, and Avogadro's number taken from wiki
        joule = 4.3597447222071 * 10**(-18) #Joule/Hartree
        Na = 6.02214076 * 10**23 # 1/mol
        return hartree * joule * Na * 0.001 
    
 
    @staticmethod
    def chk_conn(conn):
        try:
            conn.cursor()
            return True
        except Exception as ex:
            return False
        
    def calculate_score(
        self, n_cores: int = 4, envvar_scratch: str = "SCRATCH", scoring_kwargs: dict = {}
    ):
        """Calculate score of molecule, store in self.score.

        Args:
            n_cores (int, optional): Number of cores to use when run on cluster. Defaults to 1.
            envvar_scratch (str, optional): Name of environmental variable pointing to scratch directory. Defaults to 'SCRATCH'.
            scoring_kwargs (dict, optional): Additional keyword agruments parsed to scoring function. Defaults to {}.
        """
        
        self.results = compute_dH_data(database_path=self.database_path, smile=self.smiles, list_of_options=self.list_of_options)

        reactant_energy, product_energies, miscs = self.results[-1]

        dHs = compute_dH_list(smile=self.smiles, reactant_energy=reactant_energy, product_energies=product_energies, miscs_energies=miscs)

        #### Alternatively could be chosen based on the highest k value.
        self.dHabs = max(dHs, key=lambda x :x[1])

        #AmineCatalyst.hartree_to_kjmol(
        #results_list = [[reactant_energy], [ ]for prod in prods]

        # print("Reactant energy: ", self.weight_energy(reactant_confs))

        ### Decide on which product to use by k value:

        #Assign score values based on dH, k, SA

        #dH scorings alone.

        self.score = self.dHabs[1]

        #### Later a reordering code will be added here.

        #order_amine_products(dH, dG) -> top three most reactive. 

        ####
        #logP = Descriptors.MolLogP(self.mol)
        #self.score = logP

class GraphGA(GA):
    def __init__(
        self,
        mol_options,
        population_size,
        n_generations,
        mutation_rate,
        scoring_kwargs,
        db_location,
        comp_options,
        db_mols_path, 
    ):
        super().__init__(
            mol_options     = mol_options,
            population_size = population_size,
            n_generations   = n_generations,
            mutation_rate   = mutation_rate,
            db_location     = db_location,
            scoring_kwargs  = scoring_kwargs,
        )
        self.comp_options = comp_options
        self.db_mols_path = db_mols_path

    def make_initial_population(self):
        amine_pops_path = amines_csv_path
        with open(amine_pops_path, "r") as f: #data -> csv
            lines = f.readlines()[1:]
        mols = [Chem.MolFromSmiles(line.split(",")[0]) for line in lines]
        population = [AmineCatalyst(mol) for mol in mols[: self.population_size]]
        for amine in population:
            amine.list_of_options = self.comp_options
            amine.database_path = self.db_mols_path
        return population

    def crossover(self, ind1, ind2):
        mol1 = ind1.mol
        mol2 = ind2.mol
        new_mol = None
        while not new_mol:
            new_mol = graph_crossover(mol1, mol2)
        try:
            Chem.SanitizeMol(new_mol)
            ind = AmineCatalyst(new_mol)
            return ind
        except Exception:
            return None

    def mutate(self, ind):
        mol = ind.mol
        new_mol = None
        while not new_mol:
            new_mol = graph_mutate(mol)
        try:
            Chem.SanitizeMol(new_mol)
            ind = AmineCatalyst(new_mol)
            return ind
        except Exception:
            return None
        
    def run(self):
        results = []  # here the best individuals of each generation will be stored
        self.print_parameters()
        
        self.population = self.make_initial_population()
        
        # self.miscs = [get_miscs(options) for options in list_of_options]

        ### save to db.
        ####
        print("Population in run: ", self.population)

        # for pop in self.population:
        #     pop.calculate_score()
        self.population = self.calculate_scores(self.population, gen_id=0)
        for pop in self.population:
            print(Chem.MolToSmiles(pop.mol))
            print(pop.dHabs)
            insert_result_to_db(database_path=self.db_mols_path, results=pop.results, list_of_options=self.comp_options)
        
        #self.add_computed_pops_to_db()

        results = self.population
        # self.db.add_individuals(0, self.population)
        
        # self.print_population(self.population, 0)
        
        # for n in range(0, self.n_generations):
        #     print("N-generation: ", n, "\n")
        #     self.calculate_fitness(self.population)
        #     self.db.add_generation(n, self.population)
        #     self.append_results(results, gennum=n, detailed=True)
        #     children = self.reproduce(self.population, n + 1)
        #     children = self.calculate_scores(children, gen_id=n + 1)
        #     self.db.add_individuals(n + 1, children)
        #     self.population = self.prune(self.population + children)
        #     self.print_population(self.population, n + 1)
        # self.calculate_fitness(self.population)
        # self.db.add_generation(n + 1, self.population)
        # self.append_results(results, gennum=n + 1, detailed=True)
        
        return results
    
    @staticmethod
    def plot_dH_vs_dH(exp_dH, calc_dH, options, title="", figname="", xlab="", ylab=""):
        assert type(options) == dict
        plt.scatter(exp_dH, calc_dH, marker="o", color="b")
        if xlab:
            plt.xlabel(xlab)
        else:
            plt.xlabel("Experimental " + r"$ \Delta H $")

        if ylab: 
            plt.ylabel(ylab)
        else:
            plt.ylabel("Calculated "   + r"$ \Delta H $")
        plt.axline([abs(min(exp_dH + calc_dH)),abs(min(exp_dH+calc_dH))], slope=1)
        slope, intercept, r, p, se = stats.linregress(exp_dH, calc_dH)
        R2 = r**2

        def reg_pnts(x):
            return slope * x + intercept
    
        xs = np.linspace(0,100, 100)

        pnts = [reg_pnts(pnt) for pnt in xs]
        pnts = [pnt for pnt in pnts if pnt<=max(exp_dH+calc_dH)]# and pnt>min(exp_dH+calc_dH)]

        plt.plot(xs[:len(pnts)], pnts, ls="--", color="grey", label=f'slope: {slope:.2f}, intercept:+ {intercept:.2f}')
    
        plt.plot([], [], label=f'$R^{2}$: {R2:.2f}')
        plt.plot([], [], label=f'stderr: {se:.2f}')
        plt.xlim(0,max(exp_dH+calc_dH) )
        plt.ylim(0,max(exp_dH+calc_dH) )
        if title:
            plt.title(title)
        else:
            plt.title(f"{options['method']} {options['solvation']}")
        plt.legend()

        if figname:
            ## Use the specified figname.
            pass
        else:
            try:
                figname = "_".join([str(val) for val in options.values()])+".eps"
            except: 
                figname = "" + ".eps"
        plt.savefig(figname, format='eps')
        #plt.show()
        plt.close()

if __name__ == "__main__": 
    import numpy as np 
    from scipy import stats 
    #import time
    import pandas as pd 
    import matplotlib.pyplot as plt 

    ##Get paths to amines and database.
    database_path = ""
    amines_csv_path  = ""

    current_path = os.getcwd()
    if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
        amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines_with_ions.csv"
        #amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/molecules_data.db'
    elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
        amines_csv_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines_with_ions.csv"
        #amines_csv_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/molecules_data.db"
    else:
        print("Path is different than testing or running environemnt")
        print("The path is: ", current_path)

    
    amines = pd.read_csv(amines_csv_path)
    #exp_dH = amines.loc[amines['SMILES'].isin(calc_names)]['dH'].tolist()

    calc_dH, exp_dH = [], []

    cnt = 0
    names, dHs = [],[]

    list_of_options = [{"program":"xtb","method":"gfn_2", "opt":True, "solvation":"alpb", "solvent":"water"}]
                       #{"program":"orca","method":"r2SCAN-3c", "solvation":"CPCM", "solvent":"water"}]#,
                      # 

    ga = GraphGA(
        mol_options=MoleculeOptions(AmineCatalyst),
        population_size=32,
        n_generations=1,
        mutation_rate=0.0,
        db_location="organic.sqlite",
        scoring_kwargs={},
        comp_options=list_of_options,
        db_mols_path = database_path
    )

    # m = AmineCatalyst(Chem.MolFromSmiles("[NH3+]CCO"))#112.34863236070932
    # m.options = comp_options
    # m.program = comp_program
    # m.calculate_energy(n_cores=1)
    
    # print("Computed score: ", Chem.MolToSmiles(m.mol), m.score)
    #print(res)
    results = []
    #results = ga.run()

    results_xtb = ga.run()

    list_of_options = [{"program":"orca","method":"r2SCAN-3c", "solvation":"CPCM", "solvent":"water"}]#,
                      # 

    ga = GraphGA(
        mol_options=MoleculeOptions(AmineCatalyst),
        population_size=32,
        n_generations=1,
        mutation_rate=0.0,
        db_location="organic.sqlite",
        scoring_kwargs={},
        comp_options=list_of_options,
        db_mols_path = database_path
    )
    results_orca = ga.run()


    ##########################################################
    ###Temporary code for benchmarking dH computations.#######
    ##########################################################

    calc_names, calc_dH = [],[]
    for molecule in results_orca:
        print("molecuel: ", Chem.MolToSmiles(molecule.mol))
        calc_names.append(Chem.MolToSmiles(molecule.mol))
        calc_dH.append(AmineCatalyst.hartree_to_kjmol(molecule.dHabs[1]))
    ##########################################################


    ##########################################################
    ################Plotting preparation.####################
    ##########################################################
    exp_dH = amines.loc[amines['SMILES'].isin(calc_names)]['dH'].tolist()
    exp_names = amines.loc[amines['SMILES'].isin(calc_names)]['SMILES'].tolist()

    exp_df = pd.DataFrame({"SMILES":exp_names, "dH_exp":exp_dH})
    calc_df = pd.DataFrame({"SMILES":calc_names, "dH_calc":calc_dH})

    dH_df = pd.merge(calc_df, exp_df, on="SMILES")

    xtb_names, xtb_dH = [],[]
    for molecule in results_xtb:
        print("molecuel: ", Chem.MolToSmiles(molecule.mol))
        xtb_names.append(Chem.MolToSmiles(molecule.mol))
        xtb_dH.append(AmineCatalyst.hartree_to_kjmol(molecule.dHabs[1]))

    xtb_df = pd.DataFrame({"SMILES":xtb_names, "dH_xtb":xtb_dH})

    orca_names, orca_dH = [],[]
    for molecule in results_orca:
        print("molecuel: ", Chem.MolToSmiles(molecule.mol))
        orca_names.append(Chem.MolToSmiles(molecule.mol))
        orca_dH.append(AmineCatalyst.hartree_to_kjmol(molecule.dHabs[1]))

    orca_df = pd.DataFrame({"SMILES":orca_names, "dH_orca":orca_dH})


    dH_df_orca_xtb = pd.merge(xtb_df, orca_df, on="SMILES")

    # print(dH_df)

    #exp_dH = [ v[2] for v in ]

    #generations = [r[0] for r in results]
    #best_scores = [max([ind.dH for ind in res[1]]) for res in results]
    #calc_dH = [max([ind.dH for ind in res[1]]) for res in results]

    GraphGA.plot_dH_vs_dH(dH_df_orca_xtb["dH_xtb"], dH_df_orca_xtb["dH_orca"], ga.comp_options[-1], title="gfn_2(alpb) vs r2SCAN-3c(CPCM)", xlab="xtb", ylab="orca")


    GraphGA.plot_dH_vs_dH(dH_df["dH_exp"], dH_df["dH_calc"], ga.comp_options[-1])


    # fig, ax = plt.subplots()
    # ax.plot(generations, best_scores)
    # ax.set_xlabel("Generation")
    # ax.set_ylabel("Max Score")

    #plt.savefig("organic.png")