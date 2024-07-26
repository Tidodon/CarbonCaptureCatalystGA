#TODO:
#-treat tert
import math
import random
import numpy as np 
from scipy import stats 
#import time
import pandas as pd 
import matplotlib.pyplot as plt 


from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import RDConfig
from rdkit.Chem.Draw import IPythonConsole
from rdkit.ML.Cluster import Butina
import copy
import sys
import os

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

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
from k_utils import list_of_dGs, compute_k_rate_from_dG
import energy_utils
import sqlite3
import pandas as pd
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
        self.dG_lst = []
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

    @staticmethod
    def rect_conv(x, sig=1,n =1):
        a = x-sig
        b = x+sig
        return 0.5*(math.erf(n*b)-math.erf(n*a))

    @staticmethod 
    def sigmoid(x):
        return 1 / (1 + math.exp(-x))

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
       
        ###Computed transitions state dG
        self.dG_lst = list_of_dGs(Chem.RemoveHs(self.mol))

        #### Assign atom_ids to protonated amines.
        for i,sm in enumerate(dHs):
            prot_amine = Chem.MolFromSmarts("[#7;+]")
            amine_id = Chem.MolFromSmiles(sm[0]).GetSubstructMatch(prot_amine)
            dHs[i] = [*dHs[i], amine_id[0]]

        #### Match dH/k amines
        dG_ts = min(self.dG_lst, key=lambda x: x[3])
        ts_id = dG_ts[2]
        matched_dH = [] 
        for lst in dHs:
            if lst[-1]==ts_id:
                matched_dH = lst
                break

        self.kabs = compute_k_rate_from_dG(dG_ts[3], self.T_K) 
         
        #### Alternatively could be chosen based on the highest k value.
        self.dHabs = max(dHs, key=lambda x :x[1])
        self.dHabs[1] = AmineCatalyst.hartree_to_kjmol(self.dHabs[1])        
        print(f"self.dHabs: {self.dHabs}") 
        #AmineCatalyst.hartree_to_kjmol(
        #results_list = [[reactant_energy], [ ]for prod in prods]

        ### Decide on which product to use by k value:
        

        print(f"dHs: {dHs}")
        #Assign score values based on dH, k, SA
        slope_k,  intercept_k  = 0.08896083453589532, 0.7866132724304302 #Data from alpb->r2scan opt runs
        slope_dH, intercept_dH = 0.4348138758572645,  53.47577527659122  #Data from 32 amines with alpb/gfn2 opt runs
        adjusted_kabs = self.kabs * slope_k +  intercept_k
        adjusted_dH   = self.dHabs[1] * slope_dH + intercept_dH
        k_cutoff = np.log(15000)
        dH_mid = (40+75)/2
        k_score  = AmineCatalyst.sigmoid(adjusted_kabs - k_cutoff)

        dH_score = AmineCatalyst.rect_conv(adjusted_dH - dH_mid, sig=17.5, n=0.3)

        ri = self.mol.GetRingInfo().BondRings()
        ring_penalty = 0
        if len(ri)>0:
            ring_penalty  = len(max(ri, key=len))/6
        sa_score = sascorer.calculateScore(self.mol)

        print(f"Score components: \nk_score     :{k_score}\ndH_score    :{dH_score}\nring_penalty:{ring_penalty}\n sa_score    :{sa_score}")
        self.score = 5*k_score*dH_score - ring_penalty - sa_score 
        self.dHabs[1] = adjusted_dH 
        self.kabs = adjusted_kabs
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
        random.shuffle(mols) 
        population = [AmineCatalyst(mol) for mol in mols[: self.population_size]]
        for amine in population:
            amine.list_of_options = self.comp_options
            amine.database_path = self.db_mols_path
        return population

    def crossover(self, ind1, ind2):
        mol1 = ind1.mol
        mol2 = ind2.mol
        new_mol = None
        cnt = 0
        while (not new_mol) and cnt <20:
            new_mol = graph_crossover(mol1, mol2)
            cnt +=1
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

    def add_comp_db_attributes(self, pop):
        for amine in pop:
            amine.list_of_options = self.comp_options
            amine.database_path = self.db_mols_path
        return pop

    def run(self):
        results = []  # here the best individuals of each generation will be stored
        self.print_parameters()

        self.population = self.make_initial_population()

        # self.miscs = [get_miscs(options) for options in list_of_options]

        ### save to db.
        ####
        print("Population in run: ", self.population)
        #for pop in self.population:
        #    pop.calculate_score()

        self.population = self.calculate_scores(self.population, gen_id=0)
        #self.calculate_fitness(self.population)

        #children = self.reproduce(self.population, 0 + 1)
        #    #insert_result_to_db(database_path=self.db_mols_path, results=pop.results, list_of_options=self.comp_options)
        
        #self.add_computed_pops_to_db()

        self.db.add_individuals(0, self.population)
        
        self.print_population(self.population, 0)
        
        for n in range(0, self.n_generations):
            print("N-generation: ", n, "\n")
             
            self.calculate_fitness(self.population)
             
            self.db.add_generation(n, self.population)
             
            self.append_results(results, gennum=n, detailed=True)
             
            print(f"population in gen: {n} is {self.population}")

            children = self.reproduce(self.population, n + 1)
                           
            children = self.add_comp_db_attributes(children) #####Appends comp methods and db

            children = self.calculate_scores(children, gen_id=n + 1)
            self.db.add_individuals(n + 1, children)
            self.population = self.prune(self.population + children)
            self.print_population(self.population, n + 1)
        self.calculate_fitness(self.population)
        self.db.add_generation(n + 1, self.population)
        self.append_results(results, gennum=n + 1, detailed=True)
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
    ##Get paths to amines and database.
    database_path = ""
    amines_csv_path  = ""
    
    current_path = os.getcwd()
    if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
        #amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines_with_ions.csv"
        amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/molecules_data.db'
    elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
        #amines_csv_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines_with_ions.csv"
        amines_csv_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/molecules_data.db"
    else:
        print("Path is different than testing or running environemnt")
        print("The path is: ", current_path)

    def canonicalize(smile: str) -> str:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smile))

    amines_csv_path = amines_csv_path[:-10] + "conw_prepped.csv"
    amines = pd.read_csv(amines_csv_path)

    amines['smiles'] = amines['smiles'].apply(canonicalize)


    list_of_options = [{"program":"xtb","method":"gfn_2", "opt":True,  "solvation":"alpb", "solvent":"water"}]#,
    #{"program":"orca","method":"r2SCAN-3c", "solvation":"CPCM", "solvent":"water"}]#,
    
    #list_of_options = [ {"program":"orca","method":"r2SCAN-3c tightopt", "solvation":"CPCM", "solvent":"water"},
    #                    {"program":"orca","method":"CAM-B3LYP Def2-TZVPP",  "solvation":"CPCM", "solvent":"water"}]
                    #]#,
                  #
    ga = GraphGA(
        mol_options=MoleculeOptions(AmineCatalyst),
        population_size=10,
        n_generations=10,
        mutation_rate=0.5,
        db_location="organic.sqlite",
        scoring_kwargs={},
        comp_options=list_of_options,
        db_mols_path = database_path
    )

    res = []
    res = ga.run()

    
    
    #GA output code
    output_file = "ga_10_gen_outputs_007.txt"
    with open(output_file, "w") as f:
        for tp in res:
            try:
                s = tp[0]
            except:
                break
            f.write(f"gen:{tp[0]},n_mols:{len(tp[1])}\n")
            f.write(f"smiles,dH,logk,score\n")
            for ml in tp[1]:
                try:
                    f.write(f"{Chem.MolToSmiles(ml.mol)},{ml.dHabs[1]:.3f},{ml.kabs:0.3f},{ml.score:0.3f}\n")
                except:
                    continue
    print(f"Results: {res}")
    for r in res[-1][1]:
        print("r_mol:", r)
        print(f"Mol: {Chem.MolToSmiles(r.mol)}, dH: {r.dHabs[1]:.3f}, log(k): {r.kabs:0.3f}, score:{r.score:0.3f}")
    generations =[r[0] for r in res]
    best_scores = [max([ind.score for ind in r[1]]) for r in res]
    

    print(f"generations: {generations}, best_scores:{best_scores}")
    fig, ax = plt.subplots()
    ax.plot(generations, best_scores)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Max Score")
    plt.savefig("xXXXx_ga_30_gen_amines_006.png")
