import math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
import copy
import sys
sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/")
sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/catalystGA")
from catalystGA import GA
from catalystGA.reproduction_utils import graph_crossover, graph_mutate
from catalystGA.utils import MoleculeOptions
from xtb import xtb_calculate
from orca import orca_calculate
import sqlite3

#### TODOS:
#### - Remove hardcoded temperature
#### - Intermediate naming for prim/seco/tert amines in the scoring function
#### - Method based initialization of CO2/H2O/OCOO molecules.
#### - dH scoring for A.B ionic compounds. Retrieve reactant data to compute.
#### - Deal with possible duplicates in DB: non-exact but similar namings, two rows with identical methods used but different energies.
#### - Pack Misc values recovery into a single function.


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
        self.program   = "xtb" # orca or xtb 
        self.options   = {} # Specify method , solvation, opt, charge, solvent. Depending on program those will be reshaped.
        self.score     = math.nan
        self.fitness   = math.nan
        self.timing    = math.nan
        self.T_K       = 313 # Kelvin
        self.error     = ""
        self.idx       = (-1, -1)
        self.amine_type = tuple(True if mol.HasSubstructMatch(patt) else False for patt in self.patts)#Respectively primary/secondary/tertiary amine WITHOUT explicit hydrogens.
        self.dHabs = math.nan #Heat of absorbtion
        self.kabs = math.nan #k of reaction limiting step. amine->bicarbonate for tertiary amines
        
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
    
    def prepare_xtb_options(self) -> dict:
        """
        Build an input dicionary for the xtb_calculate function to perfrom xtb calculations.
        """
        xtb_options = {}
        print("optioNS in prepration method:", self.options)
        
        try:
            mtd, tp = self.options["method"].split("_")
            xtb_options[mtd]=int(tp)
        except:
            print("Unspecified QM method")

        try:
            xtb_options[self.options["solvation"]] = self.options["solvent"]
        except:
            print("Unspecified solvation")

        try:
            xtb_options["opt"] = self.options["opt"]
        except:
            print("Unspecified optimization")

        try:
            xtb_options["charge"] = self.options["charge"]
        except:
            print("Unspecified charge")

        return xtb_options #
    
    def prepare_orca_options(self) -> dict:
        orca_options = {}
        options_string = ""
        try:
            #mtd, tp = self.options["method"].split("_") ## Probably unnecessary for orca calculations.
            options_string += (f'{self.options["method"]}').upper()
        except:
            print("Unspecified QM method")

        try:
            options_string += f' {self.options["solvation"]}({self.options["solvent"]})'.upper()
        except:
            print("Unspecified solvation")

        if self.options["opt"]:
            options_string += ' OPT'
        else:
            print("Unspecified optimization")

        orca_options = {options_string:""}
        return orca_options
    
    def calculate_energy(self, n_cores, charge=0):
        ###Computes an energy for a mol object defined by its SMILES/SMARTS string. 
        # The energy is weighted by the contribution of individual conformers.
        self.mol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(self.mol)))

        _ = Chem.rdDistGeom.EmbedMultipleConfs(
                        self.mol,
                        clearConfs=True,
                        #maxAttempts = 10,
                        numConfs=200,
                        useRandomCoords=True,
                        pruneRmsThresh=0.1,
                        #randomSeed=5
                    )
        atoms = [atom.GetSymbol() for atom in self.mol.GetAtoms()]
        
        if self.program =="xtb":
            xtb_options = self.prepare_xtb_options()
            print("self options: ", self.options)
            xtb_options["charge"] = charge
            print("XTB options: ", xtb_options)
            return [xtb_calculate(atoms=atoms, coords=conformer.GetPositions(), options=xtb_options, n_cores=n_cores) for conformer in (self.mol).GetConformers()]
        elif self.program == "orca":
            charge = self.options.pop("charge") #Removes Charge key/value and returns the value to be used in orca_calculate
            orca_options = self.prepare_orca_options()
            return [orca_calculate(atoms=atoms, coords=conformer.GetPositions(), options=orca_options, n_cores=n_cores, charge=charge) for conformer in (self.mol).GetConformers()]
        else:
            raise "Incorrect specification of the QM program."
    
    def weight_energy(self, confs):
        mv = min([v[2] for v in confs])
        boltz_exponents = [(val[2]-mv)/(self.K_B * self.T_K) for val in confs ]
        boltzmann_pop_reactants = [math.exp(-boltz_expon) for boltz_expon in boltz_exponents]
        return sum([reactant_pop*conf_e[2] for reactant_pop, conf_e in zip(boltzmann_pop_reactants, confs)])/sum(boltzmann_pop_reactants)
        
    def cat_products(self, patt, repl, n_cores):
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
        self.mol = Chem.MolFromSmiles(Chem.MolToSmiles(self.mol))
        # Replacement step. It is dependenent on the whether the molecule was sanitized or not.
        print("Check recognition: ", Chem.MolToSmiles(self.mol),Chem.MolToSmarts(patt),Chem.MolToSmiles(repl))
        products = Chem.rdmolops.ReplaceSubstructs(mol=self.mol, query=patt, replacement=repl)

        for prod in products:
            print("Check products ",Chem.MolToSmiles(prod))
            cat = AmineCatalyst(prod)
            cat.options = self.options
            cat.program = self.program
            
            confs = cat.calculate_energy(n_cores=n_cores, charge=1)

            yield [Chem.MolToSmiles(cat.mol), cat.weight_energy(confs)]

    def compute_products(self, n_cores):
        if self.amine_type[0]: # If mol has primary amine
            pri_cats = [ prod for prod in self.cat_products(patt=self.patts[0], repl=self.repls[0], n_cores=n_cores)]
        else:
            pri_cats = []

        if self.amine_type[1]: # If mol has secondary amine
            sec_cats = [ prod for prod in self.cat_products(patt=self.patts[1], repl=self.repls[1], n_cores=n_cores)]
        else:
            sec_cats = []

        if self.amine_type[2]: # If mol has tertiary amine.
            ter_cats = [ prod for prod in self.cat_products(patt=self.patts[2], repl=self.repls[2], n_cores=n_cores)]
        else:
            ter_cats = []

        return pri_cats + sec_cats + ter_cats
        
    def calculate_score(
        self, n_cores: int = 1, envvar_scratch: str = "SCRATCH", scoring_kwargs: dict = {}
    ):
        """Calculate score of molecule, store in self.score.

        Args:
            n_cores (int, optional): Number of cores to use when run on cluster. Defaults to 1.
            envvar_scratch (str, optional): Name of environmental variable pointing to scratch directory. Defaults to 'SCRATCH'.
            scoring_kwargs (dict, optional): Additional keyword agruments parsed to scoring function. Defaults to {}.
        """
        # TODO: implement scoring function
        # this is just a placeholder

        ##Check if in database, if yes -> return db values for that mol instead of computing

        conn = sqlite3.connect('molecules_data.db')
        c = conn.cursor()
        
        
        
        #####
        ##### CHECK DB FOR CO2, H2O, OCOO energies given the computation method. If not found compute them and input to database.
        #####

        CO2_energy  = 0
        H2O_energy  = 0
        OCOO_energy = 0

        method, solvation = self.options['method'], self.options['solvation']

        CO2_smiles, H2O_smiles, OCOO_smiles = "O=C=O", "O", "[O-]C(=O)O"

        print("Method, solvation: ", method, solvation)
        try: 
            c.execute("SELECT smiles, energy FROM miscs WHERE method=? AND solvation=?", (method, solvation))
            miscs_data  = c.fetchall()
            names, es   = [ v[0] for v in miscs_data], [ v[1] for v in miscs_data]
            CO2_energy  = float(es[names.index(CO2_smiles)])
            H2O_energy  = float(es[names.index(H2O_smiles)])
            OCOO_energy = float(es[names.index(OCOO_smiles)])
            print("Miscs energies test: ", CO2_energy,H2O_energy,OCOO_energy)
        except:
            print("Miscallenous molecules not in database. Computing and adding to database now...")
            for name in [CO2_smiles, H2O_smiles, OCOO_smiles]:
                charge = 0
                if name == OCOO_smiles:
                    charge = -1
                
                name_mol         = AmineCatalyst(Chem.MolFromSmiles(name))
                name_mol.options = self.options
                print("name_mol options: ", name_mol.options)
                confs_e          = name_mol.calculate_energy(n_cores=n_cores, charge=charge)
                name_energy      = name_mol.weight_energy(confs_e)
                print("Miscs except: ", name, name_energy)
                params = (name, method, solvation, name_energy)
                c.execute("INSERT INTO miscs VALUES(?,?,?,?)", params)

                #### Temporary ugly block for misc energy values.
                if name == CO2_smiles:
                    CO2_energy = name_energy
                if name == H2O_smiles:
                    H2O_energy = name_energy
                if name == OCOO_smiles:
                    OCOO_energy = name_energy

        #####
        ##### Check DB for reactant energy if not computed -> compute and insert
        #####
        reactant_smiles = Chem.MolToSmiles(self.mol)
        rea_id = None
        reactant_energy = 0
        try: 
            c.execute("SELECT id, energy, product_1_id, product_2_id, product_3_id FROM reactants WHERE smiles=? AND method=? AND solvation=?", (reactant_smiles,method, solvation))
            reacs_data  = c.fetchone() #
            rea_id, reactant_energy, product_1_id, product_2_id, product_3_id = reacs_data
            print("inside try reactant energy: ", reactant_energy)
            print("Successful unpacking of reactant row")

        except:
            reactant_confs = self.calculate_energy(n_cores=n_cores,)
            reactant_energy = self.weight_energy(reactant_confs)
            print("inside except reactant energy: ", reactant_energy)

            c.execute("INSERT INTO reactants(smiles, method, solvation, energy) VALUES(?,?,?,?)",(reactant_smiles, method, solvation, reactant_energy))
            c.execute("SELECT id FROM reactants WHERE smiles=? AND method=? AND solvation=?", (reactant_smiles, method, solvation))
            rea_id = c.fetchone()[0]
            product_1_id, product_2_id, product_3_id = None, None, None
        #### Check DB for product by reactant->product ID's. if empty compute energies. Input energies to products table, 
        prods = [None for _ in range(3)]## Hardcoded number of possible products.

        ### Unpack precomputed information from SQL database.
        for n, prod_id in enumerate([product_1_id, product_2_id, product_3_id]):
            if prod_id is not None:
                c.execute("SELECT * FROM products WHERE id=?", (str(prod_id)))
                pro = c.fetchone()
                print("Pro objec : ", pro)
                # Unpack
                prods[n] = list(pro)
            else:
                continue

        ### Given none products exist already -> compute them.
        prod_ids = []
        print("prods, truth", prods, all(val is None for val in prods))
        if all(val is None for val in prods):


            ### Compute primary/secondary/tertiary amine protonation product.
            amine_products_all = self.compute_products(n_cores=n_cores)
            

            ### Just for information -> If it gets shown often I will need to introduce more product id's
            if len(amine_products_all) > 3:
                print("More than 3 possible protonation sites.")

            ### Insert computed products into database.
            ### Link ids between reactants and products.

            def get_dH (e_prod):
                products = e_prod + OCOO_energy
                reactants = reactant_energy + CO2_energy + H2O_energy
                return abs(products - reactants)
            

            prods = []

            for amine_product, prod_id_col_name in zip(amine_products_all[:3], ['product_1_id','product_2_id','product_3_id']):
                #reactants_energy = reactant_energy + CO2_energy + H2O_energy
                #products_energy = amine_product[1] + OCOO_energy
                #dH = abs(products_energy - reactants_energy)
                print("MISCS TEsts")
                print(smile)
                print(CO2_smiles, CO2_energy)
                print(H2O_smiles, H2O_energy)
                print(OCOO_smiles, OCOO_energy)
                dH = AmineCatalyst.hartree_to_kcalmol(get_dH(amine_product[1]))
                amine_product[0] = Chem.MolToSmiles(Chem.MolFromSmiles(amine_product[0]))
                print( reactant_smiles, " -> ", amine_product[0])
                print( reactant_energy, " -> ", amine_product[1])

                c.execute("INSERT INTO products(smiles,method,solvation, energy, dH, reactant_id) VALUES(?,?,?,?,?,?)", (amine_product[0], method, solvation,amine_product[1], dH, rea_id))
                c.execute("SELECT * FROM products WHERE smiles=? AND method=? AND solvation=?", (amine_product[0], method, solvation))
                prod_am = c.fetchone()
                prod_ids.append(prod_am[0])
                prods.append(prod_am)
                query = f'UPDATE reactants SET {prod_id_col_name}=? WHERE id=?'
                print("Check query: ", query)
                c.execute(query, ([prod_id, rea_id]))
        else:
            prod_ids = [prod[0] for prod in prods]

        prod_ids =[ (pid,) for pid in prod_ids ]
        print("pord_ids, ", len(prod_ids))
        #c.executemany("SELECT dH FROM products WHERE id=?", prod_ids)
        c.execute("""
            SELECT smiles,dH FROM products
            WHERE id IN ({0})""".format(', '.join(str(i[0]) for i in prod_ids)))
        out = c.fetchall()
        
        self.dHabs = max([dH[1] for dH in out])

        print("Product smiles: ",  [val[0] for val in amine_products_all])


        #product_energy = 0 # Compute for each possible product OR weight them by boltzmann
        # print("Reactant energy: ", self.weight_energy(reactant_confs))

        


        
        ### Decide on which product to use by k value:

        #Assign score values based on dH, k, SA

        #dH scorings alone.
        self.score = self.dHabs

        #c.execute("DELETE FROM products")
        #c.execute("DELETE FROM reactants")
        #c.execute("DELETE FROM miscs")


        conn.commit()
        conn.close()
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
    ):
        super().__init__(
            mol_options=mol_options,
            population_size=population_size,
            n_generations=n_generations,
            mutation_rate=mutation_rate,
            db_location=db_location,
            scoring_kwargs=scoring_kwargs,
        )

    def make_initial_population(self):
        with open("data/amines.data", "r") as f:
            lines = f.readlines()
        mols = [Chem.MolFromSmiles(line.strip(","))[0] for line in lines]
        population = [AmineCatalyst(mol) for mol in mols[: self.population_size]]
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
        
        self.population = self.calculate_scores(self.population, gen_id=0)
        
        self.db.add_individuals(0, self.population)
        
        self.print_population(self.population, 0)
        
        for n in range(0, self.n_generations):
            print("N-generation: ", n, "\n")
            self.calculate_fitness(self.population)
            self.db.add_generation(n, self.population)
            self.append_results(results, gennum=n, detailed=True)
            children = self.reproduce(self.population, n + 1)
            children = self.calculate_scores(children, gen_id=n + 1)
            self.db.add_individuals(n + 1, children)
            self.population = self.prune(self.population + children)
            self.print_population(self.population, n + 1)
        self.calculate_fitness(self.population)
        self.db.add_generation(n + 1, self.population)
        self.append_results(results, gennum=n + 1, detailed=True)
        
        return results

if __name__ == "__main__":
    import numpy as np
    from scipy import stats
    #import time
    import pandas as pd 
    import matplotlib.pyplot as plt
    amines = pd.read_csv("examples/data/amines.csv")
    calc_dH, exp_dH = [], []

    cnt = 0
    names, dHs = [],[]

    for smile, dH in zip(amines["SMILES"],amines["dH"]):
        
        if cnt == 30:
             break
        if smile == "CCCCCCCCCCCCNCCO":
            continue
        names.append(smile)
        if "." in smile:
        
            continue


        mol = AmineCatalyst(Chem.MolFromSmiles(smile))
        mol.program = "xtb"
        mol.options = {"method":"gfn_2", "opt":True, "solvation":"alpb", "solvent":"water"}

        mol.calculate_score()
        print("mol.dHabs", mol.dHabs, type(mol.dHabs))
        calc_dH.append(abs(mol.dHabs))
        exp_dH.append(dH)
        print("exp smiles: ", smile)
        print("exp dH", exp_dH)
        print("MEA dH", calc_dH)
        cnt+=1

    plt.scatter(exp_dH, calc_dH, marker="o", color="b")
    plt.xlabel("Experimental " + r"$ \Delta H $")
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
    plt.legend()
    #plt.savefig("start_pop_no_ions_dH_Calc_GFN1_alpb_water.eps", format='eps')
    plt.show()
    plt.close()
    

    #A_cat = AmineCatalyst(m)
    #start = time.time()
    #A_cat.calculate_score()
    #end = time.time()
    #print("Duration: ", end-start)
    #print("Scoring value: ", AmineCatalyst.hartree_to_kcalmol(A_cat.score))

    #population = GA.make_initial_population()





   
    """
    print("Reactant conf:", reactant_confs)
    print("KB T_K: ", reactant_confs[0][2]/(AmineCatalyst.K_B * AmineCatalyst.T_K))
   
    print("Boltzmanns: ", [conf_e[2] for conf_e in reactant_confs])
    reactant_boltzmann_dists = [math.exp(conf_e[2]/(AmineCatalyst.K_B * AmineCatalyst.T_K))  for conf_e in reactant_confs]
    print("reactant boltz dist:", reactant_boltzmann_dists)
    norm_factor= sum(reactant_boltzmann_dists)
    reactant_energy = AmineCatalyst.H2O_energy + AmineCatalyst.CO2_energy + sum([reactant_pop*conf_e[2] for reactant_pop, conf_e in zip(reactant_boltzmann_dists, reactant_confs)])/norm_factor
    print("reactant energy: ", reactant_energy)
    print(reactant_energy)
    """
    """
    import matplotlib.pyplot as plt

    ga = GraphGA(
        mol_options=MoleculeOptions(AmineCatalyst),
        population_size=5,
        n_generations=5,
        mutation_rate=0.5,
        db_location="organic.sqlite",
        scoring_kwargs={},
    )

    results = ga.run()

    generations = [r[0] for r in results]
    best_scores = [max([ind.score for ind in res[1]]) for res in results]

    fig, ax = plt.subplots()
    ax.plot(generations, best_scores)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Max Score")

    plt.savefig("organic.png")
    """