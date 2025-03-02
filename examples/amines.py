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
from xtb import xtb_calculate
from orca import orca_calculate
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
        self.program   = "xtb" # orca or xtb 
        self.options   = {} # Specify method , solvation, opt, solvent. Depending on program those will be reshaped.
        self.score     = math.nan
        self.fitness   = math.nan
        self.timing    = math.nan
        self.T_K       = 313 # Kelvin
        self.error     = ""
        self.idx       = (-1, -1)
        self.amine_type = tuple(True if mol.HasSubstructMatch(patt) else False for patt in self.patts)#Respectively primary/secondary/tertiary amine WITHOUT explicit hydrogens.
        self.dHabs = math.nan #Heat of absorbtion
        self.kabs = math.nan #k of reaction limiting step. amine->bicarbonate for tertiary amines
        self.results = {} #A dictionary with possible keys: "miscs", "reactant", "products", each will be only generated
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

        if self.options["opt"]:
            #options_string += ' OPT'
            orca_options["OPT"] = ""
        else:
            print("Unspecified optimization")
        return orca_options
    
    def calculate_energy(self, n_cores):
        ###Computes an energy for a mol object defined by its SMILES/SMARTS string. 
        # The energy is weighted by the contribution of individual conformers.

        self.mol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(self.mol)))
        print("Name of mol:", Chem.MolToSmiles(self.mol))

        threshold = 0.5
        if  Chem.RemoveHs(self.mol).GetNumAtoms() > 9:
            threshold = 1.0


        _ = Chem.rdDistGeom.EmbedMultipleConfs(
                        self.mol,
                        clearConfs=True,
                        # maxAttempts = 10,
                        numConfs=500,
                        useRandomCoords=True,
                        pruneRmsThresh=threshold,
                        #randomSeed=5
                    )
        
        optimized_confs = AllChem.MMFFOptimizeMoleculeConfs(self.mol, mmffVariant='MMFF94')


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
    
    def weight_energy(self, confs):
        # I pre-divide the exponent to get rid of the extraordinarly large exponents. Since all energies are enforced
        # to be of the similar magnitude.
        boltz_exponents = [((val[2]-confs[0][2]))/(self.K_B * self.T_K) for val in confs ]
        print("BOLTZ exponents", boltz_exponents)
        boltzmann_pop_reactants = [math.exp(-boltz_expon) for boltz_expon in boltz_exponents]
        print("Check normalization: ", boltzmann_pop_reactants, sum(boltzmann_pop_reactants))
        return sum([reactant_pop*conf_e[2] for reactant_pop, conf_e in zip(boltzmann_pop_reactants, confs)])/sum(boltzmann_pop_reactants)
        
    def compute_and_weight_energy(self, n_cores):
        mol_name = Chem.MolToSmiles(self.mol)
        if "." in mol_name:
            submols = mol_name.split(".")
            tot_e = 0

            maxsub = max(submols, key=len) #Picks the largest ionic part of the molecule.
            conf_coords = []

            for submol in submols:
                print("SUBMOL: " , submol)
                sub = AmineCatalyst(Chem.MolFromSmiles(submol))
                sub.program, sub.options = self.program, copy.copy(self.options)
                if Chem.MolFromSmiles(submol).GetNumAtoms() == 1 :
                    print("SUBOPTIONS: ", submol, sub.options)
                    _ = sub.options.pop("opt")
                confs = sub.calculate_energy(n_cores=n_cores)
                if submol == maxsub:
                    conf_coords = confs
                tot_e += sub.weight_energy(confs)
            
            #Pick largest submol -> under assumption the smaller parts are simple ions.
            return tot_e, [conf[1] for conf in conf_coords]
        else:
            confs = self.calculate_energy(n_cores=n_cores)

            return self.weight_energy(confs), [conf[1] for conf in confs]
        
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

        products = [Chem.MolFromSmiles(m) for m in set([Chem.MolToSmiles(p) for p in products])]

        print("products in cat_products: ", [Chem.MolToSmiles(p) for p in products])
        for prod in products:
            print("Check products ",Chem.MolToSmiles(prod))
            cat = AmineCatalyst(prod)
            cat.options = self.options
            cat.program = self.program
            #Yields String, tuple of energy and coords:
            yield [Chem.MolToSmiles(cat.mol), cat.compute_and_weight_energy(n_cores=n_cores)] 

    def compute_amine_products(self, n_cores):
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
        print("Inside compute products: ", pri_cats, sec_cats, ter_cats)
        return pri_cats + sec_cats + ter_cats
    
    @staticmethod
    def chk_conn(conn):
        try:
            conn.cursor()
            return True
        except Exception as ex:
            return False
        
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


        #Precomputation flags. I assume that the molecules are found in the database. 
        # If not the flags will be flipped and the computation performed.
        compute_miscs = False
        compute_reactant = False
        compute_products = False

        compute_dG1 = False
        compute_dG2 = False
        compute_dG3 = False


        results_dict = {}

        CO2_energy  = 0
        H2O_energy  = 0
        OCOO_energy = 0

        method, solvation = self.options['method'], self.options['solvation']
        CO2_smiles, H2O_smiles, OCOO_smiles = "O=C=O", "[H]O[H]", "[H]OC(=O)[O-]" # "O", "O=C([O-])O"
        


        reactant_smiles = Chem.MolToSmiles(self.mol)
        rea_id = None
        reactant_energy = 0

        product_1_id, product_2_id, product_3_id = None, None, None
        prods = []#None for _ in range(3)]## Hardcoded number of possible products.
        prod_ids = [None for _ in range(3)]

        conn = sqlite3.connect(database_path)
        print("Is calculate_score connected to database?", AmineCatalyst.chk_conn(conn))
        c = conn.cursor()
        
        
        
        #####
        ##### CHECK DB FOR CO2, H2O, OCOO energies given the computation method. If not found compute them and input to database.
        #####



        try:
            query = "SELECT smiles, energy FROM miscs WHERE method=? AND solvation=?"
            c.execute(query, (method, solvation))
            miscs_data  = c.fetchall()
            names, es   = [ v[0] for v in miscs_data], [ v[1] for v in miscs_data]
            CO2_energy  = float(es[names.index(CO2_smiles)])
            H2O_energy  = float(es[names.index(H2O_smiles)])
            OCOO_energy = float(es[names.index(OCOO_smiles)])
            print("Miscs energies test: ", CO2_energy,H2O_energy,OCOO_energy)
        except:
            compute_miscs = True

        try: 
            
            c.execute("SELECT id, energy, product_1_id, product_2_id, product_3_id FROM reactants WHERE smiles=? AND method=? AND solvation=?", (reactant_smiles,method, solvation))
            reacs_data  = c.fetchone() #
            rea_id, reactant_energy, product_1_id, product_2_id, product_3_id = reacs_data
            prod_ids = [product_1_id, product_2_id, product_3_id]
            print("Reactant recovery test: ", rea_id, reactant_energy, product_1_id, product_2_id, product_3_id)
            print("Successful unpacking of reactant row")
            ### If the data reading step fails compute the reactant value and insert it to database.

        except:
            compute_reactant = True

        #### After reading from reactants there were no pointers to any products, i.e. they are not there.
        if all(val is None for val in prod_ids):
            compute_products = True
        else:
            for prod_id in [product_1_id, product_2_id, product_3_id]:
                print("PROD ID: " , prod_id)
                if prod_id is not None:
                    query = "SELECT * FROM products WHERE id=?"
                    params = [str(prod_id)]
                    c.execute(query, params)
                    pro = c.fetchone()
                    prods.append(tuple([pro[1], pro[4], pro[5], pro[6], pro[7]])) ## [1] : smiles, [4] : computed_energy, [5:7] dG's
        print("prods after retrieve loop: ", prods)


        #####Computation blocks.
        if compute_miscs:
            print("Miscallenous molecules not in database. Computing and adding to database now...")
            miscs_lst =[]
            for name in [CO2_smiles, H2O_smiles, OCOO_smiles]:
                # charge = 0
                # if name == OCOO_smiles:
                #     charge = -1
                name_mol         = AmineCatalyst(Chem.MolFromSmiles(name))
                name_mol.program = self.program
                name_mol.options = self.options
                print("name_mol options: ", name_mol.options)
                name_energy, name_coords = name_mol.compute_and_weight_energy(n_cores=n_cores)#name_mol.weight_energy(confs_e)

                #### Temporary ugly block for misc energy values.
                if name == CO2_smiles:
                    CO2_energy = name_energy
                if name == H2O_smiles:
                    H2O_energy = name_energy
                if name == OCOO_smiles:
                    OCOO_energy = name_energy
                miscs_lst.append([name_mol, name_energy, name_coords])
            results_dict["miscs"] = miscs_lst
                
        if compute_reactant:
            
            reactant_energy, reactant_coords = self.compute_and_weight_energy(n_cores=n_cores)
            print("inside except reactant energy: ", reactant_energy)
            results_dict["reactant"] = [reactant_energy, reactant_coords]

        if compute_products:

            ### Compute primary/secondary/tertiary amine protonation product.
            ################################################
            amine_products_all = self.compute_amine_products(n_cores=n_cores)
            ################################################
            

            ### Just for information -> If it gets shown often I will need to introduce more product id's
            if len(amine_products_all) > 3:
                print("More than 3 possible protonation sites.")

            for prod in amine_products_all:
                prods.append(prod)
            results_dict["products"] = prods
            
        #### Check that all energies are not None
        ##### Compute dH value for each of the products.

        dHs = []
        for amine_product in prods:
            print("AMINE  PRODUCXT", amine_product)
            amine_product_energy = amine_product[1][0] # 1 element on second index is the list of conformer coordinates.
            if amine_product is None:
                continue
            eles = [ amine_product_energy, OCOO_energy, reactant_energy, CO2_energy, H2O_energy]
            ele_names = [ amine_product[0], OCOO_smiles, Chem.MolToSmiles(self.mol), CO2_smiles, H2O_smiles]

            #####CHECK for correct computation and assigning of values.
            for ele, ele_name in zip(eles, ele_names):
                if ele is math.nan:
                    print("This is NaN: ", ele_name, ele)

            amine_product_name = Chem.MolToSmiles(Chem.MolFromSmiles(amine_product[0]))
            products = amine_product_energy + OCOO_energy
            reactants = reactant_energy + CO2_energy + H2O_energy
            
            dHs.append([amine_product_name, AmineCatalyst.hartree_to_kjmol(abs(products - reactants))])

        ###### Compute k  #######

        #dG = compute_dG(mol, dG)
        










        #### Alternatively could be chosen based on the highest k value.
        self.dHabs = max(dHs, key=lambda x :x[1])
        
        print( reactant_smiles, " -> ", amine_product[0] )
        print( reactant_energy, " -> ", amine_product_energy )

        self.results = results_dict
        print("Results dict: ", self.results)
        #results_list = [[reactant_energy], [ ]for prod in prods]

        # print("Reactant energy: ", self.weight_energy(reactant_confs))

        ### Decide on which product to use by k value:

        #Assign score values based on dH, k, SA

        #dH scorings alone.
        self.score = self.dHabs[1]

        #### Later a reordering code will be added here.

        #order_amine_products(dH, dG) -> top three most reactive. 

        ####

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
        comp_program,
        comp_options,
    ):
        super().__init__(
            mol_options=mol_options,
            population_size=population_size,
            n_generations=n_generations,
            mutation_rate=mutation_rate,
            db_location=db_location,
            scoring_kwargs=scoring_kwargs,
        )
        self.comp_program = comp_program
        self.comp_options = comp_options

    def make_initial_population(self):
        amine_pops_path = amines_csv_path
        with open(amine_pops_path, "r") as f: #data -> csv
            lines = f.readlines()[1:]
        mols = [Chem.MolFromSmiles(line.split(",")[0]) for line in lines]
        population = [AmineCatalyst(mol) for mol in mols[: self.population_size]]
        for amine in population:
            amine.program, amine.options = self.comp_program, self.comp_options
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

        ### save to db.
        ####
        print("Population in run: ", self.population)
        for pop in self.population:
            pop.calculate_score()
        #self.population = self.calculate_scores(self.population, gen_id=0)
        for pop in self.population:
            print(Chem.MolToSmiles(pop.mol))
            print(pop.dHabs)
        
        self.add_computed_pops_to_db()

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

    def add_computed_pops_to_db(self,):

        conn = sqlite3.connect(database_path)
        c = conn.cursor()
        
        print("Did connection to database in ga.add_individuals succeed? :", AmineCatalyst.chk_conn(conn))

        prod_id_col_names = ["product_1_id","product_2_id","product_3_id"]
        method, solvation = self.population[0].options["method"], self.population[0].options["solvation"]
        for individual in self.population:
            if bool(individual.results) == False:
                continue
            if "miscs" in individual.results.keys():

                for misc in individual.results["miscs"]:
                    misc_name, misc_energy= Chem.MolToSmiles(misc[0].mol), misc[1]
                    params = (misc_name, method, solvation, misc_energy)
                    print("PARAMS", params)
                    c.execute("INSERT INTO miscs(smiles, method, solvation, energy) VALUES(?,?,?,?)", params)
 
            if "reactant" in individual.results.keys():

                reactant_smiles, reactant_energy = Chem.MolToSmiles(individual.mol) ,individual.results["reactant"]
                c.execute("INSERT INTO reactants(smiles, method, solvation, energy) VALUES(?,?,?,?)",(reactant_smiles, method, solvation, reactant_energy[0]))
                c.execute("SELECT id FROM reactants WHERE smiles=? AND method=? AND solvation=?", (reactant_smiles, method, solvation))
                rea_id = int(c.fetchone()[0])
                # rea_id, reactant_energy, product_1_id, product_2_id, product_3_id
            else:
                rea_id = None

            if "products" in individual.results.keys():

                assert "reactant" in individual.results.keys(),  "Products cannot be inserted if reactant was not."
                amine_products = individual.results["products"]
                prod_ids = []
                prod_id_col_names = prod_id_col_names[:len(amine_products)]
                for amine_product, prod_id_col_name in zip(amine_products, prod_id_col_names) :
                    c.execute("INSERT INTO products(smiles,method,solvation, energy, dH, reactant_id) VALUES(?,?,?,?,?,?)", (amine_product[0], method, solvation, amine_product[1][0], individual.dHabs[1], rea_id))
                    c.execute("SELECT id FROM products WHERE smiles=? AND method=? AND solvation=?", (amine_product[0], method, solvation))
                    prod_am = int(c.fetchone()[0])
                    prod_ids.append(prod_am)
                    query = f'UPDATE reactants SET {prod_id_col_name}=? WHERE id=?'
                    # print("Check query: ", query)
                    c.execute(query, ([prod_am, rea_id]))

        conn.commit()
        conn.close()
    
    @staticmethod
    def plot_dH_vs_dH(exp_dH, calc_dH, options):
        assert type(options) == dict
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
        plt.title(f"{options['method']} {options['solvation']}")
        plt.legend()
        try:
            figname = "_".join([str(val) for val in options.values()])+".eps"
        except: 
            figname = "" + ".eps"
        plt.savefig(figname, format='eps')
        plt.show()
        plt.close()

if __name__ == "__main__": 
    import numpy as np 
    from scipy import stats 
    #import time
    import pandas as pd 
    import matplotlib.pyplot as plt 

    ##Get paths to amines and database.
    current_path = os.getcwd()
    if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
        amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/molecules_data.db'
    elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
        amines_csv_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines.csv"
        database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/molecules_data.db"
    else:
        print("Path is different than testing or running environemnt")
        print("The path is: ", current_path)



    
    amines = pd.read_csv(amines_csv_path)
    #exp_dH = amines.loc[amines['SMILES'].isin(calc_names)]['dH'].tolist()

    calc_dH, exp_dH = [], []

    cnt = 0
    names, dHs = [],[]

    comp_program = "xtb"
    comp_options = {"method":"gfn_1", "opt":True, "solvation":"gbsa", "solvent":"water"}

    #comp_program = "orca"
    #comp_options = {"method":"r2SCAN-3c", "opt":True, "solvation":"CPCM", "solvent":"water"}


    ga = GraphGA(
        mol_options=MoleculeOptions(AmineCatalyst),
        population_size=1,
        n_generations=1,
        mutation_rate=0.0,
        db_location="organic.sqlite",
        scoring_kwargs={},
        comp_options=comp_options,
        comp_program=comp_program
    )

    # m = AmineCatalyst(Chem.MolFromSmiles("[NH3+]CCO"))#112.34863236070932
    # m.options = comp_options
    # m.program = comp_program
    # m.calculate_energy(n_cores=1)
    
    # print("Computed score: ", Chem.MolToSmiles(m.mol), m.score)
    #print(res)
    results = []
    results = ga.run()

    ##########################################################
    ###Temporary code for benchmarking dH computations.#######
    ##########################################################
    calc_names, calc_dH = [],[]
    for molecule in results:
        print("molecuel: ", Chem.MolToSmiles(molecule.mol))
        calc_names.append(Chem.MolToSmiles(molecule.mol))
        calc_dH.append(molecule.dHabs[1])
    ##########################################################


    ##########################################################
    ################Plotting preparation.####################
    ##########################################################
    exp_dH = amines.loc[amines['SMILES'].isin(calc_names)]['dH'].tolist()
    exp_names = amines.loc[amines['SMILES'].isin(calc_names)]['SMILES'].tolist()

    exp_df = pd.DataFrame({"SMILES":exp_names, "dH_exp":exp_dH})
    calc_df = pd.DataFrame({"SMILES":calc_names, "dH_calc":calc_dH})

    dH_df = pd.merge(exp_df, calc_df, on="SMILES")


    print(dH_df)

    #exp_dH = [ v[2] for v in ]

    #generations = [r[0] for r in results]
    #best_scores = [max([ind.dH for ind in res[1]]) for res in results]
    #calc_dH = [max([ind.dH for ind in res[1]]) for res in results]

    GraphGA.plot_dH_vs_dH(dH_df["dH_exp"], dH_df["dH_calc"], comp_options)

    # fig, ax = plt.subplots()
    # ax.plot(generations, best_scores)
    # ax.set_xlabel("Generation")
    # ax.set_ylabel("Max Score")

    #plt.savefig("organic.png")