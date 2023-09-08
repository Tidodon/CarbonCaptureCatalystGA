import math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import sys
sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/")
sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/catalystGA")
from catalystGA import GA
from catalystGA.reproduction_utils import graph_crossover, graph_mutate
from catalystGA.utils import MoleculeOptions
from xtb import xtb_calculate


#### TODOS:
#### - Remove hardcoded temperature
#### - Intermediate naming for prim/seco/tert amines in the scoring function
#### -


class AmineCatalyst:
    save_attributes = {}  # any other attributes to save to the database

    # Reactant energies given by GFN2 method with geom. opt.& GBSA solvation
    CO2_energy = -10.306805753197 # E_h
    H2O_energy = -5.080122224999 # E_h

    #Temperature of the runs:
    T_K = 313 # K
    ### Boltzmann constant
    K_B = 3.166811563 * math.pow(10,-6) # E_h/K

    patts = [Chem.MolFromSmarts("[D1;N;H2]"),Chem.MolFromSmarts("[D2;N;H1]"),Chem.MolFromSmarts("[D3;N;H0]")]
    repls =  [Chem.MolFromSmarts("[NH3+]"),Chem.MolFromSmarts("[NH2+]"),Chem.MolFromSmarts("[NH+]")]

    def __init__(self, mol: Chem.Mol) -> None:
        self.mol = mol
        self.score = math.nan
        self.fitness = math.nan
        self.timing = math.nan
        self.error = ""
        self.idx = (-1, -1)
        self.amine_type = tuple(True if mol.HasSubstructMatch(Chem.MolFromSmarts(patt)) else False for patt in ["[D1;N]","[D2;N]","[D3;N]"])#Respectively primary/secondary/tertiary amine
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
    
    def calculate_energy(self, n_cores, xtb_options={"gfn":2, "opt":True, "gbsa": "water"}):
        ###Computes an energy for a mol object defined by its SMILES/SMARTS string. 
        # The energy is weighted by the contribution of individual conformers.
        self.mol = Chem.AddHs(self.mol)
        _ = Chem.rdDistGeom.EmbedMultipleConfs(
                        self.mol,
                        numConfs=100,
                        useRandomCoords=True,
                        pruneRmsThresh=0.1,
                        randomSeed=3
                    )
        atoms = [atom.GetSymbol() for atom in self.mol.GetAtoms()]
        confs = []
        for conformer in self.mol.GetConformers():
            coords = conformer.GetPositions()
            opt_atoms, opt_coords, opt_energy = xtb_calculate(atoms=atoms, coords=coords, options=xtb_options, n_cores=n_cores
    )
            confs.append([opt_atoms, opt_coords, opt_energy])
        return confs
    
    def weight_energy(self, confs):
        mv = min([v[2] for v in confs])
        boltz_exponents = [(val[2]-mv)/(self.K_B * self.T_K) for val in confs ]
        boltzmann_pop_reactants = [math.exp(-boltz_expon) for boltz_expon in boltz_exponents]
        return sum([reactant_pop*conf_e[2] for reactant_pop, conf_e in zip(boltzmann_pop_reactants, confs)])/sum(boltzmann_pop_reactants)
        
    def cat_products(self, patt, repl, n_cores):

        """
        A generator method that gives smiles representation 
        of the possible product molecules given pattern(patt
        and replacement(repl). It gives energy values for each 
        of the products.

        Arguments:
        patt: recognization pattern given by a mol object
        repl: replacement of the pattern, given by a mol object.
        """

        products = Chem.rdmolops.ReplaceSubstructs(mol=self.mol, query=patt, replacement=repl)
        for prod in products:
            print("TT", Chem.MolToSmiles(prod))
            prod = Chem.AddHs(prod)

            cat = AmineCatalyst(prod)
            print(cat.smiles)
            confs = cat.calculate_energy(n_cores=n_cores)
            yield cat.smiles ,cat.weight_energy(confs)
    

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

        ##Reactant prepare:
        reactant_confs = self.calculate_energy(n_cores=n_cores, )
        reactant_energy = self.weight_energy(reactant_confs)#self.H2O_energy + self.CO2_energy + self.weight_energy(reactant_confs)
        
        product_energy = 0 # Compute for each possible product OR weight them by boltzmann

        

        pri_cats = [ prod for prod in self.cat_products(patt=self.patts[0], repl=self.repls[0], n_cores=n_cores) if self.amine_type[0]] ### List comprehension solution here?
        sec_cats = [ prod for prod in self.cat_products(patt=self.patts[1], repl=self.repls[1], n_cores=n_cores) if self.amine_type[1]]
        ter_cats = [ prod for prod in self.cat_products(patt=self.patts[2], repl=self.repls[2], n_cores=n_cores) if self.amine_type[2]]
        amine_products_all = pri_cats + sec_cats + ter_cats
        print(amine_products_all)
        
        ### Decide on which product to use by k value:


        #Assign score values based on dH, k, SA
        print("RE: ", reactant_energy)
        print("Amines: ", min([val[1] for val in amine_products_all]))
        self.score = min([val[1] for val in amine_products_all]) - reactant_energy 



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
    #import time
    m = Chem.MolFromSmiles("NCCO")
    test = AmineCatalyst(m)
    test.calculate_score()
    print("same?", test.score)
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
    