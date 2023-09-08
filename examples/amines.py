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


class AmineCatalyst:
    save_attributes = {}  # any other attributes to save to the database

    # Reactant energies given by GFN2 method with geom. opt.
    CO2_energy = -10.308452272882 # E_h
    H2O_energy = -5.0705442454 # E_h

    #Temperature of the runs:
    T_K = 313 #K
    ### Boltzmann constant
    K_B = 3.166811563 * pow(10,-6) # E_h/K

    def __init__(self, mol: Chem.Mol) -> None:
        self.mol = mol
        self.score = math.nan
        self.fitness = math.nan
        self.timing = math.nan
        self.error = ""
        self.idx = (-1, -1)
        self.amine_type = tuple(True if mol.HasSubstructMatch(Chem.MolFromSmarts(patt)) else False for patt in ["[D1;N]","[D2;N]","[D3;N]"])
        self.dHabs = math.nan #Heat of absorbtion
        self.kabs = math.nan #k of reaction limiting step. amine->bicarbonate for tertiary amines
        # 
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
    
    def calculate_energy(self, n_cores, xtb_options={"gfn":2, "opt":True}):
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
        reactant_energy = self.H2O_energy + self.CO2_energy + self.weight_energy(reactant_confs)
        
        product_energy = 0 # Compute for each possible product OR weight them by boltzmann

        ### Decide on which product to use by k value:

        prim_amines = []
        seco_amines = []
        tert_amines = []
        """
        if self.amine_type[0]:
            patt_prim = Chem.MolFromSmarts("[D1;N]")
            repl_prim = Chem.MolFromSmarts("[N+]")
            outs = Chem.rdmolops.ReplaceSubstructs(mol=m, query=patt_prim, replacement=repl_prim)
            for struct in outs:
                ###Compute k
                ###Compute dH
                prim_mol = Chem.MolFromSmiles(struct)
                Chem.AddHs(tmol2)
                AllChem.EmbedMolecule(m2)
                AllChem.MMFFOptimizeMolecule(m2)
                AllChem.Compute2DCoords(tmol2)
                Chem.MolToMolBlock(tmol2)
        """
        self.score = reactant_energy



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
        with open("data/amines.smi", "r") as f:
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
    m = Chem.MolFromSmiles("CCCC")
    A_cat = AmineCatalyst(m)
    A_cat.calculate_score()
    print("Scoring value: ", A_cat.score)

    m = Chem.AddHs(m)



   
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
    