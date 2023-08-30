from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import numpy as np

smiles = []
with open("amines.smi", "r") as f:
    for smi in f:
        smiles.append(smi.replace("\n", ""))

mols = [Chem.MolFromSmiles(smi) for smi in smiles]

fpgen = AllChem.GetRDKitFPGenerator()
fps = [fpgen.GetFingerprint(m) for m in mols]
dim = len(fps)
similarity_matrix = [ [DataStructs.DiceSimilarity(fps[i], fps[j]) for i in range(dim)] for j in range(dim) ]

def max_entropy(lst):
    s = sum(lst)
    return -np.log10(max([ val/s for val in lst]))

entropies = [max_entropy(lst) for lst in similarity_matrix]

mols = [mols[i] for i in sorted(range(len(entropies)), key=entropies.__getitem__)]
entropies = [entropies[i] for i in sorted(range(len(entropies)), key=entropies.__getitem__)]

for i in range(dim):
    print(Chem.MolToSmiles(mols[i]), entropies[i])


###Maxime entropy across different similarity measures
