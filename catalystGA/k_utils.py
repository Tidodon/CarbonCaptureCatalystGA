import copy
import os
import sys

import numpy as np

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

import sqlite3

RDLogger.DisableLog("rdApp.*")


catalyst_dir = os.path.dirname(__file__)
sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/")
sys.path.append(catalyst_dir)

from make_structures import ConstrainedEmbedMultipleConfsMultipleFrags, connect_cat_2d
from k_helper_utils import hartree2kcalmol
from xtb_utils import xtb_optimize

ts_file = os.path.join(catalyst_dir, "input_files/ts3_dummy.sdf")
ts_dummy = Chem.SDMolSupplier(ts_file, removeHs=False, sanitize=True)[0]

frag_energies = np.sum([-8.232710038092, -19.734652802142, -32.543971411432])  # 34 atoms

"""
TODO:
-Get fragment energies from db.

"""



patts = [Chem.MolFromSmarts("[ND1]"),Chem.MolFromSmarts("[ND2]"),Chem.MolFromSmarts("[ND3]")]


amine_types = {0 : "prim", 2 : "seco", 3 : "tert"}


def compute_k_rate_from_dG(dG, T):
    pass


def list_of_dGs(mol, barrier=1, ):
    """
    loop over all amines groups and compute dG for each
    Choose barrier : [1,2,3]
    outputs: a list. each element is a list/dict of : id of affected N, dG, geometry?
    """
    mol = Chem.AddHs(mol)

    outs = []

    for i, patt in enumerate(patts):
        atoms_match = [v[0] for v in list(mol.GetSubstructMatches(patt))]
        for atom_id in atoms_match:
            amine_type = amine_types[i]
            dG, geom_r, geom_ts = compute_dG(mol, atom_id, which_step=barrier, amine_type=amine_type)

            outs.append([mol, barrier, atom_id, dG, geom_r, geom_ts])

        # products = Chem.rdmolops.ReplaceSubstructs(mol=mol, query=patt, replacement=repl)

        # products = [Chem.MolFromSmiles(m) for m in set([Chem.MolToSmiles(p) for p in products])]

    return outs

def compute_dG(smile, atom_id, which_step=1, amine_type="tert"):
    """
    Depending on 
    """

    cat = Chem.MolFromSmiles(smile)

    if amine_type == "tert":
        prot_N_patt = "[#7X4;H0;D4;!+1]"
        ts_dummy = Chem.MolFromSmarts("[1*][H]O[H].C(=O)=O")
        De, ts3d_geom, cat3d_geom = ts_scoring(cat, ts_dummy, idx=(0, 0), ncpus=1, n_confs=10, cleanup=False, prot_N_patt = prot_N_patt)

    elif amine_type == "prim" or amine_type == "seco":
        """
        Add here step dependent change of ts_dummy and cat.
        """
        if amine_type   == "prim":
            prot_N_patt =  "[#7X4;H2;D2;!+1]"
        elif amine_type == "seco":
            prot_N_patt =  "[#7X4;H1;D3;!+1]"
        ts_dummy = Chem.MolFromSmarts("[1*]C(=O)=O")
        De, ts3d_geom, cat3d_geom = ts_scoring(cat, ts_dummy, idx=(0, 0), ncpus=1, n_confs=10, cleanup=False, prot_N_patt = prot_N_patt)

    else:
        raise Exception("Incorrect amine_type")

    return De, cat3d_geom, ts3d_geom



def ts_scoring(cat, ts_dummy, atom_id=0, idx=(0, 0), ncpus=1, n_confs=10, cleanup=False):
    """Calculates electronic energy difference in kcal/mol between TS and reactants

    Args:
        cat (rdkit.Mol): Molecule containing one tertiary amine
        n_confs (int, optional): Number of confomers used for embedding. Defaults to 10.
        cleanup (bool, optional): Clean up files after calculation.
                                  Defaults to False, needs to be False to work with submitit.

    Returns:
        Tuple: Contains energy difference, Geom of TS and Geom of Cat
    """


    ###########################################
    ###########################################
    ###########################################
    ########### Choose ts_dummies for  misc molecules in prim/seco and tert amine cases.
    ###########################################
    ###########################################
    ###########################################
    ts2ds = connect_cat_2d(ts_dummy, cat, atom_id, prot_N_patt = "")
    if len(ts2ds) > 1:
        print(f"{Chem.MolToSmiles(Chem.RemoveHs(cat))} contains more than one tertiary amine")
    ts2d = ts2ds[0]

    # Embed TS
    ts3d = ConstrainedEmbedMultipleConfsMultipleFrags(
        mol=ts2d,
        core=ts_dummy,
        numConfs=n_confs,
        pruneRmsThresh=0.1,
        force_constant=1e12,
    )

    # Calc Energy of TS
    ts3d_energy, ts3d_geom = xtb_optimize(
        ts3d,
        #alpb="water"
        gbsa="methanol",
        opt_level="tight",
        name=f"{idx[0]:03d}_{idx[1]:03d}_ts",
        input=os.path.join(catalyst_dir, "input_files/constr.inp"),
        numThreads=ncpus,
        cleanup=cleanup,
    )

    ###########################################
    ###########################################
    ###########################################
    # Compute energy for the amine mol -> (dis)connected case      ? ################
    ## frags i.e misc mols are precomputed. in my code in databse. In here hardcoded.
    ###########################################
    ###########################################
    ###########################################

    # Embed Catalyst
    cat3d = copy.deepcopy(cat)
    cat3d = Chem.AddHs(cat3d)
    cids = Chem.rdDistGeom.EmbedMultipleConfs(cat3d, numConfs=n_confs, pruneRmsThresh=0.1)
    if len(cids) == 0:
        raise ValueError(f"Could not embed catalyst {Chem.MolToSmiles(Chem.RemoveHs(cat))}")

    # Calc Energy of Cat
    cat3d_energy, cat3d_geom = xtb_optimize(
        cat3d,
        #alpb="water"
        gbsa="methanol",
        opt_level="tight",
        name=f"{idx[0]:03d}_{idx[1]:03d}_cat",
        numThreads=ncpus,
        cleanup=cleanup,
    )

    ####

    # Calculate electronic activation energy
    De = (ts3d_energy - frag_energies - cat3d_energy) * hartree2kcalmol
    return De, ts3d_geom, cat3d_geom


if __name__ == "__main__":
    # cat = Chem.MolFromSmiles("CN(C)C")
    # n_confs = 10

    # ts2ds = connect_cat_2d(ts_dummy, cat)
    # if len(ts2ds) > 1:
    #     print(f"{Chem.MolToSmiles(Chem.RemoveHs(cat))} contains more than one tertiary amine")
    # ts2d = ts2ds[0]

    # print("dummy", Chem.MolToSmiles(ts_dummy))

    # ts3d = ConstrainedEmbedMultipleConfsMultipleFrags(
    #     mol=ts2d,
    #     core=ts_dummy,
    #     numConfs=n_confs,
    #     pruneRmsThresh=0.1,
    #     force_constant=1e12,
    # )
    # print(Chem.MolToSmiles(ts3d))
    
    mol = Chem.MolFromSmiles("CN(C)COCF")
    print(mol)
    #prot_N_patt = Chem.MolFromSmarts("[#7X4;H0;D4;!+1]")
    mol_with_dummy = ts_dummy#Chem.MolFromSmarts("[1*][H]O[H].C(=O)=O")
    print(Chem.MolToSmiles(connect_cat_2d(mol_with_dummy, mol)[0]))#, 1, prot_N_patt = "")))
    #print(mol.GetSubstructMatch(prot_N_patt))