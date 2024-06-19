import copy
import os
import sys

import numpy as np

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Chem import rdDistGeom
from rdkit import Chem
import sqlite3
from more_itertools import consecutive_groups
RDLogger.DisableLog("rdApp.*")



#####CHANGE DEPENDING ON IF EXECUTED IN JUPYTER NOTEBOOK OR AS .py

#catalyst_dir = os.path.abspath('')
catalyst_dir = os.path.dirname(__file__)

##### TAKE CARE: USED COMPUTATION AND SOLVATION METHODS ARE HARDCODED IN RATHER THAN USED AS VARIABLES.
##### FRAGMENTS(H2O and CO2) MUST BE PRECOMPUTED AND ALREADY EXIST IN THE DATABASE.

#sys.path.append("/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/catalystGA")
sys.path.append(catalyst_dir)

from make_structures import ConstrainedEmbedMultipleConfsMultipleFrags, connect_cat_2d, frags2bonded, bonded2frags
from k_helper_utils import hartree2kcalmol
from xtb_utils import xtb_optimize
from sql_utils import csv_string_to_opt_coords
from xtb import xtb_calculate
from orca import orca_calculate

#ts_file = os.path.join(catalyst_dir, "input_files/ts3_dummy.sdf")
#ts_dummy = Chem.SDMolSupplier(ts_file, removeHs=False, sanitize=True)[0]

database_path = ""
amines_csv_path  = ""

current_path = os.getcwd()
if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/catalystGA":
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

conn = sqlite3.connect(database_path)
c = conn.cursor()

H2O_smile, CO2_smile, OCOO_smile, method, solvation = "'[H]O[H]'", "'O=C=O'", "'[H]OC(=O)[O-]'","'r2SCAN-3c'", "'CPCM'"# "'gfn_2'", "'gbsa'"

query = f"SELECT smiles, energy, opt_coords FROM miscs WHERE smiles={CO2_smile} AND method={method} AND solvation={solvation}"
print("#terdt", query)
c.execute(query)
CO2_info = c.fetchone()
query = f"SELECT smiles, energy, opt_coords FROM miscs WHERE smiles={H2O_smile} AND method={method} AND solvation={solvation}"

c.execute(query)
H2O_info = c.fetchone()
query = f"SELECT smiles, energy, opt_coords FROM miscs WHERE smiles={OCOO_smile} AND method={method} AND solvation={solvation}"

c.execute(query)
OCOO_info = c.fetchone()

conn.commit()
conn.close()




patts = [Chem.MolFromSmarts("[#7X3;H2]"),Chem.MolFromSmarts("[#7X3;H1]"),Chem.MolFromSmarts("[#7X3;H0]")]
amine_types = {0 : "prim", 1 : "seco", 2 : "tert"}

boltzmann_const_kcal = 1.987204259*10**(-3) #kcal/(mol*K)
boltzmann_const_ev = 8.617333262 * 10 **(-5) # eV/K
planck = 4.135667696 * 10 **(-15)
gas_constant_kcal = 1.9872036 * 10**-3

def write_orca_constr_string(atom_ids):
    out = '\n%geom Constraints\n'
    for atom_id in atom_ids:
        out += '    { C ' + str(atom_id) + ' C }\n'
    out+= '    end\n'
    out+= 'end'
    return out

def write_constr_file(atom_ids, method, filename, **kwargs):
    """
    This function writes a constr.inp file used for input in a constrained xtb optimization.
    atom_ids: which atom_ids to fix during computation
    method:   which method to use in computation
    filename: path to the constr.inp file.
    **kwargs: other computation options
    """

    #### Atoms constrain part
    assert len(atom_ids)>0, "Input list in write_constr_file, atom_ids, is empty!"
    atom_ids = [id+1 for id in sorted(atom_ids)]
    groups = [list(group) for group in consecutive_groups(atom_ids)]
    
    constr_str = "$fix \n    atoms:"
    for group in groups:
        if len(group)>1:
            constr_str += str(min(group))+"-"+str(max(group))+", "
        elif len(group)==1:
            constr_str += str(group[0])+", "
        else:
            pass
    constr_str = constr_str[:-2] +"\n\n"

    #### Method part
    if "gbsagrid" in kwargs:
        constr_str += f"${method} \n"
        for key, value in kwargs.items():
            constr_str += f"    {key}: {value}\n"
        constr_str += "$end"

    with open(filename, "w") as f:
        f.write(constr_str)

def compute_k_rate_from_dG(dG, T):
    """
    dG : energy difference between reactants and transition state [kcal/mol]
    T  : temperature [Kelvin] 
    
    output : ln(k)
    
    Due to possiible very high values of k depending on dG and therefore potential overflow
    I compute the natural logarithm of k instead of the actual value.
    """ 
    #a, b = np.log(planck/(boltzmann_const_ev * T)),  - (dG/(boltzmann_const_kcal*T))
    #print(f"k_rate compute check parts : {a}, {b}")
    return np.log((boltzmann_const_ev * T)/planck) - (dG/(boltzmann_const_kcal*T))


def list_of_dGs(mol, barrier=1, ):
    """
    loop over all amines groups and compute dG for each
    Choose barrier : [1,2,3]
    outputs: a list. each element is a list/dict of : id of affected N, dG, geometry?
    """
    outs = []
    mol =Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    mol = Chem.AddHs(mol)
    ###### Problem here? Outs gets reset each mol? Or pattern match issue?
    for i, patt in enumerate(patts):
        atoms_match = [v[0] for v in list(mol.GetSubstructMatches(patt))]
        for atom_id in atoms_match:
            
            dG, geom_r, geom_ts = compute_dG(mol, atom_id, which_step=barrier, amine_type=amine_types[i])
            outs.append([Chem.RemoveHs(mol), barrier, atom_id, dG, geom_r, geom_ts, amine_types[i]])
    return outs

def compute_dG(cat, atom_id, which_step=1, amine_type="tert"):
    """
    Depending on which type of amine is located at atom_id use different TS core.
    """

    ###Local machine path
    ts_core_file_path =  current_path + "/catalystGA/amine_ts_cores" #"/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/catalystGA/amine_ts_cores"
    cat_core_file_path = current_path + "/catalystGA/amine_rea_cores"

    if amine_type == "tert":
        prot_N_patt = "[#7X4;H0;D4;!+1]"
        ts_file = ts_core_file_path + "/tert_core.sdf"
        cat_file = cat_core_file_path + "/tert_rea_core.sdf"
        #ts_dummy = Chem.MolFromSmarts("[1*].[H]O[H].C(=O)=O")
    elif amine_type   == "prim":
        prot_N_patt =  "[#7X4;H2;D2;!+1]"
        ts_file = ts_core_file_path + "/prim_core.sdf"
        cat_file = cat_core_file_path + "/prim_rea_core.sdf"

    elif amine_type == "seco":
        prot_N_patt =  "[#7X4;H1;D3;!+1]"
        ts_file = ts_core_file_path + "/seco_core.sdf"
        cat_file = cat_core_file_path + "/seco_rea_core.sdf"
    else:
        raise Exception("Incorrect amine_type")
    ts_dummy = Chem.SDMolSupplier(ts_file, removeHs=False, sanitize=True)[0]
    cat_dummy = Chem.SDMolSupplier(cat_file, removeHs=False, sanitize=True)[0]

    De, ts3d_geom, cat3d_geom = ts_scoring(cat, ts_dummy, cat_dummy, atom_id=atom_id, idx=(2, 2), ncpus=6, n_confs=4, cleanup=True, prot_N_patt = prot_N_patt, amine_type=amine_type)#, orca_options={"r2SCAN-3c":"",  "CPCM":"water", "detailed_options":'%cpcm smd true \n   SMDsolvent "water" \n end'})#"OPT":"",

    return De, cat3d_geom, ts3d_geom


def ts_scoring(cat, ts_dummy, cat_dummy, atom_id=0, idx=(0, 0), ncpus=1, n_confs=10, cleanup=False, prot_N_patt = "", amine_type="", orca_options={}):
    """Calculates electronic energy difference in kcal/mol between TS and reactants

    Args:
        cat (rdkit.Mol): Molecule containing one tertiary amine
        n_confs (int, optional): Number of confomers used for embedding. Defaults to 10.
        cleanup (bool, optional): Clean up files after calculation.
                                  Defaults to False, needs to be False to work with submitit.

    Returns:
        Tuple: Contains energy difference, Geom of TS and Geom of Cat


    NOTE for finding amine in the connected mol:
    The ordering of the indices corresponds to the atom ordering
    in the query. For example, the first index is for the atom in this molecule that matches the first atom in the query.
    """
    if list(orca_options.keys()).count("detailed_options") == 1: 
        detailed_options = orca_options.pop("detailed_options")
    else:
        detailed_options = ""
    ###Just to get the smile of the complex:
    if amine_type in  ["prim", "seco"]:
        core_mol = Chem.MolFromSmiles("O=C=O.[1*]")
    elif amine_type == "tert":
        core_mol = Chem.MolFromSmiles("O=C=O.[H]O[H].[1*]")
    else:
        print("Unspecified amine type")


    ts2d = connect_cat_2d(core_mol, cat, atom_id, prot_N_patt = prot_N_patt)

    ### CANONICALIZE
    ts2d = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(ts2d))))
    ts2d = Chem.AddHs(ts2d)

    co2 = list(ts2d.GetSubstructMatch(Chem.MolFromSmarts("O=C=O")))
    n_group = ts2d.GetSubstructMatch( cat)[atom_id]

    atoms2join = ()
    if amine_type in ["prim", "seco"]:
        atoms2join = ((n_group, co2[0]),)
    elif amine_type == "tert":
        h2o = ts2d.GetSubstructMatch(Chem.MolFromSmarts("[H]-[O]-[H]"))
        atoms2join = (
                (n_group, h2o[0]), (h2o[2], co2[2]) )
    else:
        print("Incorrect amine_type in ts_scoring")
        pass

    # Embed TS
    ts3d = ConstrainedEmbedMultipleConfsMultipleFrags(
        mol=ts2d,
        core=ts_dummy,
        numConfs=n_confs,
        pruneRmsThresh=0.1,
        force_constant=1e12,
        atoms2join=atoms2join,
    )
    
    atom_ids = list(ts3d.GetSubstructMatch(ts_dummy))

    h2o_ids = ts3d.GetSubstructMatch(Chem.MolFromSmarts("[H]-O-[H]"))
    
    ### PRUNE NON H2O HYDROGENS. might be superfluous.
    for atom_id in atom_ids:
        
        #if atom_id not in h2o_ids:
        #    print(f'Not in h2o_ids: {atom_id}')
        if (ts3d.GetAtomWithIdx(atom_id).GetSymbol() == "H") and (atom_id not in h2o_ids):
            atom_ids.remove(atom_id)
        if (ts3d.GetAtomWithIdx(atom_id).GetSymbol() == "C") and (atom_id not in co2):
            atom_ids.remove(atom_id)

    filename = catalyst_dir +"/input_files/constr.inp"
    write_constr_file(atom_ids, "gbsa", filename)#, gbsagrid="normal")

    #Calc Energy of TS
    ts3d_energy, ts3d_geom = xtb_optimize(
        ts3d,
        #gbsa="water",
        alpb="water",
        #gbsa="methanol",
        opt_level="crude",
        #name=f"{idx[0]:03d}_{idx[1]:03d}_ts",
        input=os.path.join(catalyst_dir, filename),#
        numThreads=ncpus,
        cleanup=cleanup,
    )
    detailed_options_rea = copy.deepcopy(detailed_options)
    if "OPT" in orca_options.keys():
        # add 
        detailed_options += write_orca_constr_string(atom_ids) 


    if orca_options:
        print(f"orca options: {orca_options}")
        res = orca_calculate(atoms=ts3d_geom['atoms'], coords=ts3d_geom['coords'], options=orca_options, xtra_inp_str=detailed_options ,  n_cores=ncpus )
        print(f"first res: {res}")
        ts3d_energy = res['electronic_energy']  

        print("orca options!!!!!")

    ###########################################
    ###########################################
    ###########################################
    # Compute energy for the amine mol -> (dis)connected case      ? ################
    ## frags i.e misc mols are precomputed. in my code in databse. In here hardcoded.
    ###########################################
    ###########################################
    ###########################################


    cat3d = ConstrainedEmbedMultipleConfsMultipleFrags(
            mol=ts2d,
            core=cat_dummy,
            numConfs=n_confs,
            pruneRmsThresh=0.1,
            force_constant=1e12,
            atoms2join=atoms2join,
         )

    # Embed Catalyst
    #cat3d = copy.deepcopy(cat)
    #cat3d = Chem.AddHs(cat3d)
    #cids = Chem.rdDistGeom.EmbedMultipleConfs(cat3d, numConfs=n_confs, pruneRmsThresh=0.1)

    #if len(cids) == 0:
    #    raise ValueError(f"Could not embed catalyst {Chem.MolToSmiles(Chem.RemoveHs(cat))}")
    
    #Calc Energy of Cat
    cat3d_energy, cat3d_geom = xtb_optimize(
        cat3d,
        alpb="water",
        #gbsa="water",
        opt_level="crude",
        #name=f"{idx[0]:03d}_{idx[1]:03d}_cat",
        numThreads=ncpus,
        cleanup=cleanup,
    )

    if orca_options:
        res = orca_calculate(atoms=cat3d_geom['atoms'], coords=cat3d_geom['coords'], options=orca_options, xtra_inp_str=detailed_options_rea, n_cores=ncpus )
        print(f"second res: {res}")
        cat3d_energy = res['electronic_energy'] 
        print("SECOND orca options!!!!!")


    ####
    #frag_energies = [CO2_info[1]]
    #if amine_type=="tert":
    #    frag_energies.append(H2O_info[1])
    #print(f'Frag energies: {frag_energies}')
    
    # Calculate electronic activation energy
    De = (ts3d_energy - cat3d_energy) * hartree2kcalmol #-"frag_energies"
    print("REACHED END")    
    return De, ts3d_geom, cat3d_geom
