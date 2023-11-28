"""
Takes as input:
- argv[1] list_of_smiles path, each row is a reactant smile string to compute dH for
- argv[2:] each arg is an options file, where each row is a "key:value" pair corresponding to computation specifiction.

Output:
- .eps plot of between the input files.

Intended use:
plotting sets of dH data computed using different computational methods against each other.
"""




import os
import sys
import matplotlib.pyplot as plt

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


database_path = ""
amines_csv_path  = ""

current_path = os.getcwd()
if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
    amines_csv_path = "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/data/amines.csv"
    database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/molecules_data.db'
elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
    amines_csv_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines.csv"
    database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/molecules_data.db"
else:
    print("Path is different than testing or running environemnt")
    print("The path is: ", current_path)

import dH_utils

smiles_path      = sys.argv[1]
options_paths     = [path for path in sys.argv[2:]]

print(smiles_path, options_paths)
list_of_smiles        = []
options_list          = [{}]

def check_if_str_is_num(string):
    try:
        val = int(string)
        return val
    except: 
        val = float(string)
        return val
    finally:
        print(f'{string} is not a number.' )

with open(smiles_path, "r") as f:
    for smile in f:
        list_of_smiles.append(smile)

for i,options_path in enumerate(options_paths):
    with open(options_path, "r") as f:
        for line in f:
            line = line.split(":")
            key, val = line[0], line[1]
            options_list[i][key] = check_if_str_is_num(val)

dHs = [[], [] ]
for i, options in enumerate(options_list):
    for smile in list_of_smiles:
        dH_data = dH_utils.compute_dH_data(smile=smile, options=[options], database_path=database_path)
        reactant_energy, product_energies, miscs = dH_data[-1]
        dH = max(map(abs, dH_utils.compute_dH_list( smile, reactant_energy, product_energies, miscs), key=lambda x: x[1]))
        dHs[i].append(dH)

plt.plot(dHs[0], dHs[1])
plt.title(r"$ \Delta H \ vs \ \Delta H$")
plt.xlabel(f'{options_list[0]["method"]} {options_list[0]["solvation"]}')
plt.xlabel(f'{options_list[1]["method"]} {options_list[1]["solvation"]}')
plot_name = f'{options_list[0]["method"]}_{options_list[0]["solvation"]}_vs_{options_list[1]["method"]}_{options_list[1]["solvation"]}.eps'
plt.savefig(plot_name, format='eps')