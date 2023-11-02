### Functions for interacting with my sql database.

import sqlite3

def check_if_miscs_in_db(cursor, method="", solvation=""):
    query = "SELECT smiles, energy, opt_coords FROM miscs WHERE method=? AND solvation=?"
    cursor.execute(query, (method, solvation))
    miscs_data  = cursor.fetchall()
    if miscs_data:
        return recover_miscs(miscs_data)
    else:
        return False
    
def recover_miscs(miscs_data):
    names, es, coords  = [ v[0] for v in miscs_data], [ v[1] for v in miscs_data], [ v[2] for v in miscs_data]
    miscs_dict = {}
    for name, e in zip(names,es):
        miscs_dict[name] = e
    miscs_dict["opt_conformer_coords_of_"+name] = coords
    return miscs_dict


if __name__ == "__main__":
    pass