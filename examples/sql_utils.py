### Functions for interacting with my sql database.

import sqlite3
import os 
current_path = os.getcwd()

if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
   database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/molecules_data.db'

elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
    database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/molecules_data.db"
else:
   print("Outside predefined working directories.")


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
         miscs_dict[name] = tuple([e, coords])
    return miscs_dict

def build_database(c):
   c.execute("""CREATE TABLE reactants(
        id INTEGER PRIMARY KEY,
       smiles text, 
        method text,
        solvation text,
        energy real
   )""")
   c.execute("""CREATE TABLE miscs(
        smiles text, 
        method text,
        solvation text,
        energy real
   )""")
   c.execute("""CREATE TABLE products(
          id INTEGER PRIMARY KEY,
          smiles text,
          method text,
          solvation text,
          energy real DEFAULT NULL,
          dH real DEFAULT NULL,
          k1 real DEFAULT NULL,
          k2 real DEFAULT NULL,
          k3 real DEFAULT NULL,
          reactant_id integer REFERENCES reactants(rowid) ON UPDATE CASCADE
         )""")
   c.execute("ALTER TABLE reactants ADD product_1_id INTEGER DEFAULT NULL REFERENCES products(id) ON UPDATE CASCADE")
   c.execute("ALTER TABLE reactants ADD product_2_id INTEGER DEFAULT NULL REFERENCES products(id) ON UPDATE CASCADE")
   c.execute("ALTER TABLE reactants ADD product_3_id INTEGER DEFAULT NULL REFERENCES products(id) ON UPDATE CASCADE")
   c.execute("ALTER TABLE reactants ADD comments text")
   c.execute("ALTER TABLE products ADD comments text")
   c.execute("ALTER TABLE reactants ADD opt_coords BLOB DEFAULT NULL")
   c.execute("ALTER TABLE miscs ADD opt_coords BLOB DEFAULT NULL")
   c.execute("ALTER TABLE products ADD opt_coords BLOB DEFAULT NULL")


def empty_dbs(cursor, *args, **kwargs):
   ##Possible args in molecules_data.db: reactants, miscs, products.
   for arg in args:
      query = f"DELETE FROM {arg}"
      params = []
      if len(kwargs) >0:
         query += " WHERE "
         for item in kwargs.items():
            query += f" {item[0]}=? AND"
            params.append(f"{item[1]}")
         query = query[:-4]
      cursor.execute(query, params)
         
def print_table_contents(cursor, *args, **kwargs):
   for arg in args:
      query= f"SELECT * FROM {arg}"
      params = []
      if len(kwargs) >0:
         query += " WHERE "
         for item in kwargs.items():
            query += f" {item[0]}=? AND" #{item[1]} AND" 
            params.append(f"{item[1]}")
         query = query[:-4]
      cursor.execute(query, params)
      print(f"{arg} column names: ", [desc[0] for desc in cursor.description])
      v = cursor.fetchall()
      for val in v:
         print(val)
if __name__ == "__main__":
    pass