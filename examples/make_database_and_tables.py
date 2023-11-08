import sqlite3
import os 
import sql_utils
current_path = os.getcwd()

if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
   database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/molecules_data.db'

elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
    database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/molecules_data.db"
else:
   print("Outside predefined working directories.")


def build_database(c):
   c.execute("""CREATE TABLE reactants(
         id INTEGER PRIMARY KEY,
         smiles text, 
         method text,
         solvation text,
         energy real,
         opt_confs text,
         atoms_ord text
   )""")
   c.execute("""CREATE TABLE miscs(
         smiles text, 
         method text,
         solvation text,
         energy real,
         opt_confs text,
         atoms_ord text
   )""")
   c.execute("""CREATE TABLE products(
            id INTEGER PRIMARY KEY,
            smiles text,
            method text,
            solvation text,
            energy real DEFAULT NULL,
            opt_confs text,
            atoms_ord text,
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


###Code to get column names:
#build_database(c, name1,name2)
#empty_dbs(c, "reactants", "products", "miscs")#, method="gfn_1", solvation="gbsa")


#
# "miscs", "reactants", #

if __name__ == "__main__":
   conn = sqlite3.connect(database_path)
   c = conn.cursor()

   # c.execute("SELECT name FROM sqlite_master WHERE type='table';")
   # table_names = [v[0] for v in c.fetchall()]
   # print(table_names)
   #print_table_contents(c, "miscs","products", "reactants" )#"products", solvation="alpb")#, method="gfn_2", solvation="gbsa")
   #empty_dbs(c, "miscs", "products", "reactants")
   import numpy as np
   import csv


   #build_database(c)

   # c.execute("DROP table miscs")
   # c.execute("DROP table reactants")
   # c.execute("DROP table products")
   stringed_list = sql_utils.nd_list_to_csv_string(lst)
   name = "CC"
   params = (name, stringed_list)
   c.execute("INSERT INTO miscs(smiles, opt_confs) values(?,?)",params)
   print_table_contents(c, "miscs","products", "reactants" )

   #c.execute("DELETE FROM miscs WHERE smiles=?",(name,))
   c.execute("SELECT smiles, opt_confs FROM miscs WHERE smiles=?",(name,))
   lstout = c.fetchone()
   print(lstout)

   lstout =  lstout[1]
   print(lstout)

   print(lstout, type(lstout))

   reader = csv.reader(lstout.split('\n'), delimiter=',')
   for row in reader:
      print('\t'.join(row))

   conn.commit()
   conn.close()