### Functions for interacting with my sql database.

import sqlite3
import os 
import csv
import sql_utils
current_path = os.getcwd()

if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
   database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/molecules_data.db'

elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
    database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/molecules_data.db"
else:
   print("Outside predefined working directories.")
    
def check_if_in_db(cursor, smile, method="", solvation="", table=""):
   query = f"SELECT smiles, energy, opt_coords, atoms_ord FROM {table} WHERE smiles=? AND method=? AND solvation=?"
   cursor.execute(query, (smile, method, solvation))
   data  = cursor.fetchone()
   if data:
      return recover_data(data)
   else:
      return False#[]

def recover_data(data):
    name, e, coords, atoms_ord  = data
    coords     = sql_utils.csv_string_to_opt_coords(coords)
    atoms_ord  = sql_utils.csv_string_to_atoms_ord(atoms_ord)
    data_dict  = {}
    #for name, e in zip(names,es):
    data_dict[name] = tuple([e, coords, atoms_ord])
    return data_dict

def build_database(c):
   c.execute("""CREATE TABLE reactants(
        id INTEGER PRIMARY KEY,
       smiles text, 
        method text,
        solvation text,
        energy real,   
         opt_coords text,
         atoms_ord text
   )""")
   c.execute("""CREATE TABLE miscs(
         smiles text, 
         method text,
         solvation text,
         energy real,
         opt_coords text,
         atoms_ord text
   )""")
   c.execute("""CREATE TABLE products(
          id INTEGER PRIMARY KEY,
          smiles text,
          method text,
          solvation text,
          energy real DEFAULT NULL,
         opt_coords text,
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

def empty_dbs(database_path, *args, **kwargs):
   ##Possible args in molecules_data.db: reactants, miscs, products.
   conn = sqlite3.connect(database_path)
   cursor = conn.cursor()

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

   conn.commit()
   conn.close()
   
         
def print_table_contents(database_path, *args, **kwargs):
   conn = sqlite3.connect(database_path)
   cursor = conn.cursor()
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
   conn.commit()
   conn.close()

def opt_coords_to_csv_string(nd_lst):
   """
   Turns a 3d array with floats into a csv string.
   Intended is the turning of the list of conformer xyz coordinates into a saveable csv string.
   """
   out = ""

   for sub_arr in nd_lst:
      for row in sub_arr:
         row  = [str(val) for val in row]
         out += ",".join(row) + ";"

      out  = out[:-1]
      out += "_"
   return out[:-1]

def atoms_ord_to_csv_string(nd_lst):
   out = ""
   for row in nd_lst:
      row  = [str(val) for val in row]
      out += ",".join(row) + "\n"
   return out[:-1]

def csv_string_to_atoms_ord(csv_string):

   reader = csv.reader(csv_string.split('\n'))
   return [val[0] for val in reader]

def csv_string_to_opt_coords(csv_string, func=float):

   out = []

   reader = csv.reader(csv_string.split('_'))

   for subarr in reader:
      subarr  = ",".join(subarr)
      subarr  = subarr.split(";")
      sub_out = []

      for row in subarr:
         row = row.split(",")
         sub_out.append(list(map(func,[val for val in  row])))
      out.append(sub_out)

   return out

def insert_result_to_db(database_path, results, list_of_options):
   conn = sqlite3.connect(database_path)
   cursor = conn.cursor()

   for step, options in zip(results,list_of_options):
      re, pr, mi = step
      insert_mols_e_to_db(cursor, re, method=options["method"], solvation=options["solvation"], table="reactants")
      insert_mols_e_to_db(cursor, pr, method=options["method"], solvation=options["solvation"], table="products" )
      insert_mols_e_to_db(cursor, mi, method=options["method"], solvation=options["solvation"], table="miscs"    )

   conn.commit()
   conn.close()

def insert_mols_e_to_db(database_path, res, method, solvation, table):
   conn = sqlite3.connect(database_path)
   cursor = conn.cursor()

   for mol, data in res.items():
      
      if check_if_in_db(cursor, mol, method=method, solvation=solvation, table=table):
         continue
      else:
         opt_coords = opt_coords_to_csv_string(data[1])
         atoms_ord  = atoms_ord_to_csv_string( data[2])
         params     = (mol, data[0], method, solvation, opt_coords, atoms_ord)
         query      = f"INSERT INTO {table} (smiles, energy, method, solvation, opt_coords, atoms_ord) VALUES(?,?,?,?,?,?)"
         cursor.execute(query, params)

   conn.commit()
   conn.close()


if __name__ == "__main__":
   lst = [[[0,1],[23,-5]], [[3,989],[2322,-100]]]
   conn = sqlite3.connect(database_path)
   c = conn.cursor()

   # c.execute("DROP table miscs")
   # c.execute("DROP table reactants")
   # c.execute("DROP table products")

   # build_database(c)

   # print_table_contents(c, "miscs", "reactants", "products")

   # stringed_list = opt_coords_to_csv_string(lst)
   # print(stringed_list)
   # arred_list = csv_string_to_opt_coords(stringed_list)
   # print(arred_list)

   nd_lst = ["A","C", "h", "C"]
   atoms =  atoms_ord_to_csv_string(nd_lst)
   print(atoms)
   atoms_ord = csv_string_to_atoms_ord(atoms)
   print("atoms_ord:", atoms_ord)


   # arred_string = csv_string_to_arr(stringed_list)
   # print(arred_string)
   conn.commit()
   conn.close()