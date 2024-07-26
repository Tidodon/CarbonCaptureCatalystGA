### Functions for interacting with my sql database.

import sqlite3
import os 
import csv
import sql_utils as sql_utils

def check_if_in_db(database_path, smile, method="", solvation="", table=""):
   conn = sqlite3.connect(database_path)
   cursor = conn.cursor()
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
      smiles TEXT, 
      method TEXT,
      solvation TEXT,
      energy REAL,   
      opt_coords TEXT,
      atoms_ord TEXT,
      opt INTEGER  DEFAULT 0 NOT NULL,
      prev_step TEXT
   )""")
   c.execute("""CREATE TABLE miscs(
         smiles TEXT, 
         method TEXT,
         solvation TEXT,
         energy REAL,
         opt_coords TEXT,
         atoms_ord TEXT,
         opt INTEGER DEFAULT 0 NOT NULL,
         prev_step TEXT
   )""")
   c.execute("""CREATE TABLE products(
          id INTEGER PRIMARY KEY,
          smiles TEXT,
          method TEXT,
          solvation TEXT,
          energy REAL DEFAULT NULL,
          opt_coords TEXT,
          atoms_ord TEXT,
          opt INTEGER DEFAULT 0 NOT NULL,
          dH REAL DEFAULT NULL,
          dG1 REAL DEFAULT NULL,
          dG2 REAL DEFAULT NULL,
          dG3 REAL DEFAULT NULL,
          prev_step TEXT 
         )""")
   #  reactant_id integer REFERENCES reactants(rowid) ON UPDATE CASCADE
   # c.execute("ALTER TABLE reactants ADD product_1_id INTEGER DEFAULT NULL REFERENCES products(id) ON UPDATE CASCADE")
   # c.execute("ALTER TABLE reactants ADD product_2_id INTEGER DEFAULT NULL REFERENCES products(id) ON UPDATE CASCADE")
   # c.execute("ALTER TABLE reactants ADD product_3_id INTEGER DEFAULT NULL REFERENCES products(id) ON UPDATE CASCADE")
   # c.execute("ALTER TABLE reactants ADD comments text")
   # c.execute("ALTER TABLE products ADD comments text")

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
   
         
def print_table_contents(database_path, *args, print_data=True, **kwargs):
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
      if print_data:
         for val in v:
            print(val)
      else:
         return v
   conn.commit()
   conn.close()

def convert_table_output_into_dict():
   pass

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

   return [val for val in reader]

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

   steps = len(list_of_options)

   for step, options, ind in zip(results, list_of_options, range(steps)):
      re, pr, mi = step

      opt = 0
      if "opt" in options:
         opt = 1
      prev_step_string = ""
      if ind>0:
         for i in range(ind):
            prev_step_string += options_dict_to_string(list_of_options[i]) +"\n"
      prev_step_string = prev_step_string[:-1]
      insert_mols_e_to_db(database_path, re, method=options["method"], solvation=options["solvation"], opt=opt, table="reactants", prev_step_string=prev_step_string)
      insert_mols_e_to_db(database_path, pr, method=options["method"], solvation=options["solvation"], opt=opt, table="products" , prev_step_string=prev_step_string)
      insert_mols_e_to_db(database_path, mi, method=options["method"], solvation=options["solvation"], opt=opt, table="miscs"    , prev_step_string=prev_step_string)

def options_dict_to_string(options):
   options_string = ""
   for key, value in options.items():
      options_string += str(key) + "_" + str(value) + "_"
   return options_string[:-1]

def insert_mols_e_to_db(database_path, res, method, solvation, opt, table, prev_step_string):
   conn = sqlite3.connect(database_path)
   cursor = conn.cursor()

   for mol, data in res.items():
      
      if check_if_in_db(database_path, mol, method=method, solvation=solvation, table=table):
         continue
      else:
         opt_coords = opt_coords_to_csv_string(data[1])
         atoms_ord  = atoms_ord_to_csv_string( data[2])
         params     = (mol, data[0], method, solvation, opt_coords, atoms_ord, opt, prev_step_string)
         query      = f"INSERT INTO {table} (smiles, energy, method, solvation, opt_coords, atoms_ord, opt, prev_step) VALUES(?,?,?,?,?,?,?,?)"
         cursor.execute(query, params)

   conn.commit()
   conn.close()

if __name__ == "__main__":
   lst = [[[0,1],[23,-5]], [[3,989],[2322,-100]]]
  
   import pandas as pd
 
   current_path = os.getcwd()

   database_path = ""
   amines_csv_path  = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/data/amines.csv"# ""

   database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/molecules_data.db'

   if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
      database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/molecules_data.db'

   elif current_path == "/lustre/hpc/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
      database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/molecules_data.db"
   else:
      print("Outside predefined working directories.")


   database_path = '/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/molecules_data.db' 
   conn = sqlite3.connect(database_path)
   c = conn.cursor()

   # c.execute("DROP table miscs")
   # c.execute("DROP table reactants")
   # c.execute("DROP table products")

   # build_database(c)

   #print_table_contents(database_path, "products",   method="CAM-B3LYP def2-TZVPP")#, method="gfn_1")#solvation="gbsa")

   # stringed_list = opt_coords_to_csv_string(lst)
   # print(stringed_list)
   # arred_list = csv_string_to_opt_coords(stringed_list)
   # print(arred_list)
   amines = pd.read_csv(amines_csv_path)

   nd_lst = [["A","C", "h", "C"]]
   atoms  =  atoms_ord_to_csv_string(nd_lst)
   print(atoms)
   atoms_ord = csv_string_to_atoms_ord(atoms)
   print("atoms_ord:", atoms_ord)
   #empty_dbs(database_path, "reactants", "products", method="gfn_2", solvation="alpb")#)solvation="CPCM")#method="CAM-B3LYP def2-TZVP")
   print_table_contents(database_path, "miscs", "reactants", "products", print_data=True, method="r2SCAN-3c")#, solvation="CPCM")


   
   #dH_df = pd.merge(calc_df, exp_df, on="SMILES")
 
   r2_string= "0.64297,0.380819,-0.070442;1.144848,-0.923631,-0.207996;2.135425,-1.074344,-0.963625;0.522154,-1.78222,0.456882;1.192306,0.955851,-0.623052"#, 'O,C,O,O,H' 
   arred_string = csv_string_to_opt_coords(r2_string)
   print(arred_string)
   # print(arred_string)
   conn.commit()
   conn.close()
