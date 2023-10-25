import sqlite3
from pathlib import Path
current_path = str(Path.cwd())
if current_path == "/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA":
   database_path = '/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/molecules_data.db'

elif current_path == "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA":
    database_path = "/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/molecules_data.db"
else:
   print("Outside predefined working directories.")

conn = sqlite3.connect(database_path)
c = conn.cursor()

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
   

def empty_dbs(cursor, *args, **kwargs):
   ##Possible args in molecules_data.db: reactants, miscs, products.
   for arg in args:
      try: 
         query = f"DELETE FROM {arg}"
         params = []
         if len(kwargs) >0:
            query += " WHERE "
            params = []
            for item in kwargs.items():
               query += f" {item[0]}=? AND"
               params.append(f"{item[1]}")
            query = query[:-4]
         cursor.execute(query, params)
      except:
         print(f"Table with name {arg} does not exist.")
         
def print_table_contents(cursor, *args, **kwargs):
   for arg in args:
      try: 
         query= f"SELECT * FROM {arg}"
         if len(kwargs) >0:
            query += " WHERE "
            params = []
            for item in kwargs.items():
               query += f" {item[0]}=? AND" #{item[1]} AND" 
               params.append(f"{item[1]}")
            query = query[:-4]

         cursor.execute(query, params)
         print(f"{arg} column names: ", [desc[0] for desc in cursor.description])
         v = cursor.fetchall()
         for val in v:
            print(val)
      except:
         print(f"Table {arg} does not exist.")


###Code to get column names:
#build_database(c, name1,name2)
# empty_dbs(c)
print_table_contents(c, "miscs", "reactants", "products", solvation="alpb")#, method="gfn_2", solvation="gbsa")







conn.commit()
conn.close()