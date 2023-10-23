import sqlite3


#conn = sqlite3.connect('/Users/dbo/Documents/CarbonCapture/GA_playground/CarbonCaptureCatalystGA/examples/molecules_data.db')
conn = sqlite3.connect('/groups/kemi/orlowski/CarbonCapture/CarbonCaptureCatalystGA/examples/molecules_data.db')
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
   

def empty_dbs(cursor, *args):
   ##Possible args in molecules_data.db: reactants, miscs, products.
   for arg in args:
      try: 
         query = f"DELETE FROM {arg}"
         cursor.execute(query)
      except:
         print(f"Table with name {arg} does not exist.")
         
def print_table_contents(cursor, *args):
   for arg in args:
      try: 
         query= f"SELECT * FROM {arg}"
         cursor.execute(query)
         print(f"{arg} column names: ", [desc[0] for desc in cursor.description])
         v = cursor.fetchall()
         for val in v:
            print(val)
      except:
         print(f"Table {arg} does not exist.")


###Code to get column names:
#build_database(c, name1,name2)
# empty_dbs(c)




c.execute("SELECT * FROM miscs")
print("Miscs column names: ", [desc[0] for desc in c.description])
v = c.fetchall() # TUPLE object, where each element corresponds to column
for val in v:
   print(val)

c.execute("SELECT * FROM products")
print("Products column names: ", [desc[0] for desc in c.description])
v = c.fetchall() # TUPLE object, where each element corresponds to column
for val in v:
    print(val)
# for val in v:
#    print(val)

# print("Product column names: ", [desc[0] for desc in c.description])
#c.execute("DROP TABLE products")







conn.commit()
conn.close()