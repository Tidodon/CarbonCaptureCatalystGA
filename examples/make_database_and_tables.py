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
   
#build_database(c)


# # many_reactants = [ ('O','gfn2','gbsa','-5.080122224999'),
# #                   ('O=C=O','gfn2','gbsa','-10.306805408736'),
# #                   ('O','gfn2','alpb','-5.085033168595'),
# #                   ('O=C=O','gfn2','alpb','-10.305467633322'),
# #                   ('O','gfn1','alpb','-5.790291412071'),
# #                   ('O=C=O','gfn1','alpb','-11.543823125467')]

# # c.executemany("INSERT INTO miscs(smiles, method, solvation, energy) VALUES(?,?,?,?)", many_reactants)

# products = [('O=C([O-])O','gfn2','gbsa','-15.215593476193'),
#             ('O=C([O-])O','gfn2','alpb','-15.218909360801'),
#             ('O=C([O-])O','gfn1','alpb','-17.131104157996'),
#             ]
# # c.executemany("INSERT INTO miscs(smiles, method, solvation, energy) VALUES(?,?,?,?)", products)

# c.execute("""UPDATE reactants SET product_1_id=products.id
#         FROM products
#         WHERE reactants.method=products.method 
#         AND reactants.solvation=products.solvation""")

###Code to get column names:
c.execute("SELECT * FROM reactants")
v = c.fetchall()
for val in v:
    print(val)
print("Reactant column names: ", [desc[0] for desc in c.description])

def empty_dbs(connector):
   connector.execute("DELETE FROM products")
   connector.execute("DELETE FROM reactants")
   connector.execute("DELETE FROM miscs")

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