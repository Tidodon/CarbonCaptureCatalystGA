import sqlite3


conn = sqlite3.connect('molecules_data.db')
c = conn.cursor()

#c.execute(""" CREATE TABLE reactants(
#          smiles text, 
#          method text,
#          solvation text,
#          energy real
#)""")


#many_reactants = [ ('O','gfn2','gbsa','-5.080122224999'),
#                  ('O=C=O','gfn2','gbsa','-10.306805408736'),
#                  ('O','gfn2','alpb','-5.085033168595'),
#                  ('O=C=O','gfn2','alpb','-10.305467633322'),
#                  ('O','gfn1','alpb','-5.790291412071'),
#                  ('O=C=O','gfn1','alpb','-11.543823125467')]
# c.executemany("INSERT INTO reactants VALUES(?,?,?,?)", many_reactants)

#c.execute("""CREATE TABLE products (
#          smiles text,
#          method text,
##          solvation text,
 #         energy real
 #         )""")
#conn.commit()
#products = [('O=C([O-])O','gfn2','gbsa','-15.215593476193'),
#            ('O=C([O-])O','gfn2','alpb','-15.218909360801'),
#            ('O=C([O-])O','gfn1','alpb','-17.131104157996'),
#            ]
#c.executemany("INSERT INTO products VALUES(?,?,?,?)", products)
c.execute("SELECT rowid,* FROM reactants")
res = c.fetchall()

print(res)
conn.commit()
conn.close()
