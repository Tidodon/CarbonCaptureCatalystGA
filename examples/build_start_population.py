import numpy as np

outs = []

with open('../AminesFromFigure1.csv', "r") as f:
    names = f.readline()
    k2s = f.readline()
    dHs = f.readline()
    pKas = f.readline()

    print(names, k2s, dHs, pKas)
    names = names.split(";")[1:]
    k2s = k2s.split(";")[1:]
    dHs = dHs.split(";")[1:]
    pKas = pKas.split(";")[1:]



    for name, k2, dH, pKa in zip(names, k2s, dHs, pKas):
        name = name.replace( "*", "")
        name = name.replace( "\n", "")
        if name[-1] == ")":
            name = name[:name.rfind("(")]

        k2 = k2.replace(",",".")
        dH = dH.replace(",",".")
        pKa = pKa.replace(",",".")


        outs.append([name, k2, dH, pKa])


from urllib.request import urlopen
from urllib.parse import quote

def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'
    
smiles_list = [CIRconvert(val[0]) for val in outs ]
k2_list = [val[1] for val in outs]
dH_list = [val[2] for val in outs]
pKa_list = [val[3] for val in outs]

smiles_list[5] =  "C"*12 + smiles_list[0]

with open("amines.data", "w") as f:
    for smile, k2, dH, pKa in zip(smiles_list, k2_list, dH_list, pKa_list):
        f.write(smile+","+k2+","+dH+","+pKa+"\n")
