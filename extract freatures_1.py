from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem
active = Chem.SDMolSupplier('D:\\work\\actives_final.mol2\\decoy_f.sdf')
des_list = [x[0] for x in Descriptors._descList]    
renew = []
for sm in active:
    try:
        c = Descriptors.fr_sulfone(sm)
        if c != -666:
            renew.append(sm)
    #删掉错误分子
    except:
        next
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(des_list)
with open("D:\\work\\actives_final.mol2\\decoy_fdpt.csv",mode = "w") as dpt: 
    for x in des_list:
        dpt.write(str(x)+",")
    dpt.write("\n")
    for sm in renew :
        try:
            c = calculator.CalcDescriptors(sm)
            for x in c:
                dpt.write(str(x)+",")
        except:
            next
        dpt.write("\n")