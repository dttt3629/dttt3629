from multiprocessing import freeze_support

from rdkit import Chem

from mordred import Calculator, descriptors

from rdkit.Chem import Descriptors
#from mordred 文章注明
if __name__ == "__main__":
    freeze_support()

    active = Chem.SDMolSupplier("decoy_f.sdf")
    # Create Calculator
    calc = Calculator(descriptors)

    # map method calculate multiple molecules (return generator)
    # print(list(calc.map(mols)))

    # pandas method calculate multiple molecules (return pandas DataFrame)
    # print(calc.pandas(mols))
    renew = []
    for sm in active:
        try:
            c = Descriptors.fr_sulfone(sm)
            if c != -666:
                renew.append(sm)
    #删掉错误分子
        except:
            next
    c = calc.pandas(renew)
    c.to_csv("decoy_dpt.csv")

   
