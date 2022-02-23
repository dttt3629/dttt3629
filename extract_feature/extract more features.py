from multiprocessing import freeze_support

from rdkit import Chem

from mordred import Calculator, descriptors
#from mordred 文章注明
if __name__ == "__main__":
    freeze_support()

    active = Chem.SDMolSupplier("active_f.sdf")
    # Create Calculator
    calc = Calculator(descriptors)

    # map method calculate multiple molecules (return generator)
    # print(list(calc.map(mols)))

    # pandas method calculate multiple molecules (return pandas DataFrame)
    # print(calc.pandas(mols))
        c = calc.pandas(active)
        c.to_csv("active.csv")
