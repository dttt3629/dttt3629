from multiprocessing import freeze_support
import os 
from rdkit import Chem

from mordred import Calculator, descriptors
#from mordred 文章注明
os.environ["CUDA_VISIBLE_DEVICES"] = '1'
if __name__ == "__main__":
    freeze_support()

    active = Chem.SDMolSupplier("active_f.sdf")
    # Create Calculator
    calc = Calculator(descriptors)
    # result = calc(active)
    # map method calculate multiple molecules (return generator)
    # print(list(calc.map(mols)))

    # pandas method calculate multiple molecules (return pandas DataFrame)
    # print(calc.pandas(mols))
    ccc = calc.pandas(active).fill_missing()
    #fill_missing将报错的描述填0，全报错则删除
    ccc.to_csv("active_dpt1.csv")
