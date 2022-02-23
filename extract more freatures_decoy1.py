from multiprocessing import freeze_support
from rdkit import Chem
from rdkit.Chem import Descriptors 
from mordred import Calculator, descriptors
#from mordred 文章注明
import os
os.environ["CUDA_VISIBLE_DEVICES"] = '1'
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
			if c!= -666:
				renew.append(sm)
		except:
			next
	c = calc.pandas(renew).fill_missing()
	print(len(c))
	c.to_csv("decoy_dpt1.csv")
