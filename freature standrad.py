# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 21:17:54 2022

@author: d1063
"""

import pandas as pd
import numpy as np
from sklearn import preprocessing
data = pd.read_csv("D:\\work\\actives_final.mol2\\active_fdpt.csv",sep=',',header = None ,skiprows = 1)

da = (data-data.min())/(data.max()-data.min())

# new_data = pd.DataFrame(np.ones((len(data),1))*0)
# while x < data.shape[1]-1:
#     da = data.iloc[:,x]
#     da_train = da.values.reshape(-1,1)
#     MinMax = preprocessing.MinMaxScaler()
#     da_MinMax = MinMax.fit_transform(da_train)
#     da_MinMax = pd.DataFrame(da_MinMax)
#     new_data=pd.concat(new_data,da_MinMax)
#     x+=1
# print(new_data)