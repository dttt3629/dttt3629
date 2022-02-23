# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 19:40:57 2022

@author: 10638
"""
#replace函数表现奇怪慎用
import pandas as pd
c=pd.read_csv('active_dpt.csv')
c=c.replace(to_replace=r'FALSE',value=1,regex=True)
c=c.replace(to_replace=r'TRUE',value=0,regex=True)
c=c.replace(to_replace=r'[a-z,A-Z].$',value=0, regex=True)
c.to_csv("active_dpt_modfied")