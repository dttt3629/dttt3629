# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 16:34:05 2022

@author: 10638
"""
import pandas as pd
import numpy as np
from numpy import random
import os
import pickle
from multiprocessing import freeze_support
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from rdkit.Chem import Descriptors
from sklearn.preprocessing import StandardScaler, MaxAbsScaler
from sklearn.feature_selection import VarianceThreshold, SelectFromModel
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, ExtraTreesRegressor
from sklearn import svm
from sklearn.model_selection import cross_val_score, StratifiedKFold, KFold
from hyperopt import fmin, hp, tpe,Trials
from sklearn.metrics import roc_auc_score, mean_absolute_error, mean_squared_error, mean_squared_log_error,plot_roc_curve
# from xgboost import XGBClassifier
from sklearn.datasets import load_wine
from multiprocessing import Pool
import math
from sklearn.model_selection import train_test_split
from boruta import BorutaPy


def extract_features(datas,is_save = False,save_name='unnamed'):
    #用于计算分子描述符，data是SDF文件
    # datas = Chem.SDMolSupplier(path)
    # Create Calculator
    calc = Calculator(descriptors)
    # map method calculate multiple molecules (return generator)
    # print(list(calc.map(mols)))
    # pandas method calculate multiple molecules (return pandas DataFrame)
    # print(calc.pandas(mols))
    renew = []
    wrong = []
    index =0
    for sm in datas:
        index += 1
        try:
            c = Descriptors.fr_sulfone(sm)
            if c != -666:
                renew.append(sm)
    #删掉错误分子,并保留它的index
        except:
            wrong.append(index-1)
    c = calc.pandas(renew).fill_missing()
    if is_save:
        c.to_csv(f'{save_name}_dpt.csv',index=False)
    return c,wrong

def extract_cluster(file,wrong,is_save = False,save_name='unnamed'):
    #用于从文件中找到聚类，并返回一个列表
    output = []
    string = file.readline()
    while string!='':
        string = file.readline()
        if string == '> <Cluster>\n':
            string = file.readline()
            output.append(string)
    output = [x.rstrip() for x in output]
    #删去\n
    for index in wrong:
        del output[index]
    if is_save:
        with open(f"{save_name}_cluster.pkl",'wb') as f:
            pickle.dump(output,f)
    return output

def del_null(ls,per):
    #设定阈值，如果空值大于per，删去
    for col in ls:
        flag =0
        for sm in ls[col]:
            if pd.isnull(sm):
                flag+=1
        if flag >= int(len(ls)*per):
            del ls[col]
def del_zero(ls,per):
    #设定阈值，如果零值大于per，删去
    for col in ls:
        flag =0
        for sm in ls[col]:
            if sm == 0:
                flag+=1
        if flag >= int(len(ls)*per):
            del ls[col]
def change_null(c):
    #将空值换成均值
    col_num=0
    for col in c:
        col_num+=1
        k=[]
        index = 0
        for sm in c[col]:
            index+=1
            if pd.isnull(sm):
                k.append(index)
        if k!= []:
            for num in k:
                c.iloc[num-1,col_num-1]=c[col].mean()

def standard_data(ls,null_size=4,is_save = False,save_name='unnamed'):
    #传入一个csv文件，null_size限定一下最大的null数量，
    #这里建议手操，有点bug，注意文件中可能包含True和false之类的字符串要手动置换
    for col in ls:
        flag =0
        for sm in ls[col]:
            if pd.isnull(sm):
                flag+=1
        if flag > null_size:
            del ls[col]
    c = 0
    #零值插平均数，这里iloc表示self，用ls[][]是ls的副本会报错，底层不同
    for col in ls:
        c+=1
        r=0
        for sm in ls[col]:
            r+=1
            if pd.isnull(sm):
                ls.iloc[r-1,c-1]=ls[col].mean()
        ls.iloc[:,c-1]=(ls[col]-ls[col].min())/(ls[col].max()-ls[col].min())
    #标准化
    # for col in ls:
    #     ls[col]=(ls[col]-ls[col].min())/(ls[col].max()-ls[col].min())
    if is_save:
        ls.to_csv(f"{save_name}_standard_data.csv",index = False)
    return ls

def boruta_select(pd_x,pd_y,ls,is_save = False,save_name='unnamed'):
    #boruta算法筛选特征
    pd_xa = np.array(pd_x)
    pd_ya = np.array(pd_y)
    rfc = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=5,random_state=(42))
    feat_selector = BorutaPy(rfc,n_estimators='auto',random_state=1, max_iter=10)
    feat_selector.fit(pd_xa, pd_ya)
    #获取筛选到的特征
    dic_ft_select = pd.DataFrame({'name':pd_x.columns, 'select':feat_selector.support_})
    for col in dic_ft_select.T:
        if dic_ft_select.T[col][1] == False:
            del ls[dic_ft_select.T[col][0]]
    if is_save:
        ls.to_csv(f"{save_name}_boruta_data.csv",index = False)
    return ls

def split_dic(dic,train_size = 0.75):
    random.seed(42)
    train=[]
    test=[]
    for key in dic:
        d = random.choice(dic[key],size=int(train_size*len(dic[key]))+1,replace=False)
        train.extend(d)
        for num in dic[key]:
            if num in d:
                1
            else:
               test.append(num)
    return train,test

global data_n
data_n = 0
def spilit_data(data,is_random=True,is_save=True,save_name='unnamed',train_size = 0.75,cluster=None):
    if is_random:
        train_x, test_x, train_y, test_y = train_test_split(data.iloc[:,:-1], data.iloc[:,-1], train_size=train_size, 
                            stratify=data.iloc[:,-1],shuffle=True,random_state=42)
    else:
        dic = {}
        k = 0
        while k<len(cluster):
            if cluster[k] != 0:
                if cluster[k] in dic:
                    dic[cluster[k]].append(k)
                else:
                    dic[cluster[k]]=[k]
            k+=1
        train,test=split_dic(dic,train_size=train_size)
        train_x = data.iloc[train,:-1].sample(frac=1)
        train_y = data.iloc[train,-1].sample(frac=1)
        test_x =data.iloc[test,:-1].sample(frac=1)
        test_y =data.iloc[test,-1].sample(frac=1)
    if is_save:
        global data_n
        data_n+=1
        train_d=pd.concat([train_x,train_y],axis=1)
        train_d.to_csv(f"{save_name}_train{data_n}.csv",index=False)
        test_d=pd.concat([test_x,test_y],axis=1)
        test_d.to_csv(f"{save_name}_test{data_n}.csv",index = False)
    train_x=train_x.reset_index(drop=True)    
    train_y=train_y.reset_index(drop=True)
    test_x = test_x.reset_index(drop=True)
    test_y = test_y.reset_index(drop=True)
    return train_x, test_x, train_y, test_y

global calcul_n
calcul_n=0
def calcul(t,p,save_name="unmaned"):
    global calcul_n
    calcul_n+=1
    i = 0
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    while i<len(t):
        if t[i] == 1 and p[i]==1:
            tp+=1
        if t[i] == 0 and p[i]==1:
            fp+=1
        if t[i] == 0 and p[i]==0:
            tn+=1
        if t[i] == 1 and p[i]==0:
            fn+=1
        i+=1
    accuracy = (tp+tn)/len(t)
    recall = tp / (tp+fn)
    precision = tp / (tp+fp)
    try:
        F1 = 2*precision*recall /(precision+recall)
    except:
        F1='False'
    ACC = (tp+tn)/(tp+tn+fp+fn)
    MCC = (tp*tn-fp*fn)/pow((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn),0.5)
    dic = {'accuracy':accuracy,'recall':recall,'precision':precision,'F1':F1,'ACC':ACC,'MCC':MCC}
    with open(f'{save_name}_output{calcul_n}.pkl','wb') as f:
        pickle.dump(dic,f)
    return dic

global best_n
best_n =1
def svm_best(svm_train_xa,svm_train_ya,is_save=True):    
    global best_n
    hyper_parameter = {
            'C': hp.uniform('C', 0.01, 50),
            'kernel': hp.choice('kernel', ['linear', 'rbf']),
            'gamma': hp.uniform('gamma', 0.001, 1)
     }
    trials=Trials()
    def svm_model(hyper_parameter):  # 待寻优函数
        clf = svm.SVC(**hyper_parameter, class_weight='balanced', random_state=42)
        e = cross_val_score(clf, svm_train_xa, svm_train_ya, cv=StratifiedKFold(n_splits=10, shuffle=True, 
                                            random_state=42),
                                    scoring='f1', n_jobs=10).mean()
        return -e
    best = fmin(svm_model, hyper_parameter, algo=tpe.suggest, max_evals=100,
                rstate=None,trials=trials)
    if is_save:
        with open(f'svm_best{best_n}.pkl','wb') as f:
            pickle.dump(best,f,pickle.HIGHEST_PROTOCOL)
        with open(f'svm_trials{best_n}.pkl','wb') as f:
            pickle.dump(trials,f,pickle.HIGHEST_PROTOCOL)
        best_n+=1
    return best,trials

def rf_best(rf_train_xa,rf_train_ya,is_save=True):    
    global best_n
    feature_list = ['sqrt', 'log2', None]  # 特征列表
    weight_list = ['balanced', 'balanced_subsample']  # 权重列表
    hyper_parameter = {'n_estimators': hp.uniformint('n_estimators', 80, 150),  # 决策树的数量
                       'max_depth': hp.uniformint('max_depth', 3, 50),  # 树的深度 过多会过拟合
                       'max_features': hp.choice('max_features', feature_list),  # 单个决策树使用的特征
                       'class_weight': hp.choice('class_weight', weight_list)  # 正负集权重
                       }  # 选择要优化的超参数
    trials=Trials()
    def rf_model(hyper_parameter):  # 待寻优函数
        clf = RandomForestClassifier(**hyper_parameter, oob_score=False, n_jobs=5, random_state=42)
        e = cross_val_score(clf, rf_train_xa,rf_train_ya,cv=StratifiedKFold(n_splits=10, shuffle=True, random_state=42),
                                    scoring='f1', n_jobs=5).mean()
        return -e
    best = fmin(rf_model, hyper_parameter, algo=tpe.suggest, max_evals=100,trials=trials) 
    if is_save:
        with open(f'rf_best{best_n}.pkl','wb') as f:
            pickle.dump(best,f,pickle.HIGHEST_PROTOCOL)
        with open(f'rf_trials{best_n}.pkl','wb') as f:
            pickle.dump(trials,f,pickle.HIGHEST_PROTOCOL)
        best_n+=1
    return best,trials

global pred_n
pred_n=0
def set_svm(train_x, test_x, train_y, test_y,F=0,L=-1):
    global pred_n
    pred_n+=1
    #F、L用于删去索引
    ker = ['linear', 'rbf']
    train_xa = np.array(train_x.iloc[:,F:L])
    train_ya = np.array(train_y)
    best,trials=svm_best(train_xa, train_ya)
    clf = svm.SVC(C=best['C'], gamma=best['gamma'], kernel=ker[best['kernel']], class_weight='balanced',
              random_state=42,
              probability=True)
    clf.fit(train_x.iloc[:,F:L],train_y)
    a=clf.predict(test_x.iloc[:,F:L])
    test=pd.concat([test_x,test_y],axis=1)
    test_p=pd.concat([test,pd.DataFrame(a,columns=[f'Pred{pred_n}'])],axis=1)
    calcul(np.array(test_y),a,save_name='svm')
    b=clf.predict(train_x.iloc[:,F:L])
    train = pd.concat([train_x,train_y],axis=1)
    train_p=pd.concat([train,pd.DataFrame(b,columns=[f'Pred{pred_n}'])],axis=1)
    test_p.to_csv(f"svm_test_pred{pred_n}.csv",index=False)
    train_p.to_csv(f"svm_train_pred{pred_n}.csv",index=False)
    return train_p,test_p
 
def set_rf(train_x, test_x, train_y, test_y,F=0,L=-1):
    global pred_n
    pred_n+=1
    #F、L用于删去索引
    train_xa = np.array(train_x.iloc[:,F:L])
    train_ya = np.array(train_y)
    best,trials=rf_best(train_xa, train_ya)
    feature_list = ['sqrt', 'log2', None]  # 特征列表
    weight_list = ['balanced', 'balanced_subsample'] 
    clf = RandomForestClassifier(n_estimators=int(best['n_estimators']), max_depth=int(best['max_depth']),
                                     max_features=feature_list[best['max_features']],
                                     class_weight=weight_list[best['class_weight']], oob_score=False, n_jobs=5,
                                     random_state=42) 
    clf.fit(train_x,train_y)
    a=clf.predict(test_x.iloc[:,F:L])
    test=pd.concat([test_x,test_y],axis=1)
    test_p=pd.concat([test,pd.DataFrame(a,columns=[f'Pred{pred_n}'])],axis=1)
    calcul(np.array(test_y),a,save_name='rf')
    b=clf.predict(train_x.iloc[:,F:L])
    train = pd.concat([train_x,train_y],axis=1)
    train_p=pd.concat([train,pd.DataFrame(b,columns=[f'Pred{pred_n}'])],axis=1)
    test_p.to_csv("rf_test_pred{pred_n}.csv",index=False)
    train_p.to_csv("rf_train_pred{pred_n}.csv",index=False)
    return train_p,test_p
    
def train_svm(data,data_iter,is_random=True,save_name='svm'):
     it = iter(data_iter)
     F=1
     L=-1
     train_size = next(it)
     train_p,test_p=set_svm(*spilit_data(data,
            is_random=is_random,train_size=train_size,save_name=save_name,cluster=data['cluster']),
            F=F,L=L)
     while 1:
         try:
             train_size = next(it)
             L-=1
             train_x, test_x, train_y, test_y=spilit_data(test_p,
                is_random=is_random,train_size=train_size,save_name=save_name,cluster=data['cluster'])
             train_x =pd.concat([train_p.iloc[:,:-1],train_x])
             train_y =pd.concat([train_p.iloc[:,-1],train_y])
             train_p,test_p=set_svm(train_x, test_x, train_y, test_y,F=F,L=L)
         except:
             break

def train_rf(data,data_iter,is_random=True,save_name='rf'):
     it = iter(data_iter)
     F=1
     L=-1
     train_size = next(it)
     train_p,test_p=set_rf(*spilit_data(data,
            is_random=is_random,train_size=train_size,save_name=save_name,cluster=data['cluster']),
            F=F,L=L)
     while 1:
         try:
             train_size = next(it)
             L-=1
             train_x, test_x, train_y, test_y=spilit_data(test_p,
                is_random=is_random,train_size=train_size,save_name=save_name,cluster=data['cluster'])
             train_x =pd.concat([train_p.iloc[:,:-1],train_x])
             train_y =pd.concat([train_p.iloc[:,-1],train_y])
             train_p,test_p=set_rf(train_x, test_x, train_y, test_y,F=F,L=L)
         except:
             break
    
 
    
# if __name__ == "__main__":
#     #ipython做到boruta
#       data = pd.read_csv('data_boruta_f.csv')
#       F=1
#       L=-2
#       is_random=True
#       train_iter=[0.1,0.3]
#       for train_size 
#       train_x, test_x, train_y, test_y= spilit_data(data_x=data.iloc[:,:-1],data_y=data.iloc[:,-1],is_random=is_random,
#                                                     train_size=0.75,
#                                             save_name='rf',
#                                             cluster=data['cluster'])
     