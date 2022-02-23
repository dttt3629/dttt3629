    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import AllChem
    #利用二元列表录入数据并输出
    active = Chem.SDMolSupplier("decoy_f.sdf")
 #   smile_m = ["smile"]
  #  molwt_m = ["molwt"]
   # tpsa_m = ["tpsa"]
   # logp_m = ["logp"]
   # morgan1_m = ['morgan1']
   # morgan2_m = ['morgan2']
   # morgan3_m = ['morgan3']
   # maxapc_m = ["MaxAPC"]
   # minapc_m = ['MinAPC']
    morgan_m = ['morganfinger']
    
    #载入特征
    for sm in active:
        # smile_m.append(Chem.MolToSmiles(sm))
        # molwt_m.append(Descriptors.MolWt(sm))
        # tpsa_m.append(Descriptors.TPSA(sm))
        # logp_m.append(Descriptors.MolLogP(sm))
        # morgan1_m.append(Descriptors.FpDensityMorgan1(sm))
        # morgan2_m.append(Descriptors.FpDensityMorgan2(sm))
        # morgan3_m.append(Descriptors.FpDensityMorgan3(sm))
        # maxapc_m.append(Descriptors.MaxAbsPartialCharge(sm))
        # minapc_m.append(Descriptors.MinAbsPartialCharge(sm))
        try:
            morgan_m.append(list(AllChem.GetMorganFingerprintAsBitVect(sm,2)))
        except:
            next
        #这里会报错但能用
    # descript = [smile_m,molwt_m,tpsa_m,logp_m,morgan1_m,morgan2_m,morgan3_m,maxapc_m,minapc_m,morgan_m]
    descript = [morgan_m]
    index = len(descript)
    context = len(descript[0])
    i = 0
    descript_file = open("decoy_f_morgan.csv_part2",mode="w")
    while i < index :
        descript_file.write(descript[i][0]+',')
        i+=1
    descript_file.write('\n')
    j = 1
    while j < context :
        i = 0
        while i < index:
            descript_file.write(str(descript[i][j])+',')
            i += 1
        j += 1
        descript_file.write('\n')
    descript_file.close()        