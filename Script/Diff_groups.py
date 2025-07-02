import pandas as pd
nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
    "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

##Transfer data into int 
# for t in nm:
#     df=pd.read_csv('/data/haochun/fyp/Data_Pre_Primary/'+t+".RSEM_genes_normalized__data_Level_3",sep='\t',index_col=[0])
#     df.drop(labels="Label",inplace=True,axis=0)
#     for k in range(0,df.shape[0]):
#         for s in range(0,df.shape[1]):
#             df.iloc[k,s]=round(float(df.iloc[k,s]))
#             s+=1
#         k+=1
#     outputpath = '/data/haochun/fyp/New/Int/'+t+"_gene"
#     df.to_csv(outputpath, sep='\t', index=True, header=True)



#Generate 1-vs-others expression table for all cancer types

# for t in nm:
#     nm1=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
#     "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
#     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
#     df1 = pd.read_csv('/data/haochun/fyp/New/Int/' + t + '_gene', sep='\t', index_col=[0])
#     mt1 = pd.DataFrame(index=df1.columns)
#     mt1['condition'] = t
#     for k in nm1:
#         if k!=t:
#             df2 = pd.read_csv('/data/haochun/fyp/New/Int/' + k + '_gene', sep='\t', index_col=[0])
#             mt2 = pd.DataFrame(index=df2.columns)
#             mt2['condition'] = "Others"
#             mt1=pd.concat([mt1,mt2],axis=0)
#     for k in mt1.index:
#         re = k.replace('-', '.')
#         mt1.rename(index={k: re}, inplace=True)
#
#     outputpath='/data/haochun/fyp/New/Diff_condition'+t+"_Others"
#     mt1.to_csv(outputpath,sep='\t',index=True,header=True)




# round=True
# for t in nm:
#     df = pd.read_csv('/data/haochun/fyp/New/Int/' + t + '_gene', index_col=[0], sep='\t')
#     if round==True:
#         sz1=df.index
#         round=False
#     else:
#         sz2=df.index
#         sz1=[val for val in sz2 if val in sz1]
#
#
# for s in nm:
#     df1 = pd.read_csv('/data/haochun/fyp/New/Int/'+s+'_gene', sep='\t', index_col=[0])
#     df1=pd.DataFrame(df1,index=sz1)
#     nm2=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
#          "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
#         "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
#     for k in nm2:
#         if k!=s:
#             df2 = pd.read_csv('/data/haochun/fyp/New/Int/' + k + '_gene', sep='\t', index_col=[0])
#             df2=pd.DataFrame(df2,index=sz1)
#             df1 = pd.concat([df1, df2], axis=1)

#     outputpath = '/data/haochun/fyp/New/Diff_matrix/' + s + "_Others_matrix"
#     df1.to_csv(outputpath,sep='\t',index=True,header=True)







