import pandas as pd
import pickle
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict
from sklearn.metrics import accuracy_score, confusion_matrix

# Target cancer types to include
nm = ["BRCA", "CESC", "ESCA", "HNSC", "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "THCA"]

# Load MAD-selected gene list (MAD_5000)
df1 = pd.read_csv('H:/FYP/summer/Data/TOD_CUP/MAD_5000', sep='\t', index_col=0)
sz = df1.index  # Selected genes as index
round = True

# Build metastatic gene expression matrix across selected cancer types
for t in nm:
    df2 = pd.read_csv('H:/FYP/Data_Pre_move/' + t + ".RSEM_genes_normalized__data_Level_3", sep='\t', index_col=0)
    
    # Simplify gene names (remove Ensembl IDs)
    for k in df2.index:
        if k != 'Label':
            name = k.split("|")[0]
            df2.rename(index={k: name}, inplace=True)
    
    # Subset expression table to MAD-selected genes
    if round:
        pd1 = pd.DataFrame(df2, index=sz)
        round = False
    else:
        pd2 = pd.DataFrame(df2, index=sz)
        pd1 = pd.concat([pd1, pd2], axis=1)

# Export the merged matrix (rows = selected genes, columns = samples)
outputpath1 = 'H:/FYP/summer/Data/TOD_CUP/MAD_5000_TSP_Metastatic'
pd1.to_csv(outputpath1, sep='\t', index=True, header=True)


number = 200

# Load gene expression data for testing
tst = pd.read_csv("H:/FYP/summer/Data/Meta_TSP/Gene_matrix_Join", index_col=0, sep='\t').T
tst_x = tst.drop(["Label"], axis=1)
tst_y = tst["Label"]

# Apply log2 normalization
df2 = np.zeros((tst_x.shape[0], tst_x.shape[1]), dtype=np.float32)
for a in range(tst_x.shape[0]):
    for b in range(tst_x.shape[1]):
        val = float(tst_x.iloc[a, b])
        tst_x.iloc[a, b] = 1 if val <= 1 else val
        df2[a, b] = math.log2(float(tst_x.iloc[a, b]))

# Load trained logistic regression model
loaded_model = pickle.load(open("H:/FYP/summer/Model/TSP/Join_Saga.dat", "rb"))


# Accuracy on the full test set
tst_accuracy = accuracy_score(tst_y, loaded_model.predict(df2))
print("Training accuracy: ", tst_accuracy)

# Cross-validated predictions and confusion matrix
tst_train_pred = cross_val_predict(loaded_model, df2, tst_y, cv=10)
tst_conf_matrix = confusion_matrix(tst_y, tst_train_pred)
print("Training confusion matrix:\n", tst_conf_matrix)


# Confusion matrix for specified cancer types only
test_conf_matrix = confusion_matrix(tst_y, loaded_model.predict(df2), labels=nm)

# Per-class false negative counts and rates
fn = [sum(test_conf_matrix[i, :]) - test_conf_matrix[i, i] for i in range(len(nm))]
fnr = [fn[i] / sum(test_conf_matrix[i, :]) for i in range(len(nm))]
print(fnr)

# Overall false negative rate
overall_fnr = sum(fn) / np.sum(test_conf_matrix)
print(overall_fnr)


# Heatmap of confusion matrix
sns.heatmap(test_conf_matrix, annot=False, cmap='RdYlBu')
plt.xlabel('Predicted Labels')
plt.ylabel('True Labels')
plt.title('Confusion Matrix')
plt.show()
