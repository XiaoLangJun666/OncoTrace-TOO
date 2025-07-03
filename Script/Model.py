import pandas as pd
import numpy as np
import math
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, confusion_matrix, roc_auc_score, recall_score,
    precision_score, f1_score, precision_recall_curve, auc
)

# Load gene expression matrix (transposed: samples × genes)
q = "200"
file_path = f'H:/FYP/summer/Data/Chayi/Gene_matrix_{q}'
df = pd.read_csv(file_path, index_col=0, sep='\t').T

# Separate features and labels
X = df.drop(['Label'], axis=1)
y = df['Label']

# Split into training and test sets
X_train_raw, X_test_raw, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Apply log2 normalization (clip values ≤1 to avoid log(0))
def log2_normalize(df):
    df = df.copy()
    df[df <= 1] = 1
    return np.log2(df)

X_train = log2_normalize(X_train_raw)
X_test = log2_normalize(X_test_raw)

# Train one-vs-rest logistic regression
lr_model = LogisticRegression(C=10000, multi_class='ovr', solver='sag', max_iter=2000)
scores = cross_val_score(lr_model, X_train, y_train, cv=10)
lr_model.fit(X_train, y_train)

# Predictions and predicted probabilities
y_pred = lr_model.predict(X_test)
y_pred_proba = lr_model.predict_proba(X_test)

# Cancer types (label classes)
nm = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
      "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
      "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

# Confusion matrix + false negative rate (per class and overall)
conf_matrix = confusion_matrix(y_test, y_pred, labels=nm)
fnr = [(sum(conf_matrix[i, :]) - conf_matrix[i, i]) / sum(conf_matrix[i, :]) for i in range(len(nm))]
overall_fnr = sum([sum(conf_matrix[i, :]) - conf_matrix[i, i] for i in range(len(nm))]) / np.sum(conf_matrix)

# Standard evaluation metrics
train_acc = accuracy_score(y_train, lr_model.predict(X_train))
test_acc = accuracy_score(y_test, y_pred)
auc_score = roc_auc_score(y_test, y_pred_proba, multi_class='ovr')
recall_weighted = recall_score(y_test, y_pred, average='weighted')
f1_weighted = f1_score(y_test, y_pred, average='weighted')
precision_weighted = precision_score(y_test, y_pred, average='weighted')

# PRC AUC (precision-recall curve) for each class
prc_auc = {}
for i, class_name in enumerate(nm):
    precision_i, recall_i, _ = precision_recall_curve(y_test == class_name, y_pred_proba[:, i])
    prc_auc[class_name] = auc(recall_i, precision_i)
avg_prc_auc = np.mean(list(prc_auc.values()))

# Specificity (per class and average)
specificities = {}
for i, class_name in enumerate(nm):
    tn = conf_matrix.sum() - (conf_matrix[i, :].sum() + conf_matrix[:, i].sum() - conf_matrix[i, i])
    fp = conf_matrix[:, i].sum() - conf_matrix[i, i]
    specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    specificities[class_name] = specificity
avg_specificity = np.nanmean(list(specificities.values()))

# Per-class recall and precision
recall_per_class = recall_score(y_test, y_pred, average=None)
precision_per_class = precision_score(y_test, y_pred, average=None)

# Print evaluation results
print(f"Rank: {q}")
print(f"Train Accuracy: {train_acc}")
print(f"Test Accuracy: {test_acc}")
print(f"FNR per class: {fnr}")
print(f"Overall FNR: {overall_fnr}")
print(f"AUC (OvR): {auc_score}")
print(f"Sensitivity (Recall Weighted): {recall_weighted}")
print(f"F1-Score (Weighted): {f1_weighted}")
print(f"Precision (Weighted): {precision_weighted}")
print(f"Average PRC AUC: {avg_prc_auc}")
print(f"Average Specificity: {avg_specificity}")
print("Sensitivity per class:", recall_per_class)
print("Precision per class:", precision_per_class)


# import pickle
# pickle.dump(lr_model, open(f'/home/haochun/OVR/Model/cn/model_{q}.pkl', 'wb'))





# from sklearn.metrics import roc_curve, auc
# import matplotlib.pyplot as plt
# from sklearn.metrics import roc_curve, plot_roc_curve
# import numpy as np
# from sklearn.preprocessing import label_binarize
# #Calculate ROC curve and AUC for each class
# y_test_binarized = label_binarize(y_test, classes=np.unique(y))
# fpr = dict()
# tpr = dict()
# roc_auc = dict()
# n_classes = y_test_binarized.shape[1]
# for i in range(n_classes):
#     fpr[i], tpr[i], _ = roc_curve(y_test_binarized[:, i], y_pred_proba[:, i])
#     roc_auc[i] = auc(fpr[i], tpr[i])
#
# # Compute micro-average ROC curve and AUC
# fpr["micro"], tpr["micro"], _ = roc_curve(y_test_binarized.ravel(), y_pred_proba.ravel())
# roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

#Plot ROC curve for each class and micro-average ROC curve
# import matplotlib.pyplot as plt
# plt.figure(dpi=350,figsize=(20,14))
# plt.plot(fpr["micro"], tpr["micro"],
#          label='micro-average ROC curve (area = {0:0.2f})'
#                ''.format(roc_auc["micro"]),
#          color='deeppink', linestyle=':', linewidth=4)
#
# colors=['yellowgreen','yellow','peru','dodgerblue','gold','hotpink','firebrick','cyan','slateblue','forestgreen','chocolate','lawngreen','sandybrown','lightpink','lightskyblue','red','turquoise','darkcyan','darkgoldenrod','moccasin','darkkhaki','seagreen','mediumslateblue','khaki','wheat','crimson','olive','darkolivegreen','mediumvioletred','aqua','slateblue','papayawhip','sienna']
# #colors = ['blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue', 'red', 'green', 'orange', 'purple','blue','red', 'green']
# for i, color in zip(range(n_classes), colors):
#     plt.plot(fpr[i], tpr[i], color=color, lw=2,
#              label='ROC curve of class {0} (area = {1:0.4f})'
#              ''.format(i, roc_auc[i]))
#
# plt.plot([0, 1], [0, 1], 'k--', lw=2)
# plt.xlim([-0.05, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver Operating Characteristic')
# plt.legend(loc="lower right")
# plt.savefig('H:/FYP/New_Data/200_ROC_DPI_600_4')
# #plt.show()
#
#





#SHAP
# explainer = shap.Explainer(lr_model.predict_proba, X_train)
# shap_values = explainer(X_test,max_evals=1500)
# shap.summary_plot(shap_values[:, :, 1], X_test,max_display=20)
#
#
# explainer = shap.Explainer(lr_model.predict_proba, X_train)
# shap_values = explainer(X_test,max_evals=1500)
# shap.summary_plot(shap_values[:, :, 0], X_test,max_display=20)
#
# explainer = shap.Explainer(lr_model.predict_proba, X_train)
# shap_values = explainer(X_test,max_evals=1500)
# shap.summary_plot(shap_values[:, :, 2], X_test,max_display=20)





# dic1={"ACC":0,"BLCA":0,"BRCA":0,"CESC":0,"CHOL":0,"COAD":0,"DLBC":0,"ESCA":0,"GBM":0,"GBMLGG":0,"HNSC":0,"KICH":0,
#     "KIRC":0,"KIRP":0,"LAML":0,"LGG":0,"LIHC":0,"LUAD":0,"LUSC":0,"MESO":0,"OV":0,"PAAD":0,"PCPG":0,"KIPAN":0,
#     "PRAD":0,"READ":0,"SARC":0,"SKCM":0,"STAD":0,"STES":0,"TGCT":0,"THCA":0,"THYM":0,"UCEC":0,"UCS":0,"UVM":0}
# for k in y:
#     dic1[k]+=1
# print(dic1)
#
# dic2={"ACC":0,"BLCA":0,"BRCA":0,"CESC":0,"CHOL":0,"COAD":0,"DLBC":0,"ESCA":0,"GBM":0,"GBMLGG":0,"HNSC":0,"KICH":0,
#     "KIRC":0,"KIRP":0,"LAML":0,"LGG":0,"LIHC":0,"LUAD":0,"LUSC":0,"MESO":0,"OV":0,"PAAD":0,"PCPG":0,"KIPAN":0,
#     "PRAD":0,"READ":0,"SARC":0,"SKCM":0,"STAD":0,"STES":0,"TGCT":0,"THCA":0,"THYM":0,"UCEC":0,"UCS":0,"UVM":0}
#
# for s in Y_train:
#     dic2[s]+=1
# print(dic2)