import pandas as pd

# Define TCGA cancer types
nm = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
      "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
      "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

dic = {}  # Dictionary to store top genes per cancer type
new_data = pd.DataFrame(columns=["gene", "length"])  # DataFrame for storing unique gene list

# Loop through each cancer type
for t in nm:
    # Read differential expression result file
    df = pd.read_csv('H:/FYP/summer/Data/Chayi_reduced/Reduced_result/Result_' + t + "_Others", sep=',', index_col=0)
    dic1 = {}

    # Collect absolute log2FoldChange for each gene
    for k in df.index:
        dic1[k] = abs(df.loc[k, "log2FoldChange"])

    # Sort by descending importance and select top 200
    dic1 = sorted(dic1.items(), key=lambda item: item[1], reverse=True)
    dic_final = dic1[0:200]
    lst = [s[0] for s in dic_final]  # Get gene names only
    dic[t] = lst  # Store gene list

# Select genes uniquely upregulated in each cancer type (1-vs-rest)
for n in nm:
    gene1 = dic[n]  # Top genes of current cancer
    gene_rest = []
    gene_select = []

    for m in nm:
        if n != m:
            for j in dic[m]:
                if j not in gene_rest:
                    gene_rest.append(j)

    # Only keep genes not in other cancers' top list
    for i in gene1:
        if i not in gene_rest:
            gene_select.append(i)

    # Store unique gene list and count
    new_data.loc[n, "gene"] = str(gene_select)
    new_data.loc[n, "length"] = len(gene_select)

# Save the final selected gene lists
outputpath = "H:/FYP/summer/Data/Chayi_reduced/Reduced_Gene_list/Gene_list_200"
new_data.to_csv(outputpath, sep='\t', index=True, header=True)
