import pandas as pd

# List of TCGA cancer types
nm = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
      "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
      "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

# --- Step 1: Convert all expression values to integers ---
for t in nm:
    # Load preprocessed expression table
    df = pd.read_csv('/data/haochun/fyp/Data_Pre_Primary/' + t + ".RSEM_genes_normalized__data_Level_3", sep='\t', index_col=0)
    df.drop(labels="Label", inplace=True, axis=0)  # Remove label row

    # Round all expression values to the nearest integer
    for k in range(df.shape[0]):
        for s in range(df.shape[1]):
            df.iloc[k, s] = round(float(df.iloc[k, s]))

    # Save the integer-transformed matrix
    outputpath = '/data/haochun/fyp/New/Int/' + t + "_gene"
    df.to_csv(outputpath, sep='\t', index=True, header=True)

# --- Step 2: Generate 1-vs-all "condition" labels for each cancer type ---
for t in nm:
    nm1 = nm.copy()  # Clone the full list
    df1 = pd.read_csv('/data/haochun/fyp/New/Int/' + t + '_gene', sep='\t', index_col=0)
    mt1 = pd.DataFrame(index=df1.columns)
    mt1['condition'] = t  # Mark all samples from current cancer type as target

    # For all other cancer types, mark as "Others"
    for k in nm1:
        if k != t:
            df2 = pd.read_csv('/data/haochun/fyp/New/Int/' + k + '_gene', sep='\t', index_col=0)
            mt2 = pd.DataFrame(index=df2.columns)
            mt2['condition'] = "Others"
            mt1 = pd.concat([mt1, mt2], axis=0)

    # Replace '-' with '.' in sample IDs to match matrix column names
    for k in mt1.index:
        mt1.rename(index={k: k.replace('-', '.')}, inplace=True)

    # Save label file for binary classification
    outputpath = '/data/haochun/fyp/New/Diff_condition' + t + "_Others"
    mt1.to_csv(outputpath, sep='\t', index=True, header=True)

# --- Step 3: Create a full expression matrix for 1-vs-others comparison ---
# Determine common gene set across all cancer types
round = True
for t in nm:
    df = pd.read_csv('/data/haochun/fyp/New/Int/' + t + '_gene', index_col=0, sep='\t')
    if round:
        sz1 = df.index  # Initialize with gene names
        round = False
    else:
        sz2 = df.index
        sz1 = [val for val in sz2 if val in sz1]  # Keep intersection only

# Merge expression matrices: target cancer type + others
for s in nm:
    df1 = pd.read_csv('/data/haochun/fyp/New/Int/' + s + '_gene', sep='\t', index_col=0)
    df1 = pd.DataFrame(df1, index=sz1)  # Subset to common genes

    for k in nm:
        if k != s:
            df2 = pd.read_csv('/data/haochun/fyp/New/Int/' + k + '_gene', sep='\t', index_col=0)
            df2 = pd.DataFrame(df2, index=sz1)  # Subset to common genes
            df1 = pd.concat([df1, df2], axis=1)  # Merge columns

    # Save the final merged matrix for current cancer type vs others
    outputpath = '/data/haochun/fyp/New/Diff_matrix/' + s + "_Others_matrix"
    df1.to_csv(outputpath, sep='\t', index=True, header=True)







