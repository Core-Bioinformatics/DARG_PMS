import sys
import os
import pandas as pd

source_folder=
featureCountsDir=os.path.join(source_folder, "06.Feature_counts")
exprName=os.path.join(source_folder, "07.Expr_matrix/expr_matrix.csv")

if not os.path.exists(os.path.dirname(exprName)):
    os.makedirs(os.path.dirname(exprName))

firstFile=True
df_column_names = []

for file_name in os.listdir(featureCountsDir):
    current_path = os.path.join(featureCountsDir, file_name)
    if not file_name.endswith("txt") or os.path.isdir(current_path):
        continue
    
    print(current_path)
    temp_df = pd.read_csv(current_path, sep="\t", comment="#", header="infer")
    temp_df = temp_df.iloc[:, [0, -1]]
    if firstFile:
        df = temp_df
        firstFile=False
    else:
        df = pd.merge(df, temp_df, on="Geneid")
    df_column_names.append(file_name.split("_raw_counts")[0])

df.columns = ["Geneid"] + df_column_names
df.to_csv(exprName, sep=",", index=False)
