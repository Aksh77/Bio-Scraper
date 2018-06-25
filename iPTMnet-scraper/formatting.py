import csv
import pandas as pd
import numpy as np

input_file = "Extracted Data/PTMdata.csv"
df = pd.read_csv(input_file)

residues = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
df=pd.concat([df,pd.DataFrame(columns=residues)], sort=False)

for i, row in df.iterrows():
    sites = row[2].split(" ")
    for s in sites:
        residue = s[0]
        if(pd.isnull(df.at[i,residue])):
            df.at[i,residue] = s
        else:
            df.at[i,residue] = row[residue] + " " + s
del df['PTM Site']

df.to_csv("Extracted Data/PTMdata_formatted.csv", sep=',', index=False)
