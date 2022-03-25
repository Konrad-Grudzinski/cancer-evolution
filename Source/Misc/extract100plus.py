import pandas as pd
import numpy as np

data = pd.read_csv("LUNG_mc3.csv", sep="\t")
data.shape #(379732, 12)
counts_df = data.groupby("sample").agg(["count"])["chr"]
filtered_samples = counts_df.query("count >= 100")
filtered_rows = data.loc[data["sample"].isin(filtered_samples.index)]
filtered_rows.to_csv(r"C:\Users\Konrad Grudzinski\OneDrive - University of Glasgow\Computing\4th Year\Individual Project\DataSets\mc3_gene_level\Lung_mc3_filtered100plus.csv", sep="\t", index=False)








