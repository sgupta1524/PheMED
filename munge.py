import pandas as pd
import numpy as np

# Read df1
df1 = pd.read_csv("/Users/sonaligupta/Downloads/ppd2023_EUR_sumstats.tsv", sep="\t")

# Read df2 with corrected separator and clean column names
df2 = pd.read_csv("/Users/sonaligupta/Downloads/pgc.mdd.2012-04/pgc.mdd.full.2012-04.txt", sep="\s+")
df2.columns = df2.columns.str.strip()

# Print column names for verification
print("df1 columns:", df1.columns)
print("df2 columns:", df2.columns)

# Merge the DataFrames
merged_df = pd.merge(df1, df2, left_on=["SNP"], right_on=["snpid"], how="inner")
print(merged_df.head())

merged_df = merged_df[
    (merged_df["A1"] == merged_df["a1"]) & (merged_df["A2"] == merged_df["a2"])
]
merged_df = merged_df[["SNP", "CHR", "BP", "BETA", "or", "SE", "se"]]
merged_df["or"] = np.log(merged_df["or"])
print(merged_df.head())

column_renames = {
    "SNP": "SNP",
    "CHR": "CHR",
    "BP": "POS",
    "BETA": "BETA1",
    "or": "BETA2",
    "SE": "SE1",
    "se": "SE2"
}

# Rename columns in merged_df
merged_df = merged_df.rename(columns=column_renames)


merged_df.to_csv("/Users/sonaligupta/PheMED/data/sample_pgc_sumstats_munged.csv", index=False, header=True)