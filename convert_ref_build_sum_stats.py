import pandas as pd
from pyliftover import LiftOver

# Load data
df = pd.read_csv("/Users/sonali.gupta/Downloads/SCZ_sum_stats/daner_PGC_SCZ52_0513a.hq2", sep = "\t")

# Create LiftOver object
lo = LiftOver("hg19ToHg38.over.chain.gz")

# Perform liftover
with open("/Users/sonali.gupta/Downloads/SCZ_sum_stats/lifted_daner_PGC_SCZ52_0513a.hq2", "w") as outfile:
    # Write header
    outfile.write("\t".join(df.columns) + "\tnew_CHROM\tnew_BP\n")
    
    for index, row in df.iterrows():
        result = lo.convert_coordinate("chr" + str(row.iloc[0]), row.iloc[2])
        if result:
            new_chrom, new_bp = result[0][0], result[0][1]
            outfile.write("\t".join(map(str, row)) + f"\t{new_chrom}\t{new_bp}\n")