import pandas as pd
from pyliftover import LiftOver
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Convert reference build of summary statistics using LiftOver.")
    parser.add_argument("--input", type=str, required=True, help="Path to the input summary statistics file.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output file with converted coordinates.")
    parser.add_argument("--chain", type=str, required=True, help="Path to the chain file for the conversion.")
    parser.add_argument("--chrom", type=int, default=0, help="Chromosome column number in the input file.")
    parser.add_argument("--pos", type=int, default=2, help="Base pair column position number in the input file.")
    parser.add_argument("--chromosome-col-has-chr", type=bool, default=False, help="Indicate if the chromosome column has chr prefix.")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input, sep="\t")
    chain = args.chain
    chrom = args.chrom
    pos = args.pos

    # Create LiftOver object
    lo = LiftOver(chain)

    # Perform liftover
    with open(args.output, "w") as outfile:
        # Write header
        outfile.write("\t".join(df.columns))
        outfile.write("\n")
        
        for index, row in df.iterrows():
            if args.chromosome_col_has_chr:
                result = lo.convert_coordinate(str(row.iloc[chrom]), row.iloc[pos])
            else:
                result = lo.convert_coordinate("chr" + str(row.iloc[chrom]), row.iloc[pos])
            if result:
                new_bp = result[0][1]
                row.iloc[pos] = new_bp  # Update position
            outfile.write("\t".join(map(str, row)) + "\n")

if __name__ == "__main__":
    main()