import pandas as pd
import argparse

##TODO:
"""
1. if alleles are not specified by user - done, just check the output results thoroughly
2. Merging more than 2 files
3. odds ratio instead of beta
4. CI instead of SE
5. Z scores
"""


def col_lookup(col_to_lookup, columns):
    # Check if any of the values corresponding to col_to_lookup key from reverse mapping exist in columns (case-insensitive)
    reverse_mapping = {
        'SNP': ['SNP', 'MARKERNAME', 'SNPID', 'RS', 'RSID', 'RS_NUMBER', 'RS_NUMBERS', 'RSIDS'],
        'CHR': ['CHR', 'CHROM', 'CHROMOSOME', '#CHR', '#CHROM'],
        'POS': ['POS', 'BP', 'POSITION'],
        'BETA': ['BETA', 'EFFECT', 'EFFECT_SIZE', 'EFFECT_ESTIMATE'],
        'SE': ['SE', 'STDERR', 'STANDARD_ERROR', 'SEBETA', 'STD_ERR'],
        'ALLELE': ['ALLELE', 'ALLELES', 'MINOR_ALLELE', 'MAJOR_ALLELE', 'MINOR', 'MAJOR', 'A1', 'A2', 'ALLELE1', 'ALLELE2', 'REF', 'ALT', 'REF_ALLELE', 'ALT_ALLELE', 'EFFECT_ALLELE', 'OTHER_ALLELE', 'EA', 'OA']
    }
    values = reverse_mapping.get(col_to_lookup, [])
    for value in values:
        if value.lower() in [col.lower() for col in columns]:
            if value in columns:
                return value
            elif value.lower() in [col.lower() for col in columns]:
                return value.lower()
    return None



def gwas_merge(gwas1, gwas2, SNP1, SNP2, BETA1, BETA2, SE1, SE2, MINOR1, MINOR2, MAJOR1, MAJOR2):

    print(SNP1, SNP2, BETA1, BETA2, SE1, SE2, MINOR1, MINOR2, MAJOR1, MAJOR2)
    gwas1 = gwas1.rename(columns={SNP1: 'SNP1', MINOR1: 'MINOR1', MAJOR1: 'MAJOR1', 
                                    BETA1: 'BETA1', SE1: 'SE1'})
    gwas2 = gwas2.rename(columns={SNP2: 'SNP2', MINOR2: 'MINOR2', MAJOR2: 'MAJOR2', 
                                    BETA2: 'BETA2', SE2: 'SE2'})
    gwas1 = gwas1.rename(columns={col_lookup('CHR', gwas1.columns): 'CHR1', col_lookup('POS', gwas1.columns): 'POS1'})

    if MINOR1 and MINOR2 and MAJOR1 and MAJOR2:
        print(gwas1.columns)
        print(gwas2.columns)
        merged = pd.merge(gwas1, gwas2, left_on=['SNP1', 'MINOR1', 'MAJOR1'], right_on=['SNP2', 'MINOR2', 'MAJOR2'])
    
    else:
        if not MINOR1:
            MINOR1 = col_lookup('ALLELE', gwas1.columns)
        if not MAJOR1:
            MAJOR1 = col_lookup('ALLELE', [col for col in gwas1.columns if col != MINOR1])
        if not MINOR2:
            MINOR2 = col_lookup('ALLELE', gwas2.columns)
        if not MAJOR2:
            MAJOR2 = col_lookup('ALLELE', [col for col in gwas2.columns if col != MINOR2])

        gwas1 = gwas1.rename(columns={MINOR1: 'MINOR1', MAJOR1: 'MAJOR1'})
        gwas2 = gwas2.rename(columns={MINOR2: 'MINOR2', MAJOR2: 'MAJOR2'})

        print(MINOR1)
        print(MINOR2)
        print(MAJOR1)
        print(MAJOR2)

        merge_flip_neg = gwas1.merge(
            gwas2,
            left_on=['SNP1', 'MINOR1', 'MAJOR1'],
            right_on=['SNP2', 'MINOR2', 'MAJOR2'],
            how='inner'
        )
        merge_flip_neg['BETA2'] = -1 * merge_flip_neg['BETA2']

        merge_flip = gwas1.merge(
            gwas2,
            left_on=['SNP1', 'MINOR1', 'MAJOR1'],
            right_on=['SNP2', 'MINOR2', 'MAJOR2'],
            how='inner'
        )

        merge_neg = gwas1.merge(
            gwas2,
            left_on=['SNP1', 'MINOR1', 'MAJOR1'],
            right_on=['SNP2', 'MAJOR2', 'MINOR2'],
            how='inner'
        )
        merge_neg['BETA2'] = -1*merge_neg['BETA2']

        merge = gwas1.merge(
            gwas2,
            left_on=['SNP1', 'MINOR1', 'MAJOR1'],
            right_on=['SNP2', 'MAJOR2', 'MINOR2'],
            how='inner'
        )
        # Check correlations and merge based on conditions
        correlations = {
            'merge_flip_neg': merge_flip_neg['BETA1'].corr(merge_flip_neg['BETA2']),
            'merge_flip': merge_flip['BETA1'].corr(merge_flip['BETA2']),
            'merge_neg': merge_neg['BETA1'].corr(merge_neg['BETA2']),
            'merge': merge['BETA1'].corr(merge['BETA2'])
        }
        # Filter based on correlation conditions
        selected_merges = []
        for key, corr in correlations.items():
            if corr < 0:
                selected_merges.append(locals()[key])
            if len(selected_merges) == 2:  # Stop after selecting two merges
                break

        if len(selected_merges) < 2:
            # If fewer than 2 merges meet the condition, fill with other merges
            for key in correlations.keys():
                if len(selected_merges) < 2:
                    selected_merges.append(locals()[key])

        # Concatenate the selected merges
        merged = pd.concat(selected_merges, ignore_index=True)

    phemed_out = merged[['SNP1', 'CHR1', 'POS1', 'BETA1', 'BETA2', 'SE1', 'SE2']]
    phemed_out = phemed_out.rename(columns={
            'SNP1': 'SNP',
            'CHR1': 'CHR',
            'POS1': 'POS',
            'BETA1': 'STUDY1',
            'BETA2': 'STUDY2',
            'SE1': 'SE1',
            'SE2': 'SE2'
        })

    return phemed_out


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Merge GWAS summary statistics files.")
    parser.add_argument("--gwas1", required=True, help="Path to the first GWAS summary statistics file.")
    parser.add_argument("--gwas2", required=True, help="Path to the second GWAS summary statistics file.")
    parser.add_argument("--SNP1", help="SNP column name for the first file.")
    parser.add_argument("--SNP2", help="SNP column name for the second file.")
    parser.add_argument("--BETA1", help="BETA column name for the first file.")
    parser.add_argument("--BETA2", help="BETA column name for the second file.")
    parser.add_argument("--SE1", help="SE column name for the first file.")
    parser.add_argument("--SE2", help="SE column name for the second file.")
    parser.add_argument("--MINOR1", help="MINOR allele column name for the first file.")
    parser.add_argument("--MINOR2", help="MINOR allele column name for the second file.")
    parser.add_argument("--MAJOR1", help="MAJOR allele column name for the first file.")
    parser.add_argument("--MAJOR2", help="MAJOR allele column name for the second file.")

    args = parser.parse_args()

    # Load GWAS files
    gwas1 = pd.read_csv(args.gwas1, sep='\t', low_memory=False)
    gwas2 = pd.read_csv(args.gwas2, sep='\t', low_memory=False)
    SNP1 = args.SNP1 or col_lookup('SNP', gwas1.columns)
    SNP2 = args.SNP2 or col_lookup('SNP', gwas2.columns)
    BETA1 = args.BETA1 or col_lookup('BETA', gwas1.columns)
    BETA2 = args.BETA2 or col_lookup('BETA', gwas2.columns)
    SE1 = args.SE1 or col_lookup('SE', gwas1.columns)
    SE2 = args.SE2 or col_lookup('SE', gwas2.columns)
    MINOR1 = args.MINOR1
    MINOR2 = args.MINOR2
    MAJOR1 = args.MAJOR1
    MAJOR2 = args.MAJOR2

    # Merge the two dataframes
    merged = gwas_merge(
        gwas1, gwas2,
        SNP1, SNP2,
        BETA1, BETA2,
        SE1, SE2,
        MINOR1, MINOR2,
        MAJOR1, MAJOR2
    )

    # Save the merged output
    #merged.to_csv("output/merged_gwas.csv", index=False)
    print(merged.head())

if __name__ == "__main__":
    main()
    