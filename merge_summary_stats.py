import pandas as pd
import argparse
import os
import re
import numpy as np
import gzip

##TODO:
"""
Verified merging for 2 GWAS at a time with phemed preprint
OuD PGC with MVP not verified, PGC with BioVU verified, PGC with finngen verified
When all 4 were merged together, it did not work
Run David's notebook and verify
8. Remove print statements
"""

def col_lookup(col_to_lookup, columns):
    # Check if any of the values corresponding to col_to_lookup key from reverse mapping exist in columns (case-insensitive)
    mapping = {
        'SNP': ['SNP', 'MARKERNAME', 'SNPID', 'RS', 'RSID', 'RS_NUMBER', 'RS_NUMBERS', 'RSIDS', 'SNP_ID', 'VARIANT_ID'],
        'CHR': ['CHR', 'CHROM', 'CHROMOSOME', '#CHR', '#CHROM'],
        'POS': ['POS', 'BP', 'POSITION', 'BASE_PAIR_LOCATION', 'POSITION_BP', 'BP_POSITION'],
        'BETA': ['BETA', 'EFFECT', 'EFFECT_SIZE', 'EFFECT_ESTIMATE', 'LOGOR', 'LOG_ODDS', 'LOG_ODDS_RATIO', 'BETA_COEFFICIENT'],
        'OR' : ['OR', 'ODDS_RATIO', 'ODDS', 'ODDSRATIO'],
        'SE': ['SE', 'STDERR', 'STANDARD_ERROR', 'SEBETA', 'STD_ERR', 'STDERRLOGOR', 'SETest statistic'],
        'ALLELE': ['ALLELE', 'ALLELE0', 'MINOR_ALLELE', 'MAJOR_ALLELE', 'MINOR', 'MAJOR', 'A1', 'A2', 'ALLELE1', 'ALLELE2', 
                   'REF', 'ALT', 'REF_ALLELE', 'ALT_ALLELE', 'EFFECT_ALLELE', 'OTHER_ALLELE', 'EA', 'OA', 'ALLELE 1', 'ALLELE 2']
    }
    values = mapping.get(col_to_lookup, [])
    for value in values:
        for col in columns:
            if value.lower() == col.lower():
                return col
    return None
            

def detect_separator(file_path):
    """
    Detects the separator used in a file by checking the first line.
    Supports both plain text and gzipped files.
    Returns the detected separator or raises an error if no separator is found.
    """

    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt', encoding='utf-8') as f:
        first_line = f.readline().strip()
    for sep in [',', '\t', ';', ' ']:
        if sep in first_line:
            return sep
    raise ValueError(f"Unable to detect separator in {file_path}")

def merge_of_merge(merged_list):
    """
    Outer join of the elements of merged_list, which are dataframes.
    If the columns are not present in the dataframe, they will be filled with NaN.
    """
    merged = pd.merge(merged_list[0], merged_list[1], on=['SNP', 'CHR', 'POS'], how='outer', suffixes=('_dup', '')).drop_duplicates(subset=['SNP'], keep='first')
    for i in range(2, len(merged_list)):
        merged = pd.merge(merged, merged_list[i], on=['SNP', 'CHR', 'POS'], how='outer', suffixes=('_dup', '')).drop_duplicates(subset=['SNP'], keep='first')

    return fix_colnames(merged)

def fix_colnames(df):
    def extract_suffix(col_name):
        if col_name.startswith('STUDY') or col_name.startswith('SE'):
            # Remove STUDY/SE prefix and _dup suffix if present
            suffix = re.sub(r'^(STUDY|SE)', '', col_name)
            suffix = re.sub(r'_dup$', '', suffix)
            return suffix
        return col_name

    # Get all column names
    cols = df.columns.tolist()

    # Separate STUDY and SE columns
    study_cols = [col for col in cols if col.startswith('STUDY')]
    se_cols = [col for col in cols if col.startswith('SE')]

    # Group columns by suffix
    study_groups = {}
    se_groups = {}

    for col in study_cols:
        suffix = extract_suffix(col)
        if suffix not in study_groups:
            study_groups[suffix] = []
        study_groups[suffix].append(col)

    for col in se_cols:
        suffix = extract_suffix(col)
        if suffix not in se_groups:
            se_groups[suffix] = []
        se_groups[suffix].append(col)

    # Combine columns within each group and use simple naming
    # Create new dataframe with combined columns
    result_df = df[['SNP', 'CHR', 'POS']].copy()

    # Combine STUDY columns with numbered names
    study_counter = 1
    for suffix, cols in study_groups.items():
        if len(cols) > 1:   
            # Take first non-null value across duplicate columns
            result_df[f'STUDY{study_counter}'] = df[cols].bfill(axis=1).iloc[:, 0]
        else:
            result_df[f'STUDY{study_counter}'] = df[cols[0]]
        study_counter += 1

    # Combine SE columns with numbered names
    se_counter = 1
    for suffix, cols in se_groups.items():
        if len(cols) > 1:
            # Take first non-null value across duplicate columns
            result_df[f'SE{se_counter}'] = df[cols].bfill(axis=1).iloc[:, 0]
        else:
            result_df[f'SE{se_counter}'] = df[cols[0]]
        se_counter += 1

    # Reorder columns to group STUDY and SE together
    study_cols = [col for col in result_df.columns if col.startswith('STUDY')]
    se_cols = [col for col in result_df.columns if col.startswith('SE')]
    ordered_cols = ['SNP', 'CHR', 'POS'] + study_cols + se_cols
    result_df = result_df[ordered_cols]
    return result_df

def gwas_merge(gwas1, gwas2, SNP1, SNP2, BETA1, BETA2, SE1, SE2, MINOR1, MINOR2, MAJOR1, MAJOR2, 
               Z1=None, Z2=None, MAF1=None, MAF2=None, Ncases1=None, Ncases2=None, Ncontrols1=None, Ncontrols2=None,
               OR1=None, OR2=None):

    gwas1_filename = os.path.basename(gwas1)
    gwas2_filename = os.path.basename(gwas2)
    separator1 = detect_separator(gwas1)
    separator2 = detect_separator(gwas2)
    open_func1 = gzip.open if gwas1.endswith('.gz') else open
    open_func2 = gzip.open if gwas2.endswith('.gz') else open

    with open_func1(gwas1, 'rt', encoding='utf-8') as f1:
        gwas1 = pd.read_csv(f1, sep=separator1, low_memory=False)

    with open_func2(gwas2, 'rt', encoding='utf-8') as f2:
        gwas2 = pd.read_csv(f2, sep=separator2, low_memory=False)
    print(f"Loaded {gwas1_filename} with {gwas1.shape[0]} rows and {gwas1.shape[1]} columns.")
    print(f"Loaded {gwas2_filename} with {gwas2.shape[0]} rows and {gwas2.shape[1]} columns.")

    #if SNP, BETA and SE are none do a lookup
    SNP1 = col_lookup('SNP', gwas1.columns) if SNP1 is None or SNP1.strip() == "" else SNP1
    SNP2 = col_lookup('SNP', gwas2.columns) if SNP2 is None or SNP2.strip() == "" else SNP2
    SE1 = col_lookup('SE', gwas1.columns) if SE1 is None else SE1

    if SE1 is None and MAF1 and Ncases1 and Ncontrols1:
        SE1 = 'SE1'
        gwas1['Prop'] = gwas1[Ncases1] /(gwas1[Ncontrols1]+gwas1[Ncases1])
        gwas1['MAF_term'] = gwas1['MAF1'] * (1 - gwas1['MAF1'])
        gwas1[SE1] = 1 / np.sqrt(2 * gwas1['MAF_term'] * gwas1[Ncases1]*gwas1['Prop']*(1-gwas1['Prop']))

    elif SE1 is None and MINOR1 and MAJOR1 and Ncases1 and Ncontrols1:
        print(f"Calculating SE with MAF from MAFdata file for {gwas1_filename}.")
        
        gwas1[MINOR1] = gwas1[MINOR1].str.upper()
        gwas1[MAJOR1] = gwas1[MAJOR1].str.upper()

        MAFdata = pd.read_csv("/Users/sonali.gupta/Downloads/PheMEDPower/OuD/New/summary_stats_finngen_R7_F5_OPIOIDS (1) (1)", sep="\t")
        # Merge gwas1 with MAFdata to avoid row-wise operations
        merged_data = gwas1.merge(
        MAFdata[['rsids', 'alt', 'ref', 'af_alt']],
        left_on=[SNP1, MINOR1, MAJOR1],
        right_on=['rsids', 'alt', 'ref'],
        how='inner'
        )
        # Calculate N and SE
        merged_data['Prop'] = merged_data[Ncases1] /(merged_data[Ncontrols1]+merged_data[Ncases1])
        merged_data['MAF_term'] = merged_data['af_alt'] * (1 - merged_data['af_alt'])
        SE1 = 'SE1'
        merged_data[SE1] = 1 / np.sqrt(2 * merged_data['MAF_term'] * merged_data[Ncases1]*merged_data['Prop']*(1-merged_data['Prop']))
        gwas1 = merged_data  # Update gwas1 with the merged data

    elif SE1 is None:
        #Warning for missing SE1
        print(f"Warning: SE is not specified for {gwas1_filename}. Please provide # cases/ controls and effect/ non=-ffect allele to calculate SE. If you do not have these, please provide SE column in the GWAS file.")
    
    SE2 = col_lookup('SE', gwas2.columns) if SE2 is None else SE2

    if SE2 is None and MAF2 and Ncases2 and Ncontrols2:
        SE2 = 'SE2'
        gwas2['Prop'] = gwas2[Ncases2] /(gwas2[Ncontrols2]+gwas2[Ncases2])
        gwas2['MAF_term'] = gwas2['MAF2'] * (1 - gwas2['MAF2'])
        gwas2[SE2] = 1 / np.sqrt(2 * gwas2['MAF_term'] * gwas2[Ncases2]*gwas2['Prop']*(1-gwas2['Prop']))


    elif SE2 is None and MINOR2 and MAJOR2 and Ncases2 and Ncontrols2:
        print(f"Calculating SE with MAF from MAFdata file for {gwas1_filename}.")
        MAFdata = pd.read_csv("MAFdata/OUD_stringent.EUR.MVP.DrugAlcDep2021_formatted_cleaned.txt", sep="\t")

        gwas2[MINOR2] = gwas2[MINOR2].str.upper()
        gwas2[MAJOR2] = gwas2[MAJOR2].str.upper()

        # Merge gwas1 with MAFdata to avoid row-wise operations
        merged_data = gwas2.merge(
        MAFdata[['rsids', 'alt', 'ref', 'af_alt']],
        left_on=[SNP2, MINOR2, MAJOR2],
        right_on=['rsids', 'alt', 'ref'],
        how='inner'
        )
        # Calculate N and SE
        merged_data['Prop'] = merged_data[Ncases2] /(merged_data[Ncontrols2]+merged_data[Ncases2])
        merged_data['MAF_term'] = merged_data['af_alt'] * (1 - merged_data['af_alt'])
        SE2 = 'SE2'
        merged_data[SE2] = 1 / np.sqrt(2 * merged_data['MAF_term'] * merged_data[Ncases2]*merged_data['Prop']*(1-merged_data['Prop']))
        gwas2 = merged_data  # Update gwas1 with the merged data    

    elif SE2 is None:
        #Warning for missing SE1
        print(f"Warning: SE is not specified for {gwas2_filename}. Please provide # cases/ controls and effect/ non=-ffect allele to calculate SE. If you do not have these, please provide SE column in the GWAS file.")
    
    BETA1 = col_lookup('BETA', gwas1.columns) if BETA1 is None or BETA1.strip() == "" else BETA1
    OR1 = col_lookup('OR', gwas1.columns) if OR1 is None or OR1.strip() == "" else OR1
    if (BETA1 is None or BETA1.strip() == "") and (OR1 or col_lookup('OR', gwas1.columns)):
        BETA1 = 'BETA1'
        gwas1[BETA1] = gwas1[col_lookup('OR', gwas1.columns)].apply(lambda x: pd.np.log(x))
    elif BETA1 is None and Z1 and SE1:
        BETA1 = 'BETA1'
        gwas1['BETA1'] = gwas1[Z1] * gwas1[SE1]
    elif BETA1 is None:
        #Warning for missing 
        print(f"Warning: BETA is not specified for {gwas1_filename}. Provide Beta, Odds ratio or Z score.")

    BETA2 = col_lookup('BETA', gwas2.columns) if BETA2 is None else BETA2
    OR2 = col_lookup('OR', gwas2.columns) if OR2 is None or OR2.strip() == "" else OR2
    if (BETA2 is None or BETA2.strip() == "") and (OR2 or col_lookup('OR', gwas2.columns)):
        BETA2 = 'BETA2'
        gwas2[BETA2] = gwas2[col_lookup('OR', gwas2.columns)].apply(lambda x: pd.np.log(x))
    elif BETA2 is None and Z2 and SE2:
        BETA2 = 'BETA2'
        gwas2['BETA2'] = gwas2[Z2] * gwas2[SE2]
    elif BETA2 is None:
        #Warning for missing 
        print(f"Warning: BETA is not specified for {gwas2_filename}. Provide Beta, Odds ratio or Z score.")

    print(f"Using columns: {SNP1}, {SNP2}, {BETA1}, {BETA2}, {SE1}, {SE2}, {MINOR1}, {MINOR2}, {MAJOR1}, {MAJOR2}")
    gwas1 = gwas1.rename(columns={SNP1: 'SNP1', MINOR1: 'MINOR1', MAJOR1: 'MAJOR1', 
                                    BETA1: 'BETA1', SE1: 'SE1'})
    gwas2 = gwas2.rename(columns={SNP2: 'SNP2', MINOR2: 'MINOR2', MAJOR2: 'MAJOR2', 
                                    BETA2: 'BETA2', SE2: 'SE2'})
    gwas1 = gwas1.rename(columns={col_lookup('CHR', gwas1.columns): 'CHR1', col_lookup('POS', gwas1.columns): 'POS1'})
    gwas2 = gwas2.rename(columns={col_lookup('CHR', gwas2.columns): 'CHR2', col_lookup('POS', gwas2.columns): 'POS2'})

    if MINOR1 and MINOR2 and MAJOR1 and MAJOR2:
        merged_1 = pd.merge(gwas1, gwas2, left_on=['SNP1', 'MINOR1', 'MAJOR1'], right_on=['SNP2', 'MINOR2', 'MAJOR2'])
        merged_2 = pd.merge(gwas1, gwas2, left_on=['SNP1', 'MINOR1', 'MAJOR1'], right_on=['SNP2', 'MAJOR2', 'MINOR2'])
        merged_2['BETA2'] = -1 * merged_2['BETA2']
        merged = pd.concat([merged_1, merged_2], ignore_index=True).drop_duplicates(subset=['SNP1'])
    
    else:
        if not MINOR1 or MINOR1.strip() == "":
            MINOR1 = col_lookup('ALLELE', gwas1.columns)
        if not MINOR2 or MINOR2.strip() == "":
            MINOR2 = col_lookup('ALLELE', gwas2.columns)
        if not MAJOR1 or MAJOR1.strip() == "":
            MAJOR1 = col_lookup('ALLELE', [col for col in gwas1.columns if col != MINOR1])
        if not MAJOR2 or  MAJOR2.strip() == "":
            MAJOR2 = col_lookup('ALLELE', [col for col in gwas2.columns if col != MINOR2])

        gwas1 = gwas1.rename(columns={MINOR1: 'MINOR1', MAJOR1: 'MAJOR1'})
        gwas2 = gwas2.rename(columns={MINOR2: 'MINOR2', MAJOR2: 'MAJOR2'})

        gwas1['MINOR1'] = gwas1['MINOR1'].str.upper()
        gwas1['MAJOR1'] = gwas1['MAJOR1'].str.upper()
        gwas2['MINOR2'] = gwas2['MINOR2'].str.upper()
        gwas2['MAJOR2'] = gwas2['MAJOR2'].str.upper()

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
            'merge_flip_neg': merge_flip_neg.dropna(subset=['BETA1', 'BETA2'])['BETA1'].corr(merge_flip_neg.dropna(subset=['BETA1', 'BETA2'])['BETA2']),
            'merge_flip': merge_flip.dropna(subset=['BETA1', 'BETA2'])['BETA1'].corr(merge_flip.dropna(subset=['BETA1', 'BETA2'])['BETA2']),
            'merge_neg': merge_neg.dropna(subset=['BETA1', 'BETA2'])['BETA1'].corr(merge_neg.dropna(subset=['BETA1', 'BETA2'])['BETA2']),
            'merge': merge.dropna(subset=['BETA1', 'BETA2'])['BETA1'].corr(merge.dropna(subset=['BETA1', 'BETA2'])['BETA2'])
        }
        # Filter based on correlation conditions
        selected_merges = []
        for key, corr in correlations.items():
            if corr > 0:
                selected_merges.append(locals()[key])
            if len(selected_merges) == 2:  # Stop after selecting two merges
                break
        
        if len(selected_merges) < 2:
           #append the one where correlation is nan
            for df in [merge_flip_neg, merge_flip, merge_neg, merge]:   
                if isinstance(df, pd.DataFrame) and df.shape[0] == 0:
                    selected_merges.append(df)
                    if len(selected_merges) == 2:
                        break
        # Concatenate the selected merges
        merged = pd.concat(selected_merges, ignore_index=True).drop_duplicates(subset=['SNP1'], keep='first')
        #warn the user that there are duplicates in SNP1 and check data quality
        if any('SNP1' in df.columns and df['SNP1'].duplicated().any() for df in selected_merges):
            print(f"Warning: There are duplicates in SNP1 after merging {gwas1_filename} and {gwas2_filename}. Please check the data quality.")

    if 'CHR1' and 'POS1' in merged.columns:
        phemed_out = merged[['SNP1', 'CHR1', 'POS1', 'BETA1', 'BETA2', 'SE1', 'SE2']]
        phemed_out = phemed_out.rename(columns={
            'SNP1': 'SNP',
            'CHR1': 'CHR',
            'POS1': 'POS',
            'BETA1': 'STUDY'+ gwas1_filename,
            'BETA2': 'STUDY'+gwas2_filename,
            'SE1': 'SE'+gwas1_filename,
            'SE2': 'SE'+gwas2_filename
        })
    elif 'CHR1' not in merged.columns or 'POS1' not in merged.columns:
        phemed_out = merged[['SNP1', 'CHR2', 'POS2', 'BETA1', 'BETA2', 'SE1', 'SE2']]

        phemed_out = phemed_out.rename(columns={
            'SNP1': 'SNP',
            'CHR2': 'CHR',
            'POS2': 'POS',
            'BETA1': 'STUDY'+ gwas1_filename,
            'BETA2': 'STUDY'+gwas2_filename,
            'SE1': 'SE'+gwas1_filename,
            'SE2': 'SE'+gwas2_filename
        })

    else:
        raise ValueError("CHR or POS columns not found in the merged DataFrame. Please check the input files.")

    return phemed_out


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Merge GWAS summary statistics files.")
    parser.add_argument("--gwas_files", required=True, nargs='+', help="Paths to GWAS summary statistics files. Please ensure that the first GWAS has SNP, chromosome, position columns. Space-separated.")
    parser.add_argument("--SNP", nargs='+', help="SNP column names. Use None to skip specifying this argument for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--BETA", nargs='+', help="BETA column names. Use None to skip specifying this argument for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--OR", nargs='+', help="Odds Ratio column names. Use None to skip specifying this argument for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--SE", nargs='+', help="SE column names. Use None to skip specifying this argument for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--MINOR", nargs='+', help="MINOR allele column names. Use None to skip specifying this argument for a GWAS file. Space-separated." \
    "give none or both or minor and major.", default=[])
    parser.add_argument("--MAJOR", nargs='+', help="MAJOR allele column names. Use None to skip specifying this argument for a GWAS file. Space-separated." \
    "give none or both or minor and major.", default=[])
    parser.add_argument("--Z", nargs='+', help="Z Score column names. Use None to skip specifying this argument for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--MAF", nargs='+', help="Minor allele frequency column names. Use None to skip specifying this argument for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--Ncases", nargs='+', help="Number of cases column names. Use None to skip specifying this argument " \
    "for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--Ncontrols", nargs='+', help="Number of controls column names. Use None to skip specifying this " \
    "argument for a GWAS file. Space-separated.", default=[])
    parser.add_argument("--output", help="Output file path for merged results.", default="merged_gwas.csv")

    args = parser.parse_args()

    # Load GWAS files
    gwas_files = args.gwas_files
    SNP = (
    [None if x == 'None' else x for x in args.SNP]
    + [None] * (len(gwas_files) - len(args.SNP))
    if len(args.SNP) < len(gwas_files)
    else [None if x == 'None' else x for x in args.SNP]
)
    BETA = (
    [None if x == 'None' else x for x in args.BETA]
    + [None] * (len(gwas_files) - len(args.BETA))
    if len(args.BETA) < len(gwas_files)
    else [None if x == 'None' else x for x in args.BETA]
)
    SE = (
    [None if x == 'None' else x for x in args.SE]
    + [None] * (len(gwas_files) - len(args.SE))
    if len(args.SE) < len(gwas_files)
    else [None if x == 'None' else x for x in args.SE]
)
    MINOR = (
    [None if x == 'None' else x for x in args.MINOR]
    + [None] * (len(gwas_files) - len(args.MINOR))
    if len(args.MINOR) < len(gwas_files)
    else [None if x == 'None' else x for x in args.MINOR]
)
    MAJOR = (
    [None if x == 'None' else x for x in args.MAJOR]
    + [None] * (len(gwas_files) - len(args.MAJOR))
    if len(args.MAJOR) < len(gwas_files)
    else [None if x == 'None' else x for x in args.MAJOR]
)
    Z = args.Z + [None] * (len(gwas_files) - len(args.Z)) if len(args.Z) < len(gwas_files) else args.Z
    MAF = args.MAF + [None] * (len(gwas_files) - len(args.MAF)) if len(args.MAF) < len(gwas_files) else args.MAF
    Ncases = args.Ncases + [None] * (len(gwas_files) - len(args.Ncases)) if len(args.Ncases) < len(gwas_files) else args.Ncases
    Ncontrols = args.Ncontrols + [None] * (len(gwas_files) - len(args.Ncontrols)) if len(args.Ncontrols) < len(gwas_files) else args.Ncontrols

    OR = (
    [None if x == 'None' else x for x in args.OR]
    + [None] * (len(gwas_files) - len(args.OR))
    if len(args.OR) < len(gwas_files)
    else [None if x == 'None' else x for x in args.OR]
)

    OUTPUT = args.output
    # Iteratively merge GWAS files
    merged_list = []
    for i in range(len(gwas_files)):
        for j in range(i + 1, len(gwas_files)):
            merged = gwas_merge(
                gwas_files[i], gwas_files[j],
                SNP[i], SNP[j],
                BETA[i], BETA[j],
                SE[i], SE[j],
                MINOR[i], MINOR[j],
                MAJOR[i], MAJOR[j],
                Z[i], Z[j],
                MAF[i], MAF[j],
                Ncases[i], Ncases[j],
                Ncontrols[i], Ncontrols[j],
                OR[i], OR[j]
            )
            merged_list.append(merged)
    
    if len(merged_list) == 1:
        result = merged_list[0]
        #Rename result cols as SNP, CHR, POS, STUDY1, STUDY2, SE1, SE2
        result.columns = [
            'SNP',
            'CHR',
            'POS',
            'STUDY1',
            'STUDY2',
            'SE1',
            'SE2'
        ]
    else:
        result = merge_of_merge(merged_list)

    #Save the merged output
    result.to_csv(OUTPUT, index=False)

if __name__ == "__main__":
    main()