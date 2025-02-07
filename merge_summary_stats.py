import logging
import pandas as pd
import numpy as np
import argparse
import warnings


def setup_logging(log_file):
    # Create the root logger
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)  # Adjust the level to WARNING or lower

    # File handler for logging to a file
    file_handler = logging.FileHandler(log_file, mode='a')  # Append to the log file
    file_handler.setLevel(logging.WARNING)
    file_formatter = logging.Formatter('\n%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)

    # Console handler for logging to stdout
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.WARNING)
    console_formatter = logging.Formatter('\n%(levelname)s - %(message)s')  # Simpler format for console
    console_handler.setFormatter(console_formatter)

    # Add both handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)


# Mapping of various column names to standardized names
COLUMN_MAP = {
    'snp': 'SNP',
    'markername': 'SNP',
    'snpid': 'SNP',
    'rs': 'SNP',
    'rsid': 'SNP',
    'rsids': 'SNP',
    'rs_number': 'SNP',
    'rs_numbers': 'SNP',
    'chr': 'CHR',
    'chrom': 'CHR',
    'chromosome': 'CHR',
    'chromosome_number': 'CHR',
    'hg18chr': 'CHR',
    'hg19chr': 'CHR',
    '#chr': 'CHR',
    '#chrom': 'CHR',
    '#Chrom': 'CHR',
    '#Chr': 'CHR',
    'bp': 'POS',
    'position': 'POS',
    'pos': 'POS',
    'a1': 'A1',
    'allele1': 'A1',
    'allele_1': 'A1',
    'effect_allele': 'A1',
    'reference_allele': 'A1',
    'inc_allele': 'A1',
    'ref': 'A1',
    'ea': 'A1',
    'a2': 'A2',
    'allele2': 'A2',
    'allele_2': 'A2',
    'other_allele': 'A2',
    'non_effect_allele': 'A2',
    'dec_allele': 'A2',
    'nea': 'A2',
    'alt': 'A2',
    'beta': 'BETA',
    'effect': 'BETA',
    'effects': 'BETA',
    'log_odds': 'BETA',
    'or': 'OR',
    'se': 'SE',
    'standard_error': 'SE',
    'standard_errors': 'SE',
    'se_beta': 'SE',
    'se_or': 'SE',
    'se_log_odds': 'SE',
    'se_effects': 'SE',
    'se_effect': 'SE',
    'sebeta': 'SE'
}

def get_column_name(header, target):
    """
    Get the standardized column name from the header row.
    The header row is converted to lowercase to make the search case insensitive.
    """
    header_lower = [col.lower() for col in header]
    for col_name, mapped_name in COLUMN_MAP.items():
        if mapped_name == target and col_name in header_lower:
            return col_name
    return None  # Return None if the target column is not found

def process_chunk(chunk, idx, effect_allele_col, non_effect_allele_col):
    """
    Process a chunk of data, converting OR to BETA if necessary and renaming columns.
    """
    chunk.columns = [col.strip().lower() for col in chunk.columns]  # Convert header to lowercase for case-insensitive checks
    #print(chunk.columns)
    #print(non_effect_allele_col)
    # Get the required columns
    snp_col = get_column_name(chunk.columns, 'SNP')
    chrom_col = get_column_name(chunk.columns, 'CHR')
    pos_col = get_column_name(chunk.columns, 'POS')
    if effect_allele_col == "":
        ref_col = get_column_name(chunk.columns, "A1")
    else:
        ref_col = effect_allele_col.lower()
    if non_effect_allele_col == "":
        alt_col = get_column_name(chunk.columns, "A2")
    else:
        alt_col = non_effect_allele_col.lower()
    beta_col = get_column_name(chunk.columns, 'BETA')
    or_col = get_column_name(chunk.columns, 'OR')
    se_col = get_column_name(chunk.columns, 'SE')

    #print(non_effect_allele_col)
    #print(effect_allele_col,ref_col)
    #print(non_effect_allele_col,alt_col)
    # print(chunk.columns)
    # print(type(snp_col))
    # print(type(chrom_col))
    # print(type(beta_col))
    # print(type(pos_col))
    # print(type(ref_col))
    # print(type(alt_col))
    # print(type(or_col))

    if beta_col is None and or_col is None:
        raise ValueError("Either 'beta' or 'OR' column must exist in the input files.")

    # Convert OR to BETA if necessary
    if or_col is not None:
        chunk[or_col] = pd.to_numeric(chunk[or_col], errors='coerce')
        chunk[beta_col] = np.log(chunk[or_col])

    #print(chunk)
    # Rename columns to standardized names
    chunk.rename(columns={
        snp_col: 'SNP',
        chrom_col: 'CHR',
        pos_col: 'POS',
        ref_col: 'REF',
        alt_col: 'ALT',
        beta_col: 'BETA{}'.format(idx + 1),
        se_col: 'SE{}'.format(idx + 1)
    }, inplace=True)
    #print(snp_col,ref_col, alt_col)
    #print(chunk)

    return chunk[['SNP', 'CHR', 'POS', 'REF', 'ALT', 'BETA{}'.format(idx + 1), 'SE{}'.format(idx + 1)]]

def parse_dat(inputs, effect_allele_cols_list, non_effect_allele_cols_list):
    """
    Read input summary stats files in chunks, compare SNP, REF, and ALT columns,
    and merge the files to include SNP, CHROM, POS, BETA1, SE1, BETA2, SE2 columns.
    """
    input_files = inputs.split(',')
    #print(input_files)
    data_frames = []
    snp_sets = []

    for idx, (file, effect_allele_col, non_effect_allele_col) in enumerate(zip(input_files, effect_allele_cols_list, non_effect_allele_cols_list)):
        logging.info(f"Processing file {file}")
        print(f"Processing file {file}")
        # Read the first few lines to determine the separator
        with open(file, 'r') as f:
            first_line = f.readline()
            if ',' in first_line:
                sep = ','
            elif '\t' in first_line:
                sep = '\t'
            else:
                sep = r'\s+'

        chunks = pd.read_csv(file, sep=sep, engine='python', chunksize=10, comment='##')
        #print("Effect allele column:", effect_allele_col)
        #print("Non-effect allele column:", non_effect_allele_col)
        processed_chunks = [process_chunk(chunk, idx, effect_allele_col, non_effect_allele_col) for chunk in chunks]
        #print(processed_chunks)
        df = pd.concat(processed_chunks, ignore_index=True)
        data_frames.append(df)
        # Create a combination of SNP, CHR, and POS
        SNP_CHR_POS = df['SNP'].astype(str) + '_' + df['CHR'].astype(str) + '_' + df['POS'].astype(str)
        snp_sets.append(set(SNP_CHR_POS))

        # Debug: Print the DataFrame before merging
        #print("DataFrame from file {} before merging:".format(file))
        #print(df)

    # Merge data frames on SNP, CHR, POS, REF, and ALT using outer join
    if not data_frames:
        raise ValueError("No data frames were created. Please check the input files and columns.")
    
    merged_df = data_frames[0]
    merged_df = merged_df.astype(str)
    for df in data_frames[1:]:
        df = df.astype(str)
        merged_df = pd.merge(merged_df, df, on=['SNP', 'CHR', 'POS', 'REF', 'ALT'], how='outer')
        print(merged_df)
    # Debug: Print the shape of the merged dataframe before filtering
    # print("Merged DataFrame shape before filtering:", merged_df.shape)
    # print(merged_df)

    # Count SNPs based on their presence in the input files
    snp_counts = {}
    for snp_set in snp_sets:
        for snp in snp_set:
            if snp in snp_counts:
                snp_counts[snp] += 1
            else:
                snp_counts[snp] = 1

    #TODO: Merge process can be made faster
    # Filter SNPs that appear in at least two files
    snps_to_keep = [snp.split('_')[0] for snp, count in snp_counts.items() if count >= 2]
    merged_df = merged_df[merged_df['SNP'].isin(snps_to_keep)]

    # Debug: Print the shape of the merged dataframe after filtering
    #print("Merged DataFrame shape after filtering:", merged_df.shape)
    #print(merged_df)

    # Sort the merged dataframe by SNP, CHR, REF, and ALT
    merged_df.sort_values(by=['SNP', 'CHR', 'REF', 'ALT'], inplace=True)
    # Free up memory by deleting other data frames
    del data_frames
    # Replace empty values with NA
    merged_df.fillna('NA', inplace=True)

    # Log a warning if the merged DataFrame is empty or has less than 100 rows
    if len(merged_df) < 100:
        logging.warning("The merged summary stats file has less than 100 rows. The GWAS stats might be from different reference genome builds.")

    return merged_df

def merge_summary_stats():
    '''
    Parse and munge summary statistics
    '''

    parser = argparse.ArgumentParser()

    #remove n-files
    parser.add_argument("--inputs", type=str, help="comma separated paths to input summary statistics files")
    parser.add_argument("--output", type=str, help="path to output munged summary statistics file", default="output/merged_summary_stats.csv")
    parser.add_argument("--log", type=str, help="path to log file", default="output/log")
    parser.add_argument("--effect-allele-col", type=str, help="comma separated effect allele columns per file", required = False, default="")
    parser.add_argument("--non-effect-allele-col", type=str, help="comma separated non effect allele columns per file", required = False, default="")
    args = parser.parse_args()

    # Save the inputs to variables
    input_files = args.inputs
    output_file = args.output
    log_file = args.log
    num_summary_stats = len(input_files.split(",")) 
    effect_allele_cols = args.effect_allele_col
    non_effect_allele_cols = args.non_effect_allele_col

    num_summary_stats = len(input_files)  # Number of summary stats files

    # Ensure effect_allele_cols and non_effect_allele_cols are of length num_summary_stats
    if len(effect_allele_cols.split(",")) < num_summary_stats:
        effect_allele_cols_list = effect_allele_cols.split(",")
        effect_allele_cols_list += [""] * (num_summary_stats - len(effect_allele_cols_list))
    
    if len(non_effect_allele_cols.split(",")) < num_summary_stats:
        non_effect_allele_cols_list = non_effect_allele_cols.split(",")
        non_effect_allele_cols_list += [""] * (num_summary_stats - len(non_effect_allele_cols_list))

    # Setup logging to file
    setup_logging(log_file)

    # Call parse_dat with input_files
    merged_data = parse_dat(input_files, effect_allele_cols_list, non_effect_allele_cols_list)

    log_entries = []
    # Select only the required columns
    columns_to_keep = ['SNP', 'CHR', 'POS']
    for i, file in enumerate(input_files.split(',')):
        columns_to_keep.append('BETA{}'.format(i + 1))
        columns_to_keep.append('SE{}'.format(i + 1))
        log_entries.append(f'BETA{i + 1} and SE{i + 1} come from {file}')

    merged_data = merged_data[columns_to_keep]

    # Write the merged data to the output file
    merged_data.to_csv(output_file, index=False)
    with open(log_file, 'a') as log:
        log.write('\n'.join(log_entries))

if __name__ == '__main__':
    merge_summary_stats()