import pandas as pd
import numpy as np
import argparse
import warnings
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Mapping of various column names to standardized names
COLUMN_MAP = {
    'snp': 'SNP',
    'markername': 'SNP',
    'snpid': 'SNP',
    'rs': 'SNP',
    'rsid': 'SNP',
    'rs_number': 'SNP',
    'rs_numbers': 'SNP',
    'chr': 'CHR',
    'chrom': 'CHR',
    'chromosome': 'CHR',
    'chromosome_number': 'CHR',
    'hg18chr': 'CHR',
    'hg19chr': 'CHR',
    'bp': 'POS',
    'position': 'POS',
    'pos': 'POS',
    'a1': 'A1',
    'allele1': 'A1',
    'allele_1': 'A1',
    'effect_allele': 'A1',
    'reference_allele': 'A1',
    'inc_allele': 'A1',
    'ea': 'A1',
    'a2': 'A2',
    'allele2': 'A2',
    'allele_2': 'A2',
    'other_allele': 'A2',
    'non_effect_allele': 'A2',
    'dec_allele': 'A2',
    'nea': 'A2',
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
    'effect_allele_frequency': 'MAF',
    'eaf': 'MAF',
    'maf': 'MAF',
    'minor_allele_frequency': 'MAF',
    'n': 'N',
    'sample_size': 'N',
    'sample': 'N',
    'samples': 'N',
    'sample_size_per_variant': 'N',
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

def process_chunk(chunk, idx, effect_allele_col, non_effect_allele_col, n_col=None, maf_col=None):
    """
    Process a chunk of data, converting OR to BETA if necessary and renaming columns.
    """
    chunk.columns = [col.strip().lower() for col in chunk.columns]  # Convert header to lowercase for case-insensitive checks

    # Get the required columns
    snp_col = get_column_name(chunk.columns, 'SNP')
    chrom_col = get_column_name(chunk.columns, 'CHR')
    pos_col = get_column_name(chunk.columns, 'POS')
    ref_col = effect_allele_col.lower()  # Use the user-specified effect allele column as A1
    alt_col = non_effect_allele_col.lower()  # Use the user-specified non-effect allele column as A2
    beta_col = get_column_name(chunk.columns, 'BETA')
    or_col = get_column_name(chunk.columns, 'OR')
    se_col = get_column_name(chunk.columns, 'SE')
    n_col = get_column_name(chunk.columns, 'N')
    maf_col = get_column_name(chunk.columns, 'MAF')

    # Check if required columns are present
    if not all([snp_col, chrom_col, pos_col, ref_col, alt_col]):
        raise ValueError("One or more required columns (SNP, CHR, POS, REF, ALT) are missing in the input files.")

    if beta_col is None and or_col is None:
        raise ValueError("Either 'beta' or 'OR' column must exist in the input files.")

    # Convert OR to BETA if necessary
    if or_col is not None:
        chunk[beta_col] = np.log(chunk[or_col].astype(float))

    # Calculate SE if SE column is not present
    if se_col is None:
        if n_col is not None and maf_col is not None:
            chunk['SE'] = np.sqrt(1 / (2 * chunk[n_col].astype(float) * chunk[maf_col].astype(float) * (1 - chunk[maf_col].astype(float))))
            se_col = 'SE'
        else:
            raise ValueError("SE column is missing and required columns (N, MAF) to calculate SE are also missing.")

    # Rename columns to standardized names
    chunk.rename(columns={
        snp_col: 'SNP',
        chrom_col: 'CHR',
        pos_col: 'POS',
        ref_col: 'REF',  # Consider the effect allele column as A1
        alt_col: 'ALT',  # Consider the non-effect allele column as A2
        beta_col: 'BETA{}'.format(idx + 1),
        se_col: 'SE{}'.format(idx + 1)
    }, inplace=True)
    return chunk[['SNP', 'CHR', 'POS', 'REF', 'ALT', 'BETA{}'.format(idx + 1), 'SE{}'.format(idx + 1)]]

def parse_dat(inputs, effect_allele_cols, non_effect_allele_cols):
    """
    Read input summary stats files in chunks, compare SNP, REF, and ALT columns,
    and merge the files to include SNP, CHROM, POS, BETA1, SE1, BETA2, SE2 columns.
    """
    input_files = inputs.split(',')
    effect_allele_cols = effect_allele_cols.split(',') if effect_allele_cols else []
    non_effect_allele_cols = non_effect_allele_cols.split(',') if non_effect_allele_cols else []
    data_frames = []
    snp_sets = []

    for idx, file in enumerate(input_files):
        logging.info(f"Processing file {file}")
        chunks = pd.read_csv(file, sep=r'\s+|,|\t', engine='python', chunksize=100000)
        n_col = None
        maf_col = None
        for chunk in chunks:
            if n_col is None or maf_col is None:
                n_col = get_column_name(chunk.columns, 'N')
                maf_col = get_column_name(chunk.columns, 'MAF')
            
            processed_chunk = process_chunk(
                                            chunk, 
                                            idx, 
                                            effect_allele_cols[idx] if idx < len(effect_allele_cols) else "A1", 
                                            non_effect_allele_cols[idx] if idx < len(non_effect_allele_cols) else "A2", 
                                            n_col, 
                                            maf_col
                                        )
            data_frames.append(processed_chunk)
            snp_sets.append(set(processed_chunk['SNP']))

        # Debug: Print the DataFrame before merging
    print("DataFrame from file {} before merging:".format(file))
    #print(df)

    logging.info("Merging data frames")
    # Merge data frames on SNP, CHR, POS, REF, and ALT using outer join
    merged_df = data_frames[0]
    for df in data_frames[1:]:
        merged_df = pd.merge(merged_df, df, on=['SNP', 'CHR', 'POS', 'REF', 'ALT'], how='outer')

    # Debug: Print the shape of the merged dataframe before filtering
    print("Merged DataFrame shape before filtering:", merged_df.shape)
    print(merged_df)

    logging.info("Counting SNPs based on their presence in the input files")
    # Count SNPs based on their presence in the input files
    snp_counts = {}
    for snp_set in snp_sets:
        for snp in snp_set:
            if snp in snp_counts:
                snp_counts[snp] += 1
            else:
                snp_counts[snp] = 1

    logging.info("Filtering SNPs that appear in at least two files")
    # Filter SNPs that appear in at least two files
    snps_to_keep = [snp for snp, count in snp_counts.items() if count >= 2]
    merged_df = merged_df[merged_df['SNP'].isin(snps_to_keep)]

    # Debug: Print the shape of the merged dataframe after filtering
    print("Merged DataFrame shape after filtering:", merged_df.shape)
    print(merged_df)

    logging.info("Sorting the merged dataframe")
    # Sort the merged dataframe by SNP, CHR, REF, and ALT
    merged_df.sort_values(by=['SNP', 'CHR', 'REF', 'ALT'], inplace=True)

    logging.info("Replacing empty values with NA")
    # Replace empty values with NA
    merged_df.fillna('NA', inplace=True)

    return merged_df

def merge_summary_stats():
    '''
    Parse and munge summary statistics
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument("--n_files", type=str, help="Number of summary statistics to munge")
    parser.add_argument("--inputs", type=str, help="comma separated paths to input summary statistics files")
    parser.add_argument("--output", type=str, help="path to output munged summary statistics file")
    parser.add_argument("--log", type=str, help="path to log file", default = "log")
    parser.add_argument("--effect-allele-col", type=str, help="comma separated effect allele column", required = False)
    parser.add_argument("--non-effect-allele-col", type=str, help="comma separated non effect allele col", required = False)
    args = parser.parse_args()

    # Save the inputs to variables
    num_summary_stats = args.n_files
    input_files = args.inputs
    output_file = args.output
    log_file = args.log
    effect_allele_cols = args.effect_allele_col
    non_effect_allele_cols = args.non_effect_allele_col

    # Call parse_dat with input_files
    merged_data = parse_dat(input_files, effect_allele_cols, non_effect_allele_cols)
    log_entries = []
    # Select only the required columns
    columns_to_keep = ['SNP', 'CHR', 'POS']
    for i, file in enumerate(input_files.split(',')):
        columns_to_keep.append('BETA{}'.format(i + 1))
        columns_to_keep.append('SE{}'.format(i + 1))
        log_entries.append(f'BETA{i + 1} and SE{i + 1} come from {file}')

    merged_data = merged_data[columns_to_keep]

    # Check if the merged data has less than 100 rows and issue a warning
    if len(merged_data) < 100:
        warnings.warn("The merged summary stats file has less than 100 rows. The GWAS stats might be from different reference genome builds.")

    logging.info("Writing the merged data to the output file")
    # Write the merged data to the output file
    merged_data.to_csv(output_file, index=False)

    
    logging.info("Writing the log file")
    with open(log_file, 'w') as log:
        log.write('\n'.join(log_entries))

if __name__ == '__main__':
    merge_summary_stats()