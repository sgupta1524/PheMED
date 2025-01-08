import pandas as pd
import numpy as np
import argparse
import warnings

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
    'se_effect': 'SE'
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

def process_chunk(chunk, idx):
    """
    Process a chunk of data, converting OR to BETA if necessary and renaming columns.
    """
    chunk.columns = [col.strip().lower() for col in chunk.columns]  # Convert header to lowercase for case-insensitive checks

    # Get the required columns
    snp_col = get_column_name(chunk.columns, 'SNP')
    chrom_col = get_column_name(chunk.columns, 'CHR')
    pos_col = get_column_name(chunk.columns, 'POS')
    ref_col = get_column_name(chunk.columns, 'A1')
    alt_col = get_column_name(chunk.columns, 'A2')
    beta_col = get_column_name(chunk.columns, 'BETA')
    or_col = get_column_name(chunk.columns, 'OR')
    se_col = get_column_name(chunk.columns, 'SE')

    if beta_col is None and or_col is None:
        raise ValueError("Either 'beta' or 'OR' column must exist in the input files.")

    # Convert OR to BETA if necessary
    if or_col is not None:
        chunk[beta_col] = np.log(chunk[or_col].astype(float))

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

    return chunk[['SNP', 'CHR', 'POS', 'REF', 'ALT', 'BETA{}'.format(idx + 1), 'SE{}'.format(idx + 1)]]

def parse_dat(inputs):
    """
    Read input summary stats files in chunks, compare SNP, REF, and ALT columns,
    and merge the files to include SNP, CHROM, POS, BETA1, SE1, BETA2, SE2 columns.
    """
    input_files = inputs.split(',')
    data_frames = []
    snp_sets = []

    for idx, file in enumerate(input_files):
        chunks = pd.read_csv(file, sep=r'\s+|,|\t', engine='python', chunksize=100000)
        processed_chunks = [process_chunk(chunk, idx) for chunk in chunks]
        df = pd.concat(processed_chunks, ignore_index=True)
        data_frames.append(df)
        snp_sets.append(set(df['SNP']))

        # Debug: Print the DataFrame before merging
        print("DataFrame from file {} before merging:".format(file))
        print(df)

    # Merge data frames on SNP, CHR, POS, REF, and ALT using outer join
    merged_df = data_frames[0]
    for df in data_frames[1:]:
        merged_df = pd.merge(merged_df, df, on=['SNP', 'CHR', 'POS', 'REF', 'ALT'], how='outer')

    # Debug: Print the shape of the merged dataframe before filtering
    print("Merged DataFrame shape before filtering:", merged_df.shape)
    print(merged_df)

    # Count SNPs based on their presence in the input files
    snp_counts = {}
    for snp_set in snp_sets:
        for snp in snp_set:
            if snp in snp_counts:
                snp_counts[snp] += 1
            else:
                snp_counts[snp] = 1

    # Filter SNPs that appear in at least two files
    snps_to_keep = [snp for snp, count in snp_counts.items() if count >= 2]
    merged_df = merged_df[merged_df['SNP'].isin(snps_to_keep)]

    # Debug: Print the shape of the merged dataframe after filtering
    print("Merged DataFrame shape after filtering:", merged_df.shape)
    print(merged_df)

    # Sort the merged dataframe by SNP, CHR, REF, and ALT
    merged_df.sort_values(by=['SNP', 'CHR', 'REF', 'ALT'], inplace=True)

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
    parser.add_argument("--log", type=str, help="path to log file")
    args = parser.parse_args()

    # Save the inputs to variables
    num_summary_stats = args.n_files
    input_files = args.inputs
    output_file = args.output
    log_file = args.log

    # Call parse_dat with input_files
    merged_data = parse_dat(input_files)
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
        warnings.warn("The merged summary stats file has less than 100 rows. The GWAS stats might be from different builds.")


    # Write the merged data to the output file
    merged_data.to_csv(output_file, index=False)

    # Write the log file
    with open(log_file, 'w') as log:
        log.write('\n'.join(log_entries))

if __name__ == '__main__':
    merge_summary_stats()