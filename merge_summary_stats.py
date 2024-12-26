import csv
import numpy as np
import argparse
import re

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
    'hg18chrom': 'CHR',
    'hg36chr': 'CHR',
    'hg36chrom': 'CHR',
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

def get_column_index(header, target):
    """
    Get the index of the target column in the header row.
    The header row is converted to lowercase to make the search case insensitive.
    """
    header = [col.lower() for col in header]
    for col_name, mapped_name in COLUMN_MAP.items():
        if mapped_name == target and col_name in header:
            return header.index(col_name)
    raise ValueError("Column for {} not found in header".format(target))

def parse_dat(inputs):
    """
    Read input summary stats files in chunks, compare SNP, REF, and ALT columns,
    and merge the files to include SNP, CHROM, POS, BETA1, SE1, BETA2, SE2 columns.
    """
    input_files = inputs.split(',')
    data = {}

    for idx, file in enumerate(input_files):
        with open(file, 'r') as f:
            # Determine the delimiter based on the first line of the file
            first_line = f.readline()
            if '\t' in first_line:
                delimiter = '\t'
            elif ',' in first_line:
                delimiter = ','
            else:
                delimiter = None  # Use regex for space delimiter

            f.seek(0)  # Reset file pointer to the beginning

            if delimiter:
                reader = csv.reader(f, delimiter=delimiter)
            else:
                reader = (re.split(r'\s+', line.strip()) for line in f)

            header = next(reader)  # Read the header row
            header = [col.strip() for col in header]

            # Get the indices of the required columns
            snp_index = get_column_index(header, 'SNP')
            chrom_index = get_column_index(header, 'CHR')
            pos_index = get_column_index(header, 'POS')
            ref_index = get_column_index(header, 'A1')
            alt_index = get_column_index(header, 'A2')
            beta_index = get_column_index(header, 'BETA') if idx == 0 else get_column_index(header, 'OR')
            se_index = get_column_index(header, 'SE')

            for row in reader:
                snp = row[snp_index]
                ref = row[ref_index]
                alt = row[alt_index]
                chrom = row[chrom_index]
                pos = row[pos_index]
                beta = row[beta_index]
                se = row[se_index]

                if snp not in data:
                    data[snp] = {
                        'CHROM': chrom,
                        'POS': pos,
                        'REF': ref,
                        'ALT': alt,
                        'BETA1': None,
                        'SE1': None,
                        'BETA2': None,
                        'SE2': None
                    }

                # Check if REF and ALT match, then update the data dictionary
                if ref == data[snp]['REF'] and alt == data[snp]['ALT']:
                    if idx == 0:
                        data[snp]['BETA1'] = beta
                        data[snp]['SE1'] = se
                    elif idx == 1:
                        data[snp]['BETA2'] = np.log(float(beta))
                        data[snp]['SE2'] = se

    return data

def filter_alleles(a):
    '''Remove alleles that do not describe strand-unambiguous SNPs'''

def merge_summary_stats():
    '''
    Parse and munge summary statistics
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument("--n_files", type=str, help="Number of summary statistics to munge")
    parser.add_argument("--inputs", type=str, help="comma separated paths to input summary statistics files")
    parser.add_argument("--output", type=str, help="path to output munged summary statistics file")
    args = parser.parse_args()

    # Save the inputs to variables
    num_summary_stats = args.n_files
    input_files = args.inputs
    output_file = args.output

    # Call parse_dat with input_files
    merged_data = parse_dat(input_files)

    # Write the merged data to the output file
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['SNP', 'CHR', 'POS', 'BETA1', 'SE1', 'BETA2', 'SE2'])
        for snp, values in merged_data.items():
            writer.writerow([snp, values['CHROM'], values['POS'], values['BETA1'], values['SE1'], values['BETA2'], values['SE2']])

if __name__ == '__main__':
    merge_summary_stats()