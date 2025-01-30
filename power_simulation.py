import scipy.stats as stats
from scipy.stats import chi2
import math
import numpy as np
import matplotlib.pyplot as plt
import argparse
from math import exp, sqrt, pi
from scipy.special import erfinv
from scipy.special import erf

from scipy.stats import norm

def critical_value(alpha, tail='two-tailed'):
    """
    Calculate the critical value for a given alpha level and tail type.

    :param alpha: Significance level (e.g., 0.05)
    :param tail: 'one-tailed' or 'two-tailed' test
    :return: critical value corresponding to the alpha level
    """
    if tail == 'one-tailed':
        # For one-tailed test, calculate the critical value for alpha in the upper tail
        return norm.ppf(1 - alpha)
    elif tail == 'two-tailed':
        # For two-tailed test, calculate the critical value for alpha/2 in each tail
        return norm.ppf(1 - alpha / 2)
    else:
        raise ValueError("tail must be either 'one-tailed' or 'two-tailed'")

    # Example usage
    alpha = 0.05
    critical_val_one_tailed = critical_value(alpha, tail='one-tailed')
    critical_val_two_tailed = critical_value(alpha, tail='two-tailed')

    print(f"One-tailed critical value for alpha = {alpha}: {critical_val_one_tailed}")
    print(f"Two-tailed critical value for alpha = {alpha}: {critical_val_two_tailed}")

def benjamini_hochberg(pvals, alpha=0.05):
    """
    Perform Benjamini-Hochberg FDR correction.
    """
    pvals = np.array(pvals)
    n = len(pvals)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    adjusted_pvals = np.empty(n)
    cumulative_min = 1.0
    for i in range(n - 1, -1, -1):
        cumulative_min = min(cumulative_min, sorted_pvals[i] * n / (i + 1))
        adjusted_pvals[sorted_indices[i]] = cumulative_min
    return adjusted_pvals

def get_phi_value(phemed_summary_file, study):
    """
    read user input of phemed summary file and get phi for each study
    """
    return phi

def simulate_sample(num_simulations):
    """
    Simulate data for the power function.
    num_simulations = number of simulations to run
    """
    np.random.seed(42)  # For reproducibility

    # Generate random values for the parameters

    ##plausible values

    #from 1000 to 1 million
    # Sample size between 100 and 10000
    n_values = np.linspace(100, 100000, int(num_simulations)).astype(int)
    #beta between -.0.1 and 0.1
    #alpha 0.05 or 5*10e**-8
    # Minor allele frequency between 0.01 and 0.5
    # Dilution factor between 0.5 and 1.5 
    # P is typically 5% typically for neuropsychiatric traits

    return n_values

def calculate_power(n, beta, alpha, phi, P, tails, maf=None, SE=None):
    """
    Calculate the power of a GWAS study 
    n = sample size
    beta = effect size
    alpha = significance level
    phi = dilution factor
    P = proportion of cases
    maf = frequency of the risk allele (optional)
    SE = standard error (optional)
    """
    if maf is None and SE is None:
        raise ValueError("Either 'maf' or 'SE' must be provided.")
    
    # Calculate Standard error based on a logistic model if not provided
    if SE is None:
        SE = calculate_SE(n, maf, P)
    
    
    # Calculate the power of the study
    power = 1 - stats.norm.cdf(critical_value(alpha, tails) - np.divide(beta / float(phi), SE))
    print(power)
    
    # Return the power value per SNP
    return power

##Wineer's curse - bigdeli paper

##z-score=beta/SE
# z - association z-scores

# min.p - minimum p-value admitted (to avoid zero p-values/adjusted p-values which give troubles with inverse cdf)

# min.p - very large zs corresponding to min.p (i.e z > 37) are not adjusted, as their bias is essentially zero)
# FIQT <- function(z=z, min.p=10^-300){
#     pvals<-2*pnorm(abs(z),low=F)
#     pvals[pvals<min.p]<- min.p
#     adj.pvals<-p.adjust(pvals,method="fdr")
#     mu.z<-sign(z)*qnorm(adj.pvals/2,low=F)
#     mu.z[abs(z)>qnorm(min.p/2,low=F)]<-z[abs(z)>qnorm(min.p/2,low=F)]
#     mu.z
# }
#https://github.com/bacanusa/FIQT

def calculate_SE(n, maf, P):
    if P is None:
        raise ValueError("Proportion of cases (P) must be provided if SE column is not present.")
    SE = 1 / math.sqrt(2 * P   * (1 - P)*(n * (1 - maf) * maf))
    return SE

def FIQT(n,maf,P,beta, SE, min_p,alpha, tails):

    SE = calculate_SE(n, maf, P)
    # Calculate z-scores
    z = beta/SE
    # Calculate p-values
    pvals = tails * stats.norm.sf(np.abs(z))
    
    # Replace p-values less than min_p with min_p
    #Add this when dealing with range of pvalues
    #pvals[pvals < min_p] = min_p
    
    # Adjust p-values using the FDR method
    #Add this when dealing with range of pvalues
    #_, adj_pvals, _, _ = benjamini_hochberg(pvals, alpha)
    
    # Calculate adjusted z-scores
    mu_z = np.sign(z) * stats.norm.isf(pvals / tails)
    #print(mu_z)
    
    # Replace adjusted z-scores for very large z-scores
    threshold = stats.norm.isf(min_p / tails)
    if np.abs(z) > threshold:
        mu_z = z
    
    #Add this when dealing with range of pvalues
    #mu_z[np.abs(z) > threshold] = z[np.abs(z) > threshold]
    
    return mu_z

def calculate_corrected_power(n,maf,P,beta,phi,min_p_val,alpha,tails):
    """
    Calculate the power of a GWAS study 
    n = sample size
    beta = effect size
    alpha = significance level
    maf = frequency of the risk allele
    phi = dilution factor
    P = proportion of cases
    """
    if maf is None and SE is None:
        raise ValueError("Either 'maf' or 'SE' must be provided.")
    
    # Calculate Standard error based on a logistic model if not provided
    if SE is None:
        SE = calculate_SE(n, maf, P)
    
    z = beta/SE
    corrected_z = FIQT(n,maf,P,beta, SE, min_p_val,alpha, tails)
    NCP = (corrected_z/phi) ** 2
    # Calculate the power of the study
    power = 1 - stats.norm.cdf(critical_value(alpha, tails) - np.divide(beta / float(phi), SE))
    print(power)

def generate_plot(power, corrected_power, samples, beta, alpha, maf, phi, P, output_file):
    """
    Generate a plot for the given data.

    :param power: List or array of power values.
    :param corrected_power: List or array of corrected power values.
    :param samples: List or array of corresponding sample sizes or simulation indices.
    """
    if len(power) != len(samples) or len(corrected_power) != len(samples):
        raise ValueError("The lengths of 'power', 'corrected_power', and 'samples' must be the same.")
    
    plt.figure(figsize=(10, 6))
    plt.plot(samples, power, marker='o', linestyle='-', color='b', label='Power Curve', markersize=4)
    plt.plot(samples, corrected_power, marker='x', linestyle='-', color='g', label='Corrected Power Curve', markersize=4)
    plt.xlabel('Simulations or Sample Sizes', fontsize=12)
    plt.ylabel('Power', fontsize=12)
    plt.title('Power Simulation Plot', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.axhline(y=0.8, color='r', linestyle='--', label='80% Power Threshold')
    plt.legend()
    
    # Add a label on the bottom right to show the beta, alpha, maf, phi, P values
    textstr = f'beta={beta}, alpha={alpha}, maf={maf}, phi={phi}, P={P}'
    plt.gcf().text(0.95, 0.1, textstr, fontsize=10, verticalalignment='bottom', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def generate_weighted_p_thresholds(p_values=[], Power=[]):
    if len(p_values) != len(Power):
        raise ValueError("The lengths of 'p_values' and 'Power' must be the same.")
    
    mean_power = np.mean(Power)
    weighted_p = mean_power * np.array(p_values) / np.array(Power)
    return weighted_p

def calculate_power_sum_stats(summary_stats_file, alpha, P, tails, phi,n):
    '''
    Calculate power estimates from summary statistics file
    '''
    power = []
    with open(summary_stats_file, 'r') as file:
        lines = file.readlines()
    
    # Skip lines that start with '#'
    header_index = next(i for i, line in enumerate(lines) if not line.startswith('#'))
    header = lines[header_index].strip().split()
    beta_idx = next((i for i, col in enumerate(header) if col.lower() in ['beta', 'effect', 'effects', 'log_odds']), None)
    maf_idx = next((i for i, col in enumerate(header) if col.lower() in ['maf', 'minor_allele_frequency', 'allele_frequency']), None)
    se_idx = next((i for i, col in enumerate(header) if col.lower() in ['se', 'stderr', 'standard_error']), None)

    if beta_idx is None or (maf_idx is None and se_idx is None):
        raise ValueError("The summary statistics file must contain columns for beta and either MAF or SE.")
    
    for line in lines[1:]:  # Assuming the first line is a header
        columns = line.strip().split()
        beta = float(columns[beta_idx])
        maf = float(columns[maf_idx]) if maf_idx is not None else None

        if maf is not None:
            power.append(calculate_power(n, beta, alpha, maf, phi, P, tails))
        else:
            SE = float(columns[se_idx])
            power.append(calculate_power(n, beta, alpha, phi, P, tails, SE=SE))
    
    return power

def calculate_corrected_power_sum_stats(summary_stats_file, alpha, P, tails, phi,n):
    '''
    Calculate winner's curse corrected power estimates from summary statistics file
    '''
    power = []
    with open(summary_stats_file, 'r') as file:
        lines = file.readlines()
    
    # Skip lines that start with '#'
    header_index = next(i for i, line in enumerate(lines) if not line.startswith('#'))
    header = lines[header_index].strip().split()
    beta_idx = next((i for i, col in enumerate(header) if col.lower() in ['beta', 'effect', 'effects', 'log_odds']), None)
    maf_idx = next((i for i, col in enumerate(header) if col.lower() in ['maf', 'minor_allele_frequency', 'allele_frequency']), None)
    se_idx = next((i for i, col in enumerate(header) if col.lower() in ['se', 'stderr', 'standard_error']), None)

    if beta_idx is None or (maf_idx is None and se_idx is None):
        raise ValueError("The summary statistics file must contain columns for beta and either MAF or SE.")
    
    for line in lines[1:]:  # Assuming the first line is a header
        columns = line.strip().split()
        beta = float(columns[beta_idx])
        maf = float(columns[maf_idx]) if maf_idx is not None else None

        if maf is not None:
            power.append(calculate_corrected_power(n, beta, alpha, maf, phi, P, tails))
        else:
            SE = float(columns[se_idx])
            power.append(calculate_corrected_power(n, beta, alpha, phi, P, tails, SE=SE))
    
    return power

def power_simulation():
    '''
    Get power estimates and generate plots
    '''
    parser = argparse.ArgumentParser()

    #remove n-files
    parser.add_argument("--phemed-summary-file", type=str, help="Path to Phemed summary statistics file", required=False)
    parser.add_argument("--dilution-values", type=str, help="Comma separated dilution values corresponding", required=False)
    parser.add_argument("--summary-stats-file", type=str, help="Comma separated summary statistics files", required=True)
    parser.add_argument("--beta", type=float, help="Effect size", default=0.1)
    parser.add_argument("--alpha", type=float, help="Significance level", default=5*10**-8)
    parser.add_argument("--maf", type=float, help="Minor allele frequency", default=0.05)
    parser.add_argument("--tails", type=str, help="Type of tail", default="one-tailed")
    parser.add_argument("--n_samples", type=int, help="Number of samples", default=1000)
    parser.add_argument("--P", type=str, help="Comma separated proportion of cases vs controls", required=False)

    args = parser.parse_args()

    if args.phemed_summary_file and args.dilution_values:
        raise ValueError("Either phemed-summary-file or dilution-values should be provided, but not both.")

    args = parser.parse_args()

    # Save the inputs to variables
    phemed_summary_file = args.phemed_summary_file if args.phemed_summary_file else None
    dilution_values = args.dilution_values if args.dilution_values else None
    summary_stats_file = args.summary_stats_file
    #N = int(args.N)
    #output_file = args.output
    #beta = args.beta
    alpha = args.alpha
    #maf = args.maf
    P = args.P
    #min_p_val = args.min_p_val
    tails = args.tails
    n_studies = len(dilution_values.split(","))
    print(n_studies)
    n_samples = args.n_samples

    power = []
    corrected_power = []
    if phemed_summary_file != None:
        # Read the phemed summary file
        with open(phemed_summary_file, 'r') as file:
            lines = file.readlines()
    
        # Iterate over the first column to get the study
        studies = [line.split()[0] for line in lines[1:]]  # Assuming the first line is a header

        # Iterate over each study to get phi value

        for study in studies:
            phi = get_phi_value(phemed_summary_file, study)
            power += calculate_power_sum_stats(summary_stats_file, alpha, P, tails, phi,n_samples)
            #SE = 1 / (n_samples**2 * (1 - maf) * f) * math.sqrt(2 * P   * (1 - P))
            corrected_power += calculate_corrected_power_sum_stats(summary_stats_file, alpha, P, tails, phi, n_samples)

        generate_plot(power)
        generate_plot(corrected_power)
    
    #TODO : work on this elif
    elif dilution_values != None:
        for study in range(n_studies):
            phi = dilution_values.split(",")
            #Get beta, MAF from phemed summary file
            power += calculate_power_sum_stats(summary_stats_file, alpha, P, tails, phi[study],n_samples)
            #SE = 1 / (n**2 * (1 - maf) * f) * math.sqrt(2 * P   * (1 - P))
            corrected_power += calculate_corrected_power_sum_stats(summary_stats_file, alpha, P, tails, phi[study], n_samples)

        generate_plot(power)
        generate_plot(corrected_power)

    else:
       raise ValueError("Either dilution values or phemed summary file should be provided.")

    #else:
       """  #let phi be 0.5
        samples = []
        for n in simulate_sample(N):
            #print(n)
            #calculate_power(n, beta, alpha, maf, phi, P)
            phi = 1
            beta = 0.5
            alpha = 5*10**-8
            maf = 0.5
            P = 0.05
            samples.append(n) 
            power.append(calculate_power(n, beta, alpha, maf, phi, P))   
            #print(calculate_power(100, 0.5, 0.0000000000000005, 0.5, 1.5, 0.05))
            #print(calculate_power(100, 0.5, 0.0000000000000005, 0.5, 1.5, 0.05), 100)
            corrected_power.append(calculate_corrected_power(n,maf,P,beta,phi,min_p_val,alpha))

        #print(samples)
        #print(corrected_power)
        print(corrected_power == power)
        for p, cp in zip(power, corrected_power):
            if p != cp:
                print(f"Power: {p}, Corrected Power: {cp}")
        #generate_plot(power, samples, 0.5, 0.0000000000000005, 0.5, phi, 0.05, output_file)
        generate_plot(power,corrected_power, samples, beta, alpha, maf, phi, P, output_file)

    weighted_p = generate_weighted_p_thresholds(power, corrected_power)
    count_above_weighted_p = sum(cp > wp for cp, wp in zip(corrected_power, weighted_p))
    print(f"Number of corrected power values above their corresponding weighted p-values: {count_above_weighted_p}")
 """
if __name__ == '__main__':
    power_simulation()