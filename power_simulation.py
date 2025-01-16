import scipy.stats as stats
import math
import numpy as np
from statsmodels.stats.multitest import multipletests

def parse_sum_stats(inputs, effect_allele_cols, non_effect_allele_cols):
    """
    Read input summary stats files in chunks, return beta and MAF
    """

def simulate_data(num_simulations):
    """
    Simulate data for the power function.
    num_simulations = number of simulations to run
    """
    np.random.seed(42)  # For reproducibility

    # Generate random values for the parameters

    ##plausible values

    #from 1000 to 1 million
    n_values = np.random.randint(, num_simulations)  # Sample size between 100 and 10000
    #between -.0.1 and 0.1
    beta_values = np.random.normal(0, 0.01, num_simulations)  # Effect size from a normal distribution with mean 0 and std 1
    #0.05 or 5*10e**-8
    alpha_values = np.random.uniform(0.01, 0.1, num_simulations)  # Significance level between 0.01 and 0.1
    maf_values = np.random.uniform(0.01, 0.5, num_simulations)  # Minor allele frequency between 0.01 and 0.5
    phi_values = np.random.uniform(0.5, 1.5, num_simulations)  # Dilution factor between 0.5 and 1.5
    P = 0.05  # 5% typically for neuropsychiatric traits

    # Calculate power for each set of parameters
    power_values = []
    for i in range(num_simulations):
        power_value = power(n_values[i], beta_values[i], alpha_values[i], maf_values[i], phi_values[i], P_values[i])
        power_values.append(power_value)

    return power_values

def power(n, beta, alpha, maf, phi, P):
    """
    Calculate the power of a GWAS study 
    n = sample size
    beta = effect size
    alpha = significance level
    maf = frequency of the risk allele
    phi = dilution factor
    P = proportion of cases
    """
    # Calculate Standard error based on a logistic model
    SE = 1 / (n**2 * (1 - maf) * f) * math.sqrt(2 * P   * (1 - P))
    # Calculate the non-centrality parameter (NCP) for the chi-square distribution
    NCP = (beta / phi*SE) ** 2
    # Calculate the power of the study
    power = 1 - stats.chi2.cdf(stats.chi2.ppf(alpha, df=1), df=1, nc=NCP)
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

def FIQT(beta, SE, min_p):

    # Calculate z-scores
    z = beta/SE
    # Calculate p-values
    pvals = 2 * norm.sf(np.abs(z))
    
    # Replace p-values less than min_p with min_p
    pvals[pvals < min_p] = min_p
    
    # Adjust p-values using the FDR method
    _, adj_pvals, _, _ = multipletests(pvals, method='fdr_bh')
    
    # Calculate adjusted z-scores
    mu_z = np.sign(z) * norm.isf(adj_pvals / 2)
    
    # Replace adjusted z-scores for very large z-scores
    threshold = norm.isf(min_p / 2)
    mu_z[np.abs(z) > threshold] = z[np.abs(z) > threshold]
    
    return mu_z

def power_corrected(beta,SE,phi,min_p_val,alpha):
    """
    Calculate the power of a GWAS study 
    n = sample size
    beta = effect size
    alpha = significance level
    maf = frequency of the risk allele
    phi = dilution factor
    P = proportion of cases
    """
    z = beta/SE
    corrected_z = FIQT(z, min_p_val)
    NCP = (corrected_z/phi) ** 2
    # Calculate the power of the study
    power = 1 - stats.chi2.cdf(stats.chi2.ppf(alpha, df=1), df=1, nc=NCP)
    # Return the power value per SNP
    return power