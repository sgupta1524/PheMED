# __FAQ for PheMED__
Phenotypic Measurement of Effective Dilution

### Q1: Why did I get a effective dilution value smaller than 1?
This means that your reference GWAS, corresponding to the log odds ratios that appears first in your input does not have the highest deltaP value = NPV + PPV - 1.  The analysis is correct, but if you want all of your phenotypic dilution values to be greater (or equal to) 1, this issue can be remedied by changing your reference GWAS to the GWAS that achieves the smallest effective dilution value.

### Q2: Why don't I get a confidence interval for the reference study?
The (relative) effective dilution for the reference study is defined to equal 1.  Hence, a confidence interval is not meaningful.  In general, it is easier to interpret the results when the reference study refers to the study with the lowest effective dilution.  In such a case, there may be uncertainty regarding which study has the lowest effective dilution and you may see the confidence intervals for the relative effective dilution of other studies capture values lower than 1.  (See Q1 for further details.)

### Q3: I got a warning that the number of valid SNPs is low.  How do I correct this?
PheMED performs a merge on the GWAS summary stats to identify approximately independent SNPs.  If the number of retained SNPs is very low, this will output a warning.  You can try generating your own list of approximately independent SNPs based on the SNPs in your summary statistics to see if you can remedy the issue (see Q4).  Alternatively, this may mean that the number of (approximately independent) SNPs in your data set is very small.  In such a case, PheMED will produce very wide CIs for the effective dilution.  

### Q4: How can I generate my own list of randomly generated approximately independent SNPs?
For a given list of SNPs, we generate a simulated p-value for each SNP using independent uniform random variables.  We can then [clump the SNPs](https://zzz.bwh.harvard.edu/plink/clump.shtml) using the simulated p-values with a reference panel.  For our data, we used the default values listed in the Plink documentation.  By retaining only the lead SNPs from each clump, our SNPs are approximately independent.  We can repeat this process for multiple ancestries, by applying the clumping procedure iteratively on the retained SNPs for multiple reference panels of different ancestries.  

### Q5: How are missing values handled?
If there are missing values (e.g. effect sizes) for some of the studies in the data set, PheMED will impute the effect sizes to be zero with a large value for the standard deviation (e.g. 1000), such that the (meta-analyzed) output will be approximate to the case, where the missing values would be dropped completely from the analysis.

### Q6: How do I run PheMED if I have sample overlap across studies?
PheMED requires at least one study not overlap with the other studies in the analysis.  In this case, users should run multiple PheMED analyses, such that there is no overlap in any individual analysis.  For example if you want to analyze PheMED on three studies, where Studies 2 and 3 have sample overlap and Study 1 does not have sample overlap, users should run PheMED twice, where one run measures the dilution between Studies 1 and 2, whereas another run measures the dilution between Studies 1 and 3.  See our paper for details.  A future update will help streamline analyses, where studies have sample overlap.  
