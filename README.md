# BP-sex-specific-pleiotropy
## Description:
The purpose of this implementation is to identify BP-associated regions with sex-specific pleiotropy in relation to a cardiovascular disease (CVD) trait, using the GWAS data provided by the user. We employ sex-stratified genetic associations for systolic blood pressure (SBP), diastolic blood pressure (DBP), or pulse pressure (PP) from the UK Biobank (UKB) dataset, with equal sample sizes for different sexes, to conduct cross-trait GWAS colocalization. This analysis is specifically focused on the top BP loci regions identified in the BP GWAS.
Based on the comparison of colocalization posterior probabilities, we can pinpoint BP candidate regions with sex-specific and sex-biased pleiotropy with the CVD trait of interest (ML Yang, et al.).

## Usage:
bp.sex.colocalization(cvd.data,bp.trait='PP',cvd.trait='CVD',size=1000,p=0.3,Type='b',poster.p='H4',gene.pull.method='mean',wd=250000, diff=0.5, cutoff=0.5)

## Arguments:
*cvd.data: This parameter represents the user's GWAS data for the cardiovascular disease (CVD) trait of interest. The format for this input file is described below.

*bp.trait: Specify the sex-stratified GWAS data for blood pressure (BP) traits sourced from the UK Biobank (UKB), as detailed by ML Yang et al. You have the option to choose from 'SBP' (systolic blood pressure), 'DBP' (diastolic blood pressure), or 'PP' (pulse pressure).

*cvd.trait: Input the name of the user's cardiovascular disease (CVD) trait of interest. This trait name will be incorporated into the output file name for reference.

*size: Indicate the sample size used for the GWAS analysis of the provided CVD trait.

*Type: Specify the type of data in the CVD trait, which can be either "q" for quantitative data or "b" for case-control data.

*p: For case-control CVD datasets, define the proportion of samples in the CVD data that are cases.

*poster.p: This parameter determines the colocalization posterior probability used for comparison. You can choose between "H4" or "max.H3.H4." It is advisable to use "pp.H4" when the CVD GWAS data is of European ancestry and shares similar variant coverage with the UKB data used in the BP study. If the CVD trait originates from trans-ethnic data or has substantially different coverage compared to our BP study data, utilizing the maximum values of pp.H3 and pp.H4 can enhance the ability to detect sex-specific pleiotropic regions.

*gene.pull.method: Specify the method for aggregating colocalization posterior probability information from regions sharing the same nearest gene. You can opt for either 'mean' or 'max' to calculate either the average or maximum values, respectively.

*wd: Define the window size in base pair positions around the top index SNP for each BP region during colocalization. The default is +/- 250Kb of the index SNP.

*diff: Set the threshold for the absolute difference in posterior probabilities between male and female BP colocalization with CVD. This threshold is employed to identify genes that exhibit sex-biased colocalization.

*cutoff: Establish a threshold for the posterior probability to identify genes displaying sex-specific colocalization. Genes meeting this criterion will have a posterior probability greater than this cutoff in only one sex. The combination of the diff and cutoff values will be used together to screen for the top BP genes displaying sex-specific pleiotropic effects with CVD.
