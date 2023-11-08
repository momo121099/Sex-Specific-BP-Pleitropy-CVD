# BP-sex-specific-pleiotropy
## Description:
The purpose of this implementation is to identify BP-associated regions with sex-specific pleiotropy in relation to a cardiovascular disease (CVD) trait, using the GWAS data provided by the user. We employ sex-stratified genetic associations for systolic blood pressure (SBP), diastolic blood pressure (DBP), or pulse pressure (PP) from the UK Biobank (UKB) dataset, with equal sample sizes for different sexes, to conduct cross-trait GWAS colocalization. This analysis is specifically focused on the top BP loci regions identified in the BP GWAS.
Based on the comparison of colocalization posterior probabilities, we can pinpoint BP candidate regions with sex-specific and sex-biased pleiotropy with the CVD trait of interest (ML Yang, et al.).

## Usage:
bp.sex.colocalization(cvd.data,bp.trait='PP',cvd.trait='CVD',size=1000,p=0.3,Type='b',poster.p='H4',gene.pull.method='mean',wd=250000, diff=0.5, cutoff=0.5)

## Arguments:

bp.trait: This parameter specifies the sex-stratified GWAS data for blood pressure (BP) traits obtained from the UK Biobank (UKB), as described by ML Yang et al. You can choose from the following options: 'SBP' (systolic blood pressure), 'DBP' (diastolic blood pressure), or 'PP' (pulse pressure).

cvd.trait: This parameter allows you to input the name of the cardiovascular disease (CVD) trait of interest provided by the user. This trait name will be incorporated into the output file name for reference.

size: Specify the sample size used for the GWAS analysis of the provided CVD trait.

Type: Indicate the type of data in the CVD trait, which can be either "q" for quantitative data or "b" for case-control data.

p: For case-control CVD datasets, specify the proportion of samples in the CVD data that are cases.

poster.p: This parameter determines the posterior probability used for comparison. You can choose between "H4" or "max.H3.H4." It is recommended to use "pp.H4" when the CVD GWAS data is of European ancestry and has similar variant coverage as the UKB data in the BP study. If the CVD trait is from trans-ethnic data or has different coverage compared to our BP study data, using the maximum values of pp.H3 and pp.H4 can enhance the recall or power to detect sex-specific pleiotropic regions.

gene.pull.method: Specify the method for aggregating information from regions with the same nearest gene. You can choose between 'mean' or 'max' to either average or select the maximum values, respectively.

wd: Define the window size in base pair positions around the top index SNP for each BP region when running colocalization. The default is +/- 250Kb of the index SNP.

