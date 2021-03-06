# LARGE-PD CNV

This repository contains the scripts for quality control (QC) steps and CNV calling from the study: [Genome-wide Analysis of Copy Number Variation in Latin American Parkinson’s Disease Patients](https://www.medrxiv.org/content/10.1101/2020.05.29.20100859v2) using [LARGE-PD](https://large-pd.org/) cohort


**Usage**:

`Pre_QC.sh` plink command lines for the QC steps before PCA, includes admixture run, and relatedness estimates  

`PennCNV.sh` CNV calling with PennCNV

`Post_QC.sh` command lines for post QC filtering steps

`Code_QS.R` for further filtering CNV calls: calculates a quality score (QS) which estimates the probability of a CNV called by PennCNV to be confirmed by other software, adapted from [Macé et al.](https://pubmed.ncbi.nlm.nih.gov/27402902/)


**Software**:

1. [PLINK 1.9:](https://www.cog-genomics.org/plink2/) Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.

2. [ADMIXTURE:](http://software.genetics.ucla.edu/admixture/) Alexander DH, Novembre J, Lange K. (2009) Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655–1664, 2009.

3. [KING:](http://people.virginia.edu/~wc9c/KING/) Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) Robust relationship inference in genome-wide association studies. Bioinformatics 26(22):2867-2873.

4. [PennCNV:](http://penncnv.openbioinformatics.org/en/latest/) Wang K, Li M, Hadley D, Liu R, Glessner J, Grant S, Hakonarson H, Bucan M (2007) An integrated hidden Markov model designed for high-resolution copy number variation detection in whole-genome SNP genotyping data. Genome Res. 17, 1665–1674.

5. [R:](https://www.R-project.org/) R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

6. [Perl:](https://www.perl.org/get.html) The Perl programming language.

