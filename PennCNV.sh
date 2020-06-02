## Create intensity files from the initial Illumina report
perl PennCNV/1.0.4/split_illumina_report.pl -prefix PD_CNV/data/intensity_files/ PD_CNV/data/mata_grc_mega_1_all_180305_FinalReport.txt 

## Make a population frequency of B allele (PFB) file from the LargePD cohort (this step can be performed with subpopulations depending on the PCA)
# Extract necessary columns to create SNP positions file from the initial report
awk -F '\t' '{print $1,$5,$6} NR==1779829{exit}' mata_grc_mega_1_all_180305_FinalReport.txt > snppos.txt
# Rearrange columns and change field separator from space to tab delim
awk 'BEGIN {FS=" "; OFS="\t"} {print $1,$3,$2}' snppos.txt > snppos_final.txt
# Make PFB file: (int_files_list.txt includes all samples to call the CNVs from)
perl PennCNV/1.0.4/compile_pfb.pl -listfile PD_CNV/data/int_files_list.txt --snpposfile snppos_final.txt -output Large_PD.pfb




## Make GC model file to adjust for waviness
# Download from UCSC:
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gc5Base.txt.gz
sort file: sort -k 2,2 -k 3,3n  gc5Base.txt > sorted_gc5Base.txt
# Make GC model:
PennCNV-1.0.4/cal_gc_snp.pl sorted_gc5Base.txt EUR_new.pfb -output EUR.gcmodel

## CNV calling with GC model 
PennCNV-1.0.4/detect_cnv.pl -test -hmm PennCNV-1.0.4/lib/hhall.hmm -pfb EUR_new.pfb --list all_intens_files.txt -log EUR_CNVs.log -out EUR_CNVs.adjusted.rawcnv  -gcmodel EUR.gcmodel

## Merging adjacent CNV calls (the --fraction argument is set as 0.2 by default)
PennCNV-1.0.4/clean_cnv.pl combineseg EUR_CNVs.adjusted.rawcnv -signalfile EUR_new.pfb -out EUR_CNVs_clean.rawcnv
