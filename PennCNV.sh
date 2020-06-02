## Create intensity files from the initial Illumina report
perl PennCNV/1.0.5/split_illumina_report.pl -prefix PD_CNV/data/intensity_files/ mata_grc_mega_1_all_180305_FinalReport.txt 

## Make a population frequency of B allele (PFB) file from the LargePD cohort (this step can be performed with subpopulations depending on the PCA)
# Extract necessary columns to create SNP positions file from the initial report
awk -F '\t' '{print $1,$5,$6} NR==1779829{exit}' mata_grc_mega_1_all_180305_FinalReport.txt > snppos.txt
# Rearrange columns and change field separator from space to tab delim
awk 'BEGIN {FS=" "; OFS="\t"} {print $1,$3,$2}' snppos.txt > snppos_final.txt
# Make PFB file: (int_files_list.txt includes all samples to call the CNVs from)
perl PennCNV/1.0.5/compile_pfb.pl -listfile int_files_list.txt --snpposfile snppos_final.txt -output Large_PD.pfb


## Make GC model file to adjust for waviness
# Download from UCSC:
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gc5Base.txt.gz
sort file: sort -k 2,2 -k 3,3n  gc5Base.txt > sorted_gc5Base.txt
# Make GC model:
perl PennCNV/1.0.5/cal_gc_snp.pl sorted_gc5Base.txt Large_PD.pfb -output Large_PD.gcmodel

## CNV calling with GC model, with confidence scores for each CNV call
perl PennCNV/1.0.5/detect_cnv.pl -test -confidence -hmm PennCNV/1.0.4/lib/hhall.hmm -pfb Large_PD.pfb --list int_files_list.txt -log Large_PD_conf.log -out Large_PD_conf.adjusted.rawcnv -gcmodel Large_PD.gcmodel

## Merging adjacent CNV calls (the --fraction argument is set as 0.2 by default)
perl PennCNV/1.0.5/clean_cnv.pl combineseg Large_PD_conf.adjusted.rawcnv -signalfile Large_PD.pfb -out Large_PD_conf.clean.rawcnv
