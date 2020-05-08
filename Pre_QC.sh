##  Call rate filtering: filters out all samples with missing call rates exceeding 0.04 to be removed
plink --bfile LARGE_PD_final --mind 0.04 --make-bed --keep-allele-order --out LARGE_PD_final_CR_filtered


##  Merge Large PD samples with 1000 Genomes (“A global reference for human genetic variation” Nature 526, 68–74, 2015) as references for PCA
# Rename the SNP IDs --> chr_position_Allele2 (Allele2 is usually the major one) e.g. 1_58814_G
awk '{print $1"\t" $1"_" $4"_" $6"\t" $3"\t" $4"\t" $5"\t" $6 }' LARGE_PD_final_CR_filtered.bim > LARGE_PD_final_CR_filtered.bim
# Make same SNP pool: Only use SNPs in both datasets (1000 Genomes and Large PD)
cut -f2 1kG.bim > same_pool_list
plink --bfile LARGE_PD_final_CR_filtered --extract same_pool_list --keep-allele-order --make-bed --out LARGE_PD_final_CR_filtered_same
# Remove multiallelic variants
plink --bfile LARGE_PD_final_CR_filtered_same --exclude For_PCA-merge.missnp --make-bed --keep-allele-order –out LARGE_PD_final_CR_filtered_same_nonMA
# Merge 1000 Genomes and Large PD samples
plink --bfile LARGE_PD_final_CR_filtered_same_nonMA --bmerge 1kG.bed 1kG.bim 1kG.fam --make-bed --keep-allele-order --allow-no-sex --out For_PCA


##  Filtering
# MAF >5%
plink --bfile For_PCA --maf 0.05 --make-bed --keep-allele-order --out For_PCA_maf5
# SNP missinigness <2%
plink --bfile For_PCA_maf5 --geno 0.02 --make-bed --keep-allele-order --out For_PCA_maf5_geno2
# Missinigness per indiv <5%
plink --bfile For_PCA_maf5_geno2 --mind 0.05 --make-bed --keep-allele-order --out For_PCA_maf5_geno2_mind5
# HWE p-value >1e-3
plink --bfile For_PCA_maf5_geno2_mind5 --hwe 0.001 --make-bed --keep-allele-order --out For_PCA_maf5_geno2_mind5_hwe
# Pairwise IBD estimation
plink --bfile For_PCA_maf5_geno2_mind5_hwe --Z-genome --keep-allele-order --out For_PCA_maf5_geno2_mind5_hwe_IBD
# List of all samples with PI_HAT > 0.2
zless For_PCA_maf5_geno2_mind5_hwe_IBD.genome.gz | awk '$10 > 0.2' > relatives.txt
# Exclude SNPs in known complicated regions on the genome eg. MHC (chr6, 25-35Mb) and chr8 inversion region (chr8,7-13Mb)
plink --bfile For_PCA_maf5_geno2_mind5_hwe --chr 6 --from-mb 25 --to-mb 35 --make-bed --keep-allele-order --out MHC
plink --bfile For_PCA_maf5_geno2_mind5_hwe --chr 8 --from-mb 7 --to-mb 13 --make-bed --keep-allele-order --out Inver
awk '{print $2}' MHC.bim  > MHC_snps
awk '{print $2}' Inver.bim  > Inver_snps
cat MHC_snps Inver_snps > remove_snps.txt
plink --bfile For_PCA_maf5_geno2_mind5_hwe --exclude remove_snps.txt --make-bed --keep-allele-order --out For_PCA_maf5_geno2_mind5_hwe_MHCInverfree

##  LD pruning
# Randomly pick 300 individuals
awk '{print $1 "\t" $2}' For_PCA_maf5_geno2_mind5_hwe_MHCInverfree.fam > all_indiv
shuf -n 300 all_indiv > random300
# Pruning
plink --bfile For_PCA_maf5_geno2_mind5_hwe_MHCInverfree --keep random300 --make-bed --keep-allele-order --out Random300
plink --bfile Random300 --indep-pairwise 200 100 0.2 --keep-allele-order --out Random300_LDpruned
plink --bfile For_PCA_maf5_geno2_mind5_hwe_MHCInverfree --exclude Random300_LDpruned.prune.out --make-bed --keep-allele-order --out For_PCA_maf5_geno2_mind5_hwe_MHCInverfree_LDpruned

##  PCA
plink --bfile For_PCA_maf5_geno2_mind5_hwe_MHCInverfree_LDpruned --pca header tabs --keep-allele-order --out data_pca
