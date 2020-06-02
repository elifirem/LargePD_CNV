for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv For_PCA_maf5_geno2_mind5_hwe_MHCInverfree_LDpruned.bed $K -j10 | tee log${K}.out; done
