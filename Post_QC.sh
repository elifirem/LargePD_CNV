## Set log R ratio SD (LRRSD) threshold so stringent that no sample passes (to get summary stats): 
perl PennCNV/1.0.5/filter_cnv.pl Large_PD_conf.clean.rawcnv -qclogfile Large_PD_conf.log -qclrrsd 0.0001 -qcpassout Large_PD_conf.clean.qcpass -qcsumout Large_PD_conf.clean.qcsum -out Large_PD_conf.clean.goodcnv

## In R:
# Read in the .qcsum file
# Plot the distribution of number of CNVs per sample in a histogram and decide the cutoff
# Exclude all samples which have more CNVs than the cutoff

# By using the samples that passed the number of CNV cutoff:
# Calculate the median+3*SD of the LRRSD for qclrrsd
# Calculate the median+3*SD of the B allele frequency (BAF) drift for qcbafdrift
# Calculate the median+3*SD of the Waviness for qcwf


## Implement the calculated values above in this command
perl PennCNV/1.0.5/filter_cnv.pl Large_PD_conf.clean.rawcnv -qclogfile Large_PD_conf.log -qcnumcnv 1639 -qclrrsd 0.2720423 -qcwf 0.03110594 -qcbafdrift 0.001405433 -qcpassout Large_PD_conf2.clean.qcpass -qcsumout Large_PD_conf2.clean.qcsum -out Large_PD_conf2.clean.goodcnv

# Filter out spurious regions: 
# For telomeric regions, one can treat the 100kb or 500kb region within start or end of chromosome as telomeric region.
# For centromeric regions: provided by pennCNV add 500000 to each end of region
# For immunogloubin regions liftover to transfer coordinated from hg18 to hg19 (hg18 coordinates are provided by PennCNV)
perl PennCNV/1.0.5/scan_region.pl Large_PD_conf2.clean.goodcnv spurious_regions.txt -minqueryfrac 0.5 > Large_PD_conf2.spurious; fgrep -v -f Large_PD_conf2.spurious Large_PD_conf2.clean.goodcnv > Large_PD_conf2.clean.goodcnv_regionfiltered

