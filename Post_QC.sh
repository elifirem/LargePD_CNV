## Set lrrsd threshold so stringent that no sample passes(to get summary stats): 
perl PennCNV/1.0.5/filter_cnv.pl Large_PD_conf.clean.rawcnv -qclogfile Large_PD_conf.log -qclrrsd 0.0001 -qcpassout Large_PD_conf.clean.qcpass -qcsumout Large_PD_conf.clean.qcsum -out Large_PD_conf.clean.goodcnv

## In R:
# Read in the .qcsum file
# Plot the distribution of number of CNVs per sample in a histogram and decide the cutoff


perl PennCNV/1.0.5/filter_cnv.pl Large_PD_conf.clean.rawcnv -qclogfile Large_PD_conf.log -qcnumcnv 1639 -qclrrsd 0.2720423 -qcwf 0.03110594 -qcbafdrift 0.001405433 -qcpassout Large_PD_conf2.clean.qcpass -qcsumout Large_PD_conf2.clean.qcsum -out Large_PD_conf2.clean.goodcnv
perl PennCNV/1.0.5/scan_region.pl Large_PD_conf2.clean.goodcnv spurious_regions.txt -minqueryfrac 0.5 > Large_PD_conf2.spurious; fgrep -v -f Large_PD_conf2.spurious Large_PD_conf2.clean.goodcnv > /home/sarihae/isilon/Irem/PD_CNV/data/

# Move the output Large_PD_conf2.clean.goodcnv_regionfiltered to R and 
#follow 1_Code_QS.R to produce the file Large_PD_QS2.QC2
