library(stringr)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)

cnv_data    <- as.matrix(read.table("data/Large_PD_conf2.clean.goodcnv_regionfiltered", sep="\t", header=F, stringsAsFactors = F))
sample_data <- read.table("data/Large_PD_conf2.clean.qcsum", sep = '\t', header=T, stringsAsFactors = F)

cnv_data   <- t(apply(cnv_data, 1, FUN =  function(CNV.line){
  #browser()
  CNV.lines.tp   <- CNV.line
  CNV.line       <- strsplit(CNV.lines.tp, ' ')[[1]]
  snp.idx        <- grep('numsnp', CNV.line)
  cn.idx         <- grep('state', CNV.line)
  id.idx         <- grep('/', CNV.line)
  conf.idx       <- grep('conf', CNV.line)
  tmp            <- str_extract_all(CNV.line[1], "[0-9]+")
  tmp2           <- str_extract_all(CNV.line[cn.idx], "[0-9]+")
  tmp3           <- strsplit(CNV.line[id.idx ], "/")
  tmp4           <- str_extract_all(CNV.line[conf.idx], "[0-9]+")
  tmp5           <- str_extract_all(CNV.line[snp.idx], "[0-9]+")
  tmp6           <- str_extract_all(CNV.line[3], "[0-9]+")
  chr            <- as.integer(tmp[[1]][1])
  start          <- as.integer(tmp[[1]][2])
  end            <- as.integer(tmp[[1]][3])
  if(length(cn.idx) > 0){cn <- as.integer(tmp2[[1]][2])}else{cn  <- NA}
  if(length(conf.idx) > 0){conf <- as.integer(tmp4[[1]][1]) + as.integer(tmp4[[1]][2]) / 1000}else{conf  <- 0; print('No conf score')}
  if(!is.na(as.numeric(CNV.line[id.idx+1]))){
    id      <- paste(tmp3[[1]][length(tmp3[[1]])], CNV.line[id.idx+1], sep = '_')
  }else{
    id      <- tmp3[[1]][length(tmp3[[1]])]
  }
  if(length(snp.idx) > 0){numSNP  <- tmp5[[1]][1]}else{numSNP  <- NA}
  
  l.val   <- end - start
  return(c('Chromosome' = chr, 'Start_Position_bp' = start, 'End_Position_bp' = end, 'Copy_Number' = cn, 'Max_Log_BF' = conf, 'Sample_Name' = id, 'No_Probes' = numSNP, 'Length_bp' = l.val))
}))


cnv_data    <- cnv_data[,c('Sample_Name', 'Chromosome', 'Start_Position_bp', 'End_Position_bp', 'Copy_Number', 'Max_Log_BF', 'Length_bp', 'No_Probes')]
cnv_data    <- as.data.frame(cnv_data, stringAsFactors = FALSE)

cnv_data$Sample_Name        <- as.character(cnv_data$Sample_Name)
cnv_data$Chromosome         <- as.numeric(as.character(cnv_data$Chromosome))
cnv_data$Start_Position_bp  <- as.numeric(as.character(cnv_data$Start_Position_bp))
cnv_data$End_Position_bp    <- as.numeric(as.character(cnv_data$End_Position_bp))
cnv_data$Copy_Number        <- as.numeric(as.character(cnv_data$Copy_Number))
cnv_data$Max_Log_BF         <- as.numeric(as.character(cnv_data$Max_Log_BF))
cnv_data$Length_bp          <- as.numeric(as.character(cnv_data$Length_bp))
cnv_data$No_Probes          <- as.numeric(as.character(cnv_data$No_Probes))


# remove CNV with NA confidence score
if(length(which(is.na(cnv_data$Max_Log_BF))) > 0){
  cnv_data<- cnv_data[-which(is.na(cnv_data$Max_Log_BF)),]
}

#  Clean the Sample ID for the sample data
sample_data$File    <- sapply(as.character(sample_data$File), FUN = function(data){
  val <- strsplit(data, '/')
  return(val[[1]][length(val[[1]])])
})

#  Merge CNV and sample data
cnv_data_global  <- merge(cnv_data, sample_data, by.x = 'Sample_Name', by.y = 'File')

#  Count number of samples and CNV before cleaning
log_data <- c(nb_sample_before = length(unique(cnv_data_global$Sample_Name)), nb_cnv_before = nrow(cnv_data_global))
log_data 


###########################################################################################
#-CNV filters, Post QC Step 2
#-Create density columns
cnv_data_global <- mutate(cnv_data_global, density = No_Probes / Length_bp)

#-Passing filter
largePD_pass <- subset(cnv_data_global, No_Probes >=20 & Length_bp>=20000 & Length_bp<1000000 & density>=0.0001)
largePD_pass$Flag <- "In"
largePD_1mb <- subset(cnv_data_global, Length_bp>=1000000 & No_Probes >=20)
largePD_1mb$Flag <- "In"

#-Not passing filter
largePD_fail1 <- subset(cnv_data_global, !(No_Probes >=20 & Length_bp>=20000 & Length_bp<1000000 & density>=0.0001))
largePD_fail2 <- subset(largePD_fail1, !(Length_bp>=1000000 & No_Probes >=20))
largePD_fail2$Flag <- "Out"          

#-Combine flagged CNVs
largePD_all <- rbind(largePD_pass,largePD_1mb,largePD_fail2)

#-Count CNVs per Sample only using "In" flag
largePD_in <- subset(largePD_all, Flag == "In")
largePD_NumCNVs <- as.data.frame(table(largePD_in$Sample_Name))
largePD_out <- subset(largePD_all, Flag == "Out")
largePD_fail_NumCNVs<-as.data.frame(table(largePD_out$Sample_Name))

colnames(largePD_fail_NumCNVs)<-c("Sample","NumCNV")
colnames(largePD_NumCNVs) <- c("Sample", "NumCNV")

#-Add a second flag by using mean+3*sd
summary(largePD_NumCNVs$NumCNV)
hist(largePD_NumCNVs$NumCNV, breaks = 100)
ggplot(largePD_NumCNVs, aes(largePD_NumCNVs$NumCNV)) + 
  geom_histogram(binwidth = 0.05)+
  xlab("\nLog(#CNV per sample)") +
  ylab("Counts\n")+
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed")+
  geom_vline(xintercept = 3, colour = "red", linetype = "dashed")+
  geom_vline(xintercept = 5, colour = "red", linetype = "dashed")+
  geom_vline(xintercept = 10, colour = "red", linetype = "dashed")+
  geom_vline(xintercept = 13, colour = "red", linetype = "dashed")+
  geom_vline(xintercept = 250, colour = "red", linetype = "dashed")+
  scale_x_log10()

#median(largePD_NumCNVs$NumCNV)+3*sd(largePD_NumCNVs$NumCNV)
mean(largePD_NumCNVs$NumCNV)+3*sd(largePD_NumCNVs$NumCNV)

largePD_NumCNVs$Flag2 <- ifelse(largePD_NumCNVs$NumCNV > 96.8, "Out8.2", "In")
colnames(largePD_NumCNVs)[1] <- "Sample_Name"
largePD_in$Sample_Name <- as.factor(largePD_in$Sample_Name)
largePD_in_new <- left_join(largePD_in, largePD_NumCNVs, "Sample_Name")

#-To merge excluded ones
largePD_fail_NumCNVs$Flag2 <- "Out8.1"
colnames(largePD_fail_NumCNVs)[1] <- "Sample_Name"
largePD_out$Sample_Name <- as.factor(largePD_out$Sample_Name)
largePD_out_new <- left_join(largePD_out, largePD_fail_NumCNVs, "Sample_Name") 

cnv_data_global_all <- rbind(largePD_in_new, largePD_out_new) # All CNVs with flags
cnv_data_global_ins <- subset(cnv_data_global_all, Flag2 == "In") # To count passed CNVs so far

#  Count number of samples and CNV after cleaning
log_data    <- c(log_data, nb_sample_after = length(unique(droplevels(cnv_data_global_ins$Sample_Name))), 
                 nb_cnv_after = nrow(cnv_data_global_ins), nb_deletion = length(which(cnv_data_global_ins$Copy_Number < 2)), 
                 nb_duplication = length(which(cnv_data_global_ins$Copy_Number > 2)))
log_data
log_data    <- c(log_data, nb_sample_after = length(unique(droplevels(cnv_data_global_all$Sample_Name))), 
                 nb_cnv_after = nrow(cnv_data_global_all), nb_deletion = length(which(cnv_data_global_all$Copy_Number < 2)), 
                 nb_duplication = length(which(cnv_data_global_all$Copy_Number > 2)))
log_data


####################################################################################################

#  Save the data - QC passed version
save(cnv_data_global_ins, file = 'data/cnv_data_global_ins.rdata')
save(cnv_data_global_all, file = 'data/cnv_data_global_all.rdata')


#################################### QS calculation ################################################
# Define function
calculate_QS <- function(cnv_data, model_data){
  QS_val <- rep(model_data['(Intercept)', 'Estimate'], nrow(cnv_data))
  for(para in 2:nrow(model_data)){
    para_name   <- rownames(model_data)[para]
    QS_val      <- QS_val + model_data[para_name, 'Estimate'] * scale(cnv_data[, para_name])
  }
  QS_val  <- 1 / (1 + exp(-QS_val))
  return(QS_val)
}


#  Create separate variable for the deletions and duplications
cnv_data_global_del <- subset(cnv_data_global_all, subset = Copy_Number < 2)
cnv_data_global_dup <- subset(cnv_data_global_all, subset = Copy_Number > 2)

#  define the QS coefficients
QS_coefficient_duplication              <- cbind(c(-0.81529185, 3.90248885, -0.90184846, 0.34376765, -0.06351849, -0.01497046, -0.30665878, -0.10354156, -0.44829901, 0.06039380, 0.03645913), rep(0, 11))
rownames(QS_coefficient_duplication)    <- c('(Intercept)', 'Max_Log_BF', 'No_Probes', 'LRR_SD', 'BAF_drift', 'WF', 'LRR_mean', 'NumCNV', 'BAF_SD', 'Length_bp', 'BAF_mean')
colnames(QS_coefficient_duplication)    <- c('Estimate', 'Pr(>|z|)')

QS_coefficient_deletion             <- cbind(c(-1.75648377, 4.64898485, -2.50150285, -0.47552224, 0.28876320, 0.10205302, 0.14363692, 0.02959571, 0.00000000, -0.00000000), rep(0, 10))
rownames(QS_coefficient_deletion)   <- c('(Intercept)', 'Max_Log_BF', 'No_Probes', 'NumCNV', 'LRR_SD', 'Length_bp', 'LRR_mean', 'WF', 'BAF_drift', 'BAF_SD')
colnames(QS_coefficient_deletion)   <- c('Estimate', 'Pr(>|z|)')


#  Calculate the Quality Score
QS_deletion     <- calculate_QS(cnv_data = cnv_data_global_del, model_data = QS_coefficient_deletion)
QS_duplication  <- calculate_QS(cnv_data = cnv_data_global_dup, model_data = QS_coefficient_duplication)

QS                                          <- matrix(NA, nrow = nrow(cnv_data_global_all), ncol = 1)
QS[which(cnv_data_global_all$Copy_Number < 2),] <- QS_deletion * -1
QS[which(cnv_data_global_all$Copy_Number > 2),] <- QS_duplication
colnames(QS)                                <- c('Quality_Score')
cnv_data_global_QS    <- cbind(cnv_data_global_all, QS)


#  Save the new CNV dataset with the Quality Score for deletions and duplications
save(cnv_data_global_QS, file = 'data/cnv_data_global_QS.rdata') 

