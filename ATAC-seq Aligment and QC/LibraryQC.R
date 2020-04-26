# Analysis of BMDM ATAC seq data from Ashwood Lab
# AVC 8/9/2018 a.ciernia@gmail.com

##################################################################################
# Using esATAC package version 1.0.23:
#   This package provides a framework and complete preset pipeline for
# quantification and analysis of ATAC-seq Reads. It covers raw sequencing
# reads preprocessing (FASTQ files), reads alignment (Rbowtie2), aligned reads
# file operations (SAM, BAM, and BED files), peak calling (F-seq), genome
# annotations (Motif, GO, SNP analysis) and quality control report. The package
# is managed by dataflow graph. It is easy for user to pass variables seamlessly
# between processes and understand the workflow. Users can process FASTQ files
# through end-to-end preset pipeline which produces a pretty HTML report for
# quality control and preliminary statistical results, or customize workflow
# starting from any intermediate stages with esATAC functions easily and flexibly.
# 
# https://www.bioconductor.org/packages/release/bioc/vignettes/esATAC/inst/doc/esATAC-Introduction.html#3_customized_pipeline
# 
# If you need to use fseq, we recommend to set max memory size for java (8G, 8000M in the example). Or rJava will use the default parameter for fseq.
# 
# And ATAC QC 1.2.10
# https://bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html
#
##################################################################################
library(magrittr)
library(R.utils)
library(ChIPQC)
library(DiffBind)
library(GenomicAlignments)
library(GenomicRanges)
#library(genomation)


# bowtie2 alignment logs ---------------------------------------------------------------

#set path to folder containing alignments in subfolders
path = "/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/mmarge/bowtielogs"
setwd(path)
#get input files:
input.files <-dir(path, pattern ='*_bowtie2log.log')

#read in files and store information in a list named for each file
library(stringr)
bowtie_logs <- NULL
for(i in 1:length(input.files)){
  file <- read.table(input.files[i],header=F, sep="\t", stringsAsFactors=FALSE)
  reads <- file[,1]
  totalpereads <- stringi::stri_extract_first_words(reads[5])
  pe_concordant1 <- stringi::stri_extract_first_words(reads[8])
  Align_percent <- stringi::stri_extract_first_words(reads[19])
  
  sampleID <-word(input.files[i],1,sep = "_")
  tmp <- cbind(sampleID,totalpereads,pe_concordant1, Align_percent)
  
  bowtie_logs <- rbind(tmp, bowtie_logs)
}

bowtie_logs <- as.data.frame(bowtie_logs)
# flagstat data from deduplicated, filtered bams ---------------------------------------------------------------

#set path to folder containing alignments in subfolders
path = "/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/mmarge/alignment_output"
setwd(path)

#get input files:
input.files <-dir(path, pattern ='*_dedupFlagstat.txt')

#read in files and store information in a list named for each file
library(stringr)
depth <- NULL
for(i in 1:length(input.files)){
  file <- read.table(input.files[i],header=F, sep="\t", stringsAsFactors=FALSE)
  reads <- file[1,]
  reads <- stringi::stri_extract_first_words(reads)
  reads <- as.numeric(reads)
  sampleID <-word(input.files[i],1,sep = "_")
  tmp <- cbind(sampleID,reads)
  
  depth <- rbind(tmp, depth)
}
depth <- as.data.frame(depth)
colnames(depth) <- c("sampleID","PE_deduplreads")


# flagstat data from NFR bams (<147bp)---------------------------------------------------------------

#add in samtools flagstat info from post filtered NFR files:

flagstat <- read.csv("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/mmarge/Final_Flagstat/FilteredLibrarySizes.csv")

flagstat <- flagstat[,c(2,3)]
colnames(flagstat) <- c("sampleID","PE_NFRreads")

# peaks per sample ---------------------------------------------------------------
peaks <- read.csv("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/DiffBind/Peaks_per_sample.csv")
peaks <- peaks[,c(7,8)]
colnames(peaks) <- c("sampleID","HOMER_peaks")

#Frip in consensus peaks ---------------------------------------------------------------
Frip <- read.csv("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/DiffBind/Frip_allconsensuspeaks.csv")
Frip <- Frip[,c(7,8)]
colnames(Frip) <- c("sampleID","FRIP")

# merge dfs---------------------------------------------------------------
setwd("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/QC")

mdf2 <- merge(mdf,bowtie_logs, by="sampleID")

mdf3 <- merge(mdf2,depth, by="sampleID")

mdf4 <- merge(mdf3,flagstat, by="sampleID")

target_df<- read.xlsx("QCdataBMDM_ATACseq.xlsx",sheetIndex="Sheet1")

mdf5 <- merge(target_df,peaks, by = "sampleID", all.x=T)

merge2 <- merge(mdf5,Frip, by = "sampleID", all.x=T)

merge2$Tissue <- as.factor(merge2$Tissue)
merge2$Tissue <- relevel(merge2$Tissue, ref = "C57")
merge2$Treatment <- as.factor(merge2$Treatment)
merge2$Treatment <- relevel(merge2$Treatment, ref= "media")
merge2$sampleID <- as.character(merge2$sampleID)
#clean up
target_df <- merge2

#make numeric
target_df[,8:14] <- lapply(target_df[,8:14], function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
  

write.xlsx(target_df,"QCdataBMDM_ATACseq2.xlsx")

#############################################################################################
#RM ANOVA
#############################################################################################library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(lsmeans)
library(nlme)

variables <- colnames(target_df[c(8:12)])

outanova <- NULL
for (i in variables){
  print(i) 
  #dftmp <- mdf[,c(i,"Tissue","Treatment","batch","SampleID")]
  variable <- target_df[,i]
  tmp <- lme(variable ~ Tissue*Treatment, ~1|sampleID, data = target_df)
  anova <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
  anova$factor <- rownames(anova)
  anova <- anova[- grep("(Intercept)", anova$factor),]
  anova$comparison <- paste(i)
  outanova <- rbind(outanova,anova)
}

anova1 <- outanova

variables2 <- colnames(target_df[c(13:14)])
target_df2 <- target_df[complete.cases(target_df),]
outanova <- NULL
for (i in variables2){
  print(i) 
  #dftmp <- mdf[,c(i,"Tissue","Treatment","batch","SampleID")]
  variable <- target_df2[,i]
  tmp <- lme(variable ~ Tissue*Treatment, ~1|sampleID, data = target_df2)
  anova <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
  anova$factor <- rownames(anova)
  anova <- anova[- grep("(Intercept)", anova$factor),]
  anova$comparison <- paste(i)
  outanova <- rbind(outanova,anova)
}

anova2 <- outanova

anovadf <- rbind(anova1,anova2)

write.csv(anovadf,"RMANOVA_ATACqc.csv")

#Fragment distributions ---------------------------------------------------------------
#set path to folder containing alignments in subfolders
path = "/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/mmarge/alignment_output"
setwd(path)

#get input files:
input.files <-dir(path, pattern ='*_insertmetric.tsv')

#read in files and store information in a list named for each file
library(stringr)
df <- as.data.frame(seq(1:1000))
colnames(df) <- c("insert_size")

for(i in 1:length(input.files)){
  file <- read.table(input.files[i],skip=10, stringsAsFactors=FALSE, header=T)
  colnames(file)[2] <- word(input.files[i],1,sep = "_")
  
  df <- merge(df,file, all.x=T)

}

#replace NA with 0
df <- as.data.frame(df)
df[is.na(df)] <- 0

#read in and match column names
samples <- read.csv("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/libQCreport/sampleinfoBMDM_ATACseq2.csv")
samples$Factor <- as.character(samples$Factor)
#samples$SampleID <- as.character(samples$SampleID)
#rename
#colnames(df)[2:ncol(df)] <- plyr::mapvalues(colnames(df)[2:ncol(df)], from = samples$Factor, to = samples$SampleID)

library(dplyr)
df2 <- df %>% group_by(insert_size) %>%
  gather(sampleID,count,2:ncol(df))

df3 <- merge(df2, samples, by.x = "sampleID", by.y = "Factor") 

df3$Condition <- factor(df3$Condition, levels <-c("C57 media",  "C57 LPS1",   "C57 LPS2",   "BTBR media", "BTBR LPS1",  "BTBR LPS2" ))
df3$Replicate <- factor(df3$Replicate)
df4 <- df3 %>% arrange(Condition)

setwd("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/libQCreport")
#plot
ggplot(df4, aes(insert_size,count,color = Replicate)) +
  facet_wrap(~Condition) +
  geom_path() +
  xlab("Insert Size (bp)")+
  ylab("Read counts") + 
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  geom_vline(xintercept = c(147,294), linetype=2) + 
  labs(title = "Insert Size Distribution") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(0, 500))
  #scale_y_continuous(limits = c(0, 500))+
  #theme(legend.position="none")

ggsave(filename = "InsertSizeDistributions.pdf", plot = last_plot(), width = 11, height = 8.5, units = c("in"), dpi = 300)



#normalize read counts by library size

lib.size <- colSums(df[,2:ncol(df)])

scaled_DF <- mapply(`/`, df[,2:ncol(df)], lib.size)
scaled_DF <- as.data.frame(scaled_DF)
scaled_DF$insert_size <- df$insert_size

df2 <- scaled_DF %>% group_by(insert_size) %>%
  gather(sampleID,count,1:(ncol(df)-1))

df3 <- merge(df2, samples, by.x = "sampleID", by.y = "Factor") 

df3$Condition <- factor(df3$Condition, levels <-c("C57 media",  "C57 LPS1",   "C57 LPS2",   "BTBR media", "BTBR LPS1",  "BTBR LPS2" ))
df3$Replicate <- factor(df3$Replicate)
df4 <- df3 %>% arrange(Condition)

setwd("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/libQCreport")
#plot
ggplot(df4, aes(insert_size,count,color = Replicate)) +
  facet_wrap(~Condition) +
  geom_path() +
  xlab("Insert Size (bp)")+
  ylab("Read Counts/Library Size") + 
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  geom_vline(xintercept = c(147,294), linetype=2) + 
  labs(title = "Insert Size Distribution") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(0, 500))
#scale_y_continuous(limits = c(0, 500))+
#theme(legend.position="none")

ggsave(filename = "InsertSizeDistributions_LibraryScaled.pdf", plot = last_plot(), width = 11, height = 8.5, units = c("in"), dpi = 300)

