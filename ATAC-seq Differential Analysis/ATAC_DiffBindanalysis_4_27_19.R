##########################################################################
#AVC a.ciernia@gmail.com
#8/2018
#ATACseq BMDM +/- 1 or 2 doses LPS in vitro
##########################################################################
#Analysis of Differential Accessibility in ATAC seq Data
#Peaks called on Nucleosome Free Reads (reads < 147bp) that were filtered
#to remove mitochonrida reads, PCR duplicates, low quality, and for only properly paired reads
#Uses DiffBind to build a DBA object and and identify consensus peaks amoung replicates (3/4)
#Statistically differentially accessible regions are then found from consensus peaks 
#using EdgeR+limma-voom with factors for treatment x genotype and including a corrleation weight matrix to account for the repeated measures design
##########################################################################
#DiffBind analysis package: https://rdrr.io/bioc/DiffBind/f/inst/doc/DiffBind.pdf
## try http:// if https:// URLs are not supported
#source("http://bioconductor.org/biocLite.R")
#biocLite("DiffBind") # install ChIPQC

library(DiffBind)
library(edgeR)
library(ChIPseeker)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)
library(tools)
library(rlist)
library(qdapRegex)
library(valr)
#######################################################################################################################
#load bams and check headers: PE bams with duplicates removed
#######################################################################################################################

# setwd("/Volumes/ANNIE_DATA/ATACseq_4_2018/sorted/")
# 
# dataDir <- "/Volumes/ANNIE_DATA/ATACseq_4_2018/sorted"
# bamFiles <- dir(file.path(dataDir), pattern="*.bam$", full.name=T)
# 
# bamFiles
# bfList <- BamFileList(bamFiles)
# bfList
# path(bfList)
# #check sam headers
# for (i in 1:length(bamFiles)){
#   samHeader <- scanBamHeader(bamFiles[i])
#   str(samHeader,max.level=2)
# }
#   


#######################################################################################################################
#load sample sheet for first pair data for QC
#######################################################################################################################

#server:
#samples <- read.csv("/share/lasallelab/Annie/BMDM/ATACanalysis2019/DiffBind/sampleinfoBMDM_4_27_19.csv")

setwd("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/DiffBind")
samples <- read.csv("sampleinfoBMDM_4_27_19.csv")

#######################################################################################################################
#DiffBind Analysis
#######################################################################################################################

# calculate a binding matrix with scores based on read counts for every
# sample (affinity scores), rather than confidence scores for only those peaks called in a specific
# sample (occupancy scores). These reads are obtained using the dba.count function.1 As this
# example is based on a transcription factor that binds to the DNA, resulting in "punctate",
# narrow peaks, it is advisable to use the "summits" option to re-center each peak around
# the point of greatest enrichment. This keeps the peaks at a consistent width (in this case,
# with summits=250, the peaks will be 500bp, extending 250bp up- and down- stream of the summit):
# This shows that all the samples are using the same, 2845 length consensus peakset. Also, a
# new column has been added, called FRiP, which stands for Fraction of Reads in Peaks. This
# is the proportion of reads for that sample that overlap a peak in the consensus peakset, and
# can be used to indicate which samples show more enrichment overall.                                                                              

#PE data: https://support.bioconductor.org/p/53291/
# if the bLowMem parameter is set to TRUE in dba.count,
# DiffBind will use summarizeOverlaps. Changing the config value
# DBA$config$singleEnd to FALSE allows it to use paired-end BAM files.
# The documentation for summarizeOverlaps in the GenomicRanges package
# explains how it handles paired-end data.

#set minOverlap = 3 >only include peaks in at least this many peaksets when generating consensus peakset 
#set bScaleControl =FALSE >logical indicating if the Control reads should be scaled based on relative library sizes. If TRUE, and there are more reads in the Control library than in the ChIP library, the number of Control reads for each peak will be multiplied by a scaling factor determined by dividing the total number of reads in the ChIP library by the total number of reads in the C
#set bUseSummarizeOverlaps = TRUE >ogical indicating that summarizeOverlaps should be used for counting instead of the built-in counting code. This option is slower but uses the more standard counting function. If TRUE, all read files must be BAM (.bam extension), with associated index files (.bam.bai extension). The insertLength parameter must absent.
#set DBA$config$singleEnd =FALSE > true for single end
#set bParallel = TRUE >use multicore to get counts for each read file in parallel

#Likewise, if you are using summarizeOverlaps to count PE data, you should set the configuration parameter prior to the call to dba.count():
#myDBA$config$singleEnd <- FALSE
#myDBA <- dba.read(myDBA, bUseSummarizeOverlaps=TRUE)

dba <- dba(sampleSheet=samples)

#save(dba,file = "/share/lasallelab/Annie/BMDM/ATACanalysis2019/DiffBind/BMDM_ATACDBA.RData")

# correlation heatmap can be generated which gives an initial clustering of the samples using the cross-correlations of each row of the binding matrix
pdf("/share/lasallelab/Annie/BMDM/ATACanalysis2019/DiffBind/CCheatmap_AllPossiblepeaks.pdf",width=8, height=8,useDingbats=FALSE)
plot(dba)
dev.off()

pdf("/share/lasallelab/Annie/BMDM/ATACanalysis2019/DiffBind/PCA_AllPossiblepeaks.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotPCA((dba),DBA_CONDITION,label=DBA_ID,labelSize=1)
dev.off()


#for PE data:
dba$config$singleEnd <- FALSE

#raw and TMM CPM counts for downstream analysis
dba_count_raw <- dba.count(dba,minOverlap = 3, bScaleControl=FALSE,summits=75, score=DBA_SCORE_READS)
dba_count_CPM <- dba.count(dba,minOverlap = 3, bScaleControl=FALSE,summits=75, score=DBA_SCORE_TMM_READS_FULL_CPM)

save(samples,dba,dba_count_raw,dba_count_CPM,file = "BMDM_ATACDBA.RData")
load("BMDM_ATACDBA.RData")

pdf("initialheatmap_Allpossible_peaks.pdf",width=8, height=8,useDingbats=FALSE)
plot(dba_count_CPM)
dev.off()


#PCA normalized read counts for all the binding sites
pdf("PCA_normcounts_Allpossible_peaks.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotPCA(dba_count_CPM,DBA_CONDITION,label=DBA_ID,labelSize=1)
dev.off()


#write out consensus peaks to file for background use
allpeaks_gr <- dba.peakset(dba_count_CPM,bRetrieve=T)
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(allpeaks_gr),
                 starts=start(allpeaks_gr)-1,
                 ends=end(allpeaks_gr))

write.table(df, file="Allpossible_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


#######################################################################################################################
#consensus peaks across Replicates
#######################################################################################################################

#The overlap rate for samples can be isolated using a sample mask. A set
#of sample masks are automatically associated with a DBA object in the $masks field:

names(dba$masks)

#The returned data in olap.rate is a vector containing the number of peaks that appear in at least one, two, three, and so on up to all eleven peaksets.
olap.rate <- dba.overlap(dba, mode=DBA_OLAP_RATE)


pdf("Overlaps_Allpossible_peaks.pdf",width=8, height=8,useDingbats=FALSE)
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
dev.off()

#Arbitrary masks can be generated using the dba.mask function, or simply by specifying a vector of peakset numbers. 
#Value is a vector whose length is the number of peaksets, containing the number of overlapping peaks at the corresponding minOverlaps threshold 
#(i.e., Value[1] is the total number of unique sites, 
#Value[2] is the number of unique sites appearing in at least two peaksets, Value[3] the number of sites overlapping in at least three peaksets, etc.).

pdf("C57media_overlap_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(dba,dba$masks$`C57 media`,main="C57 BMDM media Open Chromatin Peaks")
dev.off()

pdf("C57LPSx1_overlap_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(dba,dba$masks$`C57 LPS1`,main="C57 BMDM LPSx1 Open Chromatin Peaks")
dev.off()

pdf("C57LPSx2_overlap_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(dba,dba$masks$`C57 LPS2`,main="C57 BMDM LPSx2 Open Chromatin Peaks")
dev.off()

pdf("BTBRmedia_overlap_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(dba,dba$masks$`BTBR media`,main="BTBR BMDM media Open Chromatin Peaks")
dev.off()

pdf("BTBRLPSx1_overlap_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(dba,dba$masks$`BTBR LPS1`,main="BTBR BMDM LPSx1 Open Chromatin Peaks")
dev.off()

pdf("BTBRLPSx2_overlap_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(dba,dba$masks$`BTBR LPS2`,main="BTBR BMDM LPSx2 Open Chromatin Peaks")
dev.off()


#######################################################################################################################
#consensus peaks 
#######################################################################################################################
# When forming the global binding matrix consensus peaksets, DiffBind first identifies all unique
# peaks amongst the relevant peaksets. As part of this process, it merges overlapping peaks,
# replacing them with a single peak representing the narrowest region that covers all peaks
# that overlap by at least one base. There are at least two consequences of this that are worth
# noting.
# First, as more peaksets are included in analysis, the average peak width tends to become
# longer as more overlapping peaks are detected and the start/end points are adjusted outward
# to account for them. Secondly, peak counts may not appear to add up as you may expect
# due to merging. For example, if one peakset contains two small peaks near to each other,
# while a second peakset includes a single peak that overlaps both of these by at least one
# base, these will all be replaced in the merged matrix with a single peak. A s more peaksets
# are added, multiple peaks from multiple peaksets may be merged together to form a single,
# wider peak. Use of the "summits" parameter is recommended to control for this widening
# effect.


#adds a new consensus peakset for each set of samples that share the same Tissue and Condition values.
#peak must be in 3/4 replicats for each tissu x treatment combination
consensus_peaks <- dba.peakset(dba, consensus=c(DBA_TISSUE,DBA_TREATMENT), minOverlap=0.75)

#save(consensus_peaks,file ="consensuspeaks.Rdata")
#setwd("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/atac\ analysis\ 2018/Diffbind2")
load("consensuspeaks.Rdata")
dba.show(consensus_peaks,mask=consensus_peaks$masks$Consensus)

# ID Tissue Factor  Condition Treatment  Replicate Caller Intervals
# 23  BTBR:LPS1   BTBR    1-2  BTBR LPS1      LPS1  5-8-14-21    bed     18097
# 24  BTBR:LPS2   BTBR    1-2  BTBR LPS2      LPS2  6-9-15-22    bed     16966
# 25 BTBR:media   BTBR    1-2 BTBR media     media  4-7-13-20    bed     22213
# 26   C57:LPS1    C57    1-2   C57 LPS1      LPS1 2-11-16-18    bed     23506
# 27   C57:LPS2    C57    1-2   C57 LPS2      LPS2    3-12-19    bed     20124
# 28  C57:media    C57    1-2  C57 media     media    1-10-17    bed      9860

#C57 overlaps of consensus peaks:
pdf("C57_consensuspeaks_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(consensus_peaks,26:28)
dev.off()

#BTBR overlaps of consensus peaks:
pdf("BTBR_consensuspeaks_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(consensus_peaks,23:25)
dev.off()

#media C57vsBTBR
pdf("media_consensuspeaks_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(consensus_peaks,consensus_peaks$masks$`BTBR media` & consensus_peaks$masks$`C57 media` )
dev.off()

#LPS1 C57vsBTBR
pdf("LPS1_consensuspeaks_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(consensus_peaks,consensus_peaks$masks$LPS1 & consensus_peaks$masks$`Replicate.1-2-3-4`)
dev.off()

#LPS2 C57vsBTBR
pdf("LPS2_consensuspeaks_venn.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotVenn(consensus_peaks,consensus_peaks$masks$LPS2 & consensus_peaks$masks$`Replicate.1-2-3-4`)
dev.off()

# create a new DBA object by adding three peaksets
#retrieve the consensus peakset as RangedData object
C57media_consensus_gr <- dba.peakset(consensus_peaks,consensus_peaks$masks$`Replicate.1-10-17` &consensus_peaks$masks$`C57 media`,bRetrieve=TRUE)

#write to bed
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(C57media_consensus_gr),
                 starts=start(C57media_consensus_gr)-1,
                 ends=end(C57media_consensus_gr),
                 names=c(rep("C57media_peaks", length(C57media_consensus_gr))),
                 scores=elementMetadata(C57media_consensus_gr)$C57.media,
                 strands=strand(C57media_consensus_gr))

write.table(df, file="C57media_consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


#retrieve the consensus peakset as RangedData object
C57LPSx1_consensus_gr <- dba.peakset(consensus_peaks,consensus_peaks$masks$`Replicate.2-11-16-18` &consensus_peaks$masks$`C57 LPS1`,bRetrieve=TRUE)

#write to bed
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(C57LPSx1_consensus_gr),
                 starts=start(C57LPSx1_consensus_gr)-1,
                 ends=end(C57LPSx1_consensus_gr),
                 names=c(rep("C57LPSx1_peaks", length(C57LPSx1_consensus_gr))),
                 scores=elementMetadata(C57LPSx1_consensus_gr)$C57.LPS1,
                 strands=strand(C57LPSx1_consensus_gr))

write.table(df, file="C57LPSx1_consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

#retrieve the consensus peakset as RangedData object
C57LPSx2_consensus_gr <- dba.peakset(consensus_peaks,consensus_peaks$masks$`Replicate.3-12-19` &consensus_peaks$masks$`C57 LPS2`,bRetrieve=TRUE)

#write to bed
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(C57LPSx2_consensus_gr),
                 starts=start(C57LPSx2_consensus_gr)-1,
                 ends=end(C57LPSx2_consensus_gr),
                 names=c(rep("C57LPSx2_peaks", length(C57LPSx2_consensus_gr))),
                 scores=elementMetadata(C57LPSx2_consensus_gr)$C57.LPS2,
                 strands=strand(C57LPSx2_consensus_gr))

write.table(df, file="C57LPSx2_consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


#retrieve the consensus peakset as RangedData object
BTBRmedia_consensus_gr <- dba.peakset(consensus_peaks,consensus_peaks$masks$`Replicate.4-7-13-20` &consensus_peaks$masks$`BTBR media`,bRetrieve=TRUE)

#write to bed
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(BTBRmedia_consensus_gr),
                 starts=start(BTBRmedia_consensus_gr)-1,
                 ends=end(BTBRmedia_consensus_gr),
                 names=c(rep("BTBRmedia_peaks", length(BTBRmedia_consensus_gr))),
                 scores=elementMetadata(BTBRmedia_consensus_gr)$BTBR.media,
                 strands=strand(BTBRmedia_consensus_gr))

write.table(df, file="BTBRmedia_consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


#retrieve the consensus peakset as RangedData object
BTBRLPSx1_consensus_gr <- dba.peakset(consensus_peaks,consensus_peaks$masks$`Replicate.5-8-14-21` &consensus_peaks$masks$`BTBR LPS1`,bRetrieve=TRUE)

#write to bed
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(BTBRLPSx1_consensus_gr),
                 starts=start(BTBRLPSx1_consensus_gr)-1,
                 ends=end(BTBRLPSx1_consensus_gr),
                 names=c(rep("BTBRLPSx1_peaks", length(BTBRLPSx1_consensus_gr))),
                 scores=elementMetadata(BTBRLPSx1_consensus_gr)$BTBR.LPS1,
                 strands=strand(BTBRLPSx1_consensus_gr))

write.table(df, file="BTBRLPSx1_consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

#retrieve the consensus peakset as RangedData object
BTBRLPSx2_consensus_gr <- dba.peakset(consensus_peaks,consensus_peaks$masks$`Replicate.6-9-15-22` &consensus_peaks$masks$`BTBR LPS2`,bRetrieve=TRUE)

#write to bed
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(BTBRLPSx2_consensus_gr),
                 starts=start(BTBRLPSx2_consensus_gr)-1,
                 ends=end(BTBRLPSx2_consensus_gr),
                 names=c(rep("BTBRLPSx2_peaks", length(BTBRLPSx2_consensus_gr))),
                 scores=elementMetadata(BTBRLPSx2_consensus_gr)$BTBR.LPS2,
                 strands=strand(BTBRLPSx2_consensus_gr))

write.table(df, file="BTBRLPSx2_consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


##################################################################################################################
#count in just the consensus peaks for all datasets in original DBA 
##################################################################################################################
#https://support.bioconductor.org/p/78918/
#The result of this is a binding affinity matrix containing a (normalized) read count for each sample at every potential binding site.

dba$config$singleEnd <- FALSE

#all consensus peaks:
consensus_allcond <- dba.peakset(consensus_peaks, minOverlap=1,bRetrieve=TRUE)

#write out consensus peaks to file
#subtract 1 to put in bed coordinates:
df <- data.frame(seqnames=seqnames(consensus_allcond),
                 starts=start(consensus_allcond)-1,
                 ends=end(consensus_allcond))

write.table(df, file="Allconsensuspeaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


#remove meta data
mcols(consensus_allcond) <- NULL


dba_count_consensus <- dba.count(dba,peaks = consensus_allcond, bScaleControl=FALSE,summits=TRUE, score=DBA_SCORE_TMM_READS_FULL_CPM)

#dba_count_consensus <- dba.count(dba, peaks = consensus_allcond, bScaleControl=FALSE,score=DBA_SCORE_TMM_READS_FULL_CPM,bParallel=TRUE)
save(consensus_peaks,dba_count_consensus,consensus_allcond,file ="consensuspeaks.Rdata")

pdf("heatmap_allconsensuspeaks.pdf", width=8, height=8, useDingbats=FALSE)
plot(dba_count_consensus,label=DBA_CONDITION)
dev.off()

#PCA normalized read counts for  consensuspeaks
pdf("PCA_normcountsAllconsensuspeaks.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotPCA(dba_count_consensus,DBA_CONDITION,label=DBA_ID,labelSize=1)
dev.off()


#######################################################################################################################
#contrasts and DE analysis using edgeR on raw counts normalized for total library size in original BAM input
#######################################################################################################################
#load data
#load("consensuspeaks_counts.Rdata")
#load("BMDM_ATACDBA.RData")

#library(genomation)
library(edgeR)
library(limma)
library(DiffBind)
#consensus_allcond <- readBed("Allconsensuspeaks.bed")

#get peak counts for all sampels in original dba
peakcounts <- dba
peaks <- peakcounts$peaks
peaknumber <- rapply(peaks, length, how="list")
peaknumber <- as.data.frame(do.call(rbind,peaknumber))
peaknumber <- unlist(peaknumber$V1)

peakcounts_out <- as.data.frame(cbind(peakcounts$samples[,1:6],peaknumber))

write.xlsx(peakcounts_out, file ="Peaks_per_sample.xlsx")
write.csv(peakcounts_out, file ="Peaks_per_sample.csv")

#get Frip for all consensus peaks
dba_count_consensus_raw <- dba.count(dba,peaks = consensus_allcond, bScaleControl=FALSE,summits=TRUE, score=DBA_SCORE_READS)
save(dba_count_consensus_raw,file="dba_count_consensus_raw.Rdata")
load("dba_count_consensus_raw.Rdata")

Frip <- dba_count_consensus_raw
Frip_out <- as.data.frame(cbind(Frip$samples[,1:6],Frip$SN))

write.xlsx(Frip_out, file ="Frip_allconsensuspeaks.xlsx")
write.csv(Frip_out, file ="Frip_allconsensuspeaks.csv")
#expore counts to datafile
#https://support.bioconductor.org/p/71301/
#raw scores: DBA_SCORE_READS
counts <- dba.peakset(dba_count_consensus_raw, bRetrieve=TRUE, score=DBA_SCORE_READS, DataType=DBA_DATA_FRAME)
#write.table(counts,"RAWcounts_consensuspeaks.txt",col.names = T, row.names = F, quote=F)

rownames(counts) <- paste(counts$CHR,counts$START,counts$END,sep="_")
counts <- counts[,4:ncol(counts)]
save(counts,file="Rawcounts_consensuspeaks.Rdata")
load("Rawcounts_consensuspeaks.Rdata")
#######################################################################################################################
#get library sizes from samtools flagstat and import
#######################################################################################################################

#librarysized calcualted on NFR Mito removed bam files using samtools flagstat
path2 <- "/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/QC/NFR_flagstat"
setwd(path2)
#get input files:
input.files <-dir(path2, pattern ='*NFRflagstat.txt')

#read in files and store information in a list named for each file
library(stringr)
depth <- NULL
for(i in 1:length(input.files)){
  file <- read.table(input.files[i],header=F, sep="\t", stringsAsFactors=FALSE)
  reads <- file[1,]
  reads <- stringi::stri_extract_first_words(reads)
  reads <- as.numeric(reads)
  filename <-word(input.files[i],1,sep = "_")
  tmp <- cbind(filename,reads)
  
  depth <- rbind(tmp, depth)
}
depth <- as.data.frame(depth)
depth$reads <- as.numeric(as.character(depth$reads))/2 #for paired end reads

setwd("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/DiffBind")

#reorder depth by columns in samples
depth2 <- merge(samples,depth, by.x="Replicate",by.y="filename")
depth2$SampleID <- as.character(depth2$SampleID)

target <- colnames(counts)
depth3 <- depth2[match(target, depth2$SampleID),]

write.csv(depth3,"FilteredLibrarySizes.csv")
#######################################################################################################################
#dgelist object with raw counts for NFR peaks with library size set for entire filtered bam file
dgList <- DGEList(counts=counts, genes=rownames(counts),lib.size =depth3$reads )

# Get log2 counts per million before RLE normalized data
logcounts <- cpm(dgList,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million Pre Normalization",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs Pre Normalization")
ggsave("BoxplotlogCPMsPreNormalization.pdf", w=8,h=5)

#Normalize by RLE : https://davetang.org/muse/2011/01/24/normalisation-methods-for-dge-data/
#“relative log expression”, as median library is calculated from the geometric mean of all columns
#and the median ratio of each sample to the median library is taken as the scale factor.
#set library size to full sequencing depth size
#also tried TMM but over-normalized
dgList <- calcNormFactors(dgList,method="RLE",lib.size =depth3$reads)

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(dgList$samples$lib.size,names=colnames(dgList),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# Get log2 counts per million from RLE normalized data
logcounts <- cpm(dgList,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs RLE Normalized")
ggsave("BoxplotlogCPMsPostNormalization.pdf", w=8,h=5)

#MDS
group <- word(rownames(dgList$samples),1,sep = "_")
group <- as.numeric(as.factor(group))
plotMDS(dgList, col=group)
ggsave("MDSplotPostNormalization.pdf", w=5,h=5)


# Get counts per million from RLE normalized data
CPM_RLEcounts <- cpm(dgList,log=FALSE)
CPM_RLEcounts <- as.data.frame(CPM_RLEcounts)
save(CPM_RLEcounts ,file="CPM_RLEcounts.RData")
load("CPM_RLEcounts.RData")
#######################################################################################################################
#setup model
#######################################################################################################################

#genotype:
#group <- word(rownames(dgList$samples),1,sep = "_")

Genotype <- rep("C57", ncol(dgList)) #C57 or BTBR
Genotype[grep("BTBR", colnames(dgList))] <- "BTBR"

#treatment:
Treatment <- rep("media", ncol(dgList)) #media, LPS1 or LPS2
Treatment[grep("LPS1", colnames(dgList))] <- "LPS1"
Treatment[grep("LPS2", colnames(dgList))] <- "LPS2"

#batch > can't use in model b/c redundant with sampleReplicate
#https://support.bioconductor.org/p/68092/
batch <- rep("batch1", ncol(dgList)) #1 = repl1 and 2, 2= repl3 and 4
batch[grep("_3", colnames(dgList))] <- "batch2"
batch[grep("_4", colnames(dgList))] <- "batch2"

#sample replicates
sampleReplicate <- rep("S1", ncol(dgList))
sampleReplicate[grep("_1", colnames(dgList))] <- "S1"
sampleReplicate[grep("_2", colnames(dgList))] <- "S2"
sampleReplicate[grep("_3", colnames(dgList))] <- "S3"
sampleReplicate[grep("_4", colnames(dgList))] <- "S4"

sampleReplicate <- paste(Genotype,sampleReplicate)
#Section 9.7 "Multi-level Experiments" of the limma User's Guide
#https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

#to make comparisons both within and between subjects
#make combined factor
y = dgList$samples 
d <- data.frame(Sample=rownames(y),Genotype,Treatment,sampleReplicate,batch)

Treat <- factor(paste(d$Genotype,d$Treatment,sep="."))
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)

#######################################################################################################################
#EdgeR analysis
#######################################################################################################################

#calculate dispersion = need all three for complex designs
dgList <- estimateGLMCommonDisp(dgList, design=design)
dgList <- estimateGLMTrendedDisp(dgList, design=design)
dgList <- estimateGLMTagwiseDisp(dgList, design=design)

#plot dispersion
plotBCV(dgList)

#https://support.bioconductor.org/p/59700/
#What the voom function does is essentially to provide weights for the regression fit. 
#Now, in a situation where you want to estimate a correlation you first obtain "working weights"
#that are estimated WITHOUT the correlation. 
#You use them to estimate the correlation and then use the voom function again, 
#this time taking the estimated correlation into account in the estimation of the  weights.
#v <- voom(dgList, design, plot=TRUE)
#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
# voom
v <- voom(dgList, design, plot=T) 
corfit <- duplicateCorrelation(v, design, block = d$sampleReplicate) #find weights due to sample correlation
corfit$consensus.correlation #0.1847

fit <- lmFit(v, design, block = d$sampleReplicate, correlation = corfit$consensus.correlation)

save(counts, depth, dgList,v,fit,design,corfit, file="EdgeRlimmaVoomModles.RData")
#box plots for the voom normalised data to compare to before normalisation (only RLE)
#v$E already log2CPM

boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2, ylim = c(-5, 15))
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Voom transformed logCPM")
ggsave("BoxplotlogCPMsPostNormalizationVoom.pdf", w=8,h=5)


#get log2CPM counts from voom and put in dataframe:
library(plotrix)
#average log2 CPM and sem
countdf <- as.data.frame(v$E)
countdf$regions <- rownames(counts)
countdf2 <- countdf %>% group_by(regions) %>% gather(condition,log2CPM, 1:22) %>% 
  separate(condition,into=c("treatment","replicate"),sep="_") %>%
  group_by(regions,treatment) %>% summarize(meanlog2CPM = mean(log2CPM),SEM = std.error(log2CPM))

save(countdf,countdf2,file="RLE&voomLog2CPM.RData")

#build contrast matrix for comparisons of interest:
cont_matrix <- makeContrasts(BTBRmediavsC57media = BTBR.media - C57.media,
                             BTBRLPS1vsBTBRmedia = BTBR.LPS1 - BTBR.media, 
                             BTBRLPS2vsBTBRmedia = BTBR.LPS2 - BTBR.media, 
                             C57LPS1vsC57media = C57.LPS1 - C57.media, 
                             C57LPS2vsC57media = C57.LPS2 - C57.media,
                             BTBRLPS1vsC57LPS1 = BTBR.LPS1 - C57.LPS1,
                             BTBRLPS2vsC57LPS2 = BTBR.LPS2 - C57.LPS2,
                             levels = design)


##loop for each comparison:
comparisons=colnames(cont_matrix)

comp_out <- NULL
comp_out_sig <- NULL
for(i in 1:length(comparisons)){
  #comparison name
  comp=comparisons[i]
  #make comparisons 
  diff <- contrasts.fit(fit,contrast=cont_matrix[,comp])
  tmp <- eBayes(diff)
  tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
  tmp2$comparison <- paste(comp)
  
  pdf(file = paste(comp,"_Volcano.pdf", sep=""), wi = 9, he = 6, useDingbats=F)
  
  #volcano plot
  with(tmp2, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(tmp2, adj.P.Val<0.05), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  # Label points with the textxy function from the calibrate plot
  #library(calibrate)
  #with(subset(tmp2, adj.P.Val<0.05 & abs(logFC)>10), textxy(logFC., -log10(adj.P.Val), labs=Column.ID, cex=.8))
  dev.off()
  
  
  #get conditions
  condition <- unique(comp)
  condition1 <- strsplit(condition, "vs")[[1]][1]
  condition2 <- strsplit(condition, "vs")[[1]][2]
  
  targets <- c(condition1,condition2)
  #get means and sem
  counttmp1 <- countdf2[(countdf2$treatment %in% targets),]
  counttmp <- counttmp1 %>% dplyr::select(-SEM) %>% group_by(regions) %>% spread(treatment,meanlog2CPM)
  counttmp <- as.data.frame(counttmp)
  counttmp$CPM_cond1 <- 2^counttmp[,2]
  counttmp$CPM_cond2 <- 2^counttmp[,3]
  
  counttmp$CPMdifference <- counttmp$CPM_cond1 - counttmp$CPM_cond2
  
  #merge data
  mergedf <- merge(tmp2,counttmp, by.x="genes",by.y="regions")
  
  
  write.csv(mergedf,file=paste(comp,"allpeaksdata.csv"))
  
  #change colnames for output
  mergedf2 <- mergedf
  colnames(mergedf2)[9] <- c("Log2CPMcondition1")
  colnames(mergedf2)[10] <- c("Log2CPMcondition2")
  comp_out <- rbind(comp_out,mergedf2)
  
  
  #write to bed file for each set of significant peaks for incr or decr
  bed <- mergedf %>% separate(genes,into = c("chr","start","end"),sep="_") %>% filter(adj.P.Val < 0.05)
  
  bed$direction <- c("null")
  
  bed$direction[which(bed$CPMdifference<0.00)] <- paste(condition1,"<",condition2, sep=" ")
  bed$direction[which(bed$CPMdifference>0.00)] <- paste(condition1,">",condition2, sep=" ")
  
  write.csv(bed,file=paste(comp,"sigpeaksdata.csv"))
  
  #change col names for binding
  
  bed2 <- bed
  colnames(bed2)[11] <- c("Log2CPMcondition1")
  colnames(bed2)[12] <- c("Log2CPMcondition2")
  comp_out_sig <- rbind(comp_out_sig,bed2)
  
  
  df_up_bed <- bed %>% filter(CPMdifference < 0.00) %>% dplyr::select(c(1:3,10))
  df_down_bed <- bed %>% filter(CPMdifference > 0.00) %>% dplyr::select(c(1:3,10))
  
  write.table(df_up_bed, file=paste(condition1,"<",condition2,"DEpeaks.mm10.bed" ,sep="_"),col.names = F, row.names = F, quote=F, sep="\t")
  write.table(df_down_bed, file=paste(condition1,">",condition2,"DEpeaks.mm10.bed" ,sep="_"),col.names = F, row.names = F, quote=F, sep="\t")
  
}

#write out all peaks
outallpeaks <- comp_out %>% separate(genes,into = c("chr","start","end"),sep="_")
write.csv(outallpeaks, file="Allpeaks_alldata.csv")

#bed
allpeaks <- outallpeaks[,c(1:3,10)]
write.table(allpeaks, file="Allpeaks_alldata.bed",col.names = F, row.names = F, quote=F, sep="\t")

#write out all sig peaks
write.csv(comp_out_sig, file="AllSignificantpeaks_alldata.csv")
write.table(comp_out_sig[,c(1:3,10)], file="AllSignificantpeaks_alldata.bed",col.names = F, row.names = F, quote=F, sep="\t")


#######################################################################################################################
#summary table of significant peaks per condition
#######################################################################################################################

#43434
AllSignificantpeaks <- comp_out_sig[,c(1:3,10)]

#make GRanges object for each bed and save to list
path <- getwd()

#https://stackoverflow.com/questions/27911604/loop-through-all-files-in-a-directory-read-and-save-them-in-an-object-r


files <- list.files(path=".", pattern="*.mm10.bed", all.files=T, full.names=T)
filelist <- lapply(files, readPeakFile)


#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files)))
names_list <- gsub("_DEpeaks.mm10", "", names_list)
names_list <- gsub("_", " ", names_list)

names(filelist) <- names_list

#times in each list within the filelist
number_regions <- lapply(filelist, function(x) length(x))

number_regions <- as.data.frame(do.call(rbind, number_regions))

colnames(number_regions) <- c("count")
number_regions$peaklist <- rownames(number_regions)

openxlsx::write.xlsx(number_regions, file="Counts_per_DEpeaklist.xlsx")



#######################################################################################################################
#counts for all significant peaks in comparisons:
#######################################################################################################################

#load("BMDM_ATACDBA.RData")

AllSignificantpeaks <- readPeakFile("AllSignificantpeaks_alldata.bed")

#summit height (maximum read pileup value), normalized to relative library size
dba_count_sig <- dba.count(dba,peaks = AllSignificantpeaks, bScaleControl=FALSE,summits=TRUE, score=DBA_SCORE_SUMMIT_ADJ)

save(AllSignificantpeaks,dba_count_sig,file ="BMDM_ATAC_SigPeakCount.RData")

pdf("heatmap_significnatconsensuspeaks_summitheight.pdf",width=8, height=8,useDingbats=FALSE)
plot(dba_count_sig)
dev.off()

#PCA normalized read counts for significnatconsensuspeaks
pdf("PCA_significnatconsensuspeaks_summitheight.pdf",width=8, height=8,useDingbats=FALSE)
dba.plotPCA(dba_count_sig,DBA_CONDITION,label=DBA_ID,labelSize=1)
dev.off()

#######################################################################################################################
#look to see if BTBR are > C57 at all DE peaks
#######################################################################################################################
setwd("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/peakfiles_5_5_19")

load("BMDM_ATAC_SigPeakCount.RData")

files <- list.files(path=".", pattern="*.mm10.bed", all.files=T, full.names=T)
filelist <- lapply(files, readPeakFile)

#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files)))
names_list <- gsub(".mm10", "", names_list)
names_list <- gsub(" ", "", names_list)

names(filelist) <- names_list

save(filelist,file="peaklists.RData")

#summit height (maximum read pileup value), normalized to relative library size
#dba_count_sig_BTBRmediaLC57media <- dba.count(dba,peaks = filelist[9], bScaleControl=FALSE,summits=TRUE, score=DBA_SCORE_SUMMIT_ADJ)

setwd("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/DiffBind")

#run on server:
for (i in 1:length(names(filelist))) {
  print(i)
  #summit height (maximum read pileup value), normalized to relative library size
  dba_count_sig2 <- dba.count(dba,peaks = filelist[[i]], bScaleControl=FALSE,summits=TRUE, score=DBA_SCORE_SUMMIT_ADJ)
  
  #write out consensus peaks to file for background use
  peaks <- dba.peakset(dba_count_sig2,bRetrieve=T, numCols=6)
  
  filename <- gsub(" ","",names(filelist[i]))
  
  write.csv(peaks,paste(filename, "_peaksummitheights.csv", sep=""))
}

files <- list.files(path=".", pattern="*_peaksummitheights.csv", all.files=T, full.names=F)

peakcount <- NULL
for (i in 1:length(files)) {
  
  tmp <- read.csv(files[i])
  
  name <- paste(files[i])
  name <- gsub("_peaksummitheights.csv","",name)
  
  tmp$peakset <- name
  
  peakcount <- rbind(peakcount,tmp)
}


peakcount2 <- peakcount[,7:ncol(peakcount)]

library(dplyr)
library(tidyr)

peakcount2$peakset <- factor(peakcount2$peakset)
peakcount3 <- peakcount2 %>% group_by(peakset) %>%
  gather(Replicate,AdjSummitHeight,1:(ncol(peakcount2)-1)) %>%
  separate(Replicate,into=c("condition","animal"),sep="_")

library(ggplot2)
library(cowplot)

levels <- c("C57media","C57LPS1" ,  "C57LPS2" ,"BTBRmedia","BTBRLPS1" , "BTBRLPS2")
peakcount3$condition <- factor(peakcount3$condition,levels=levels)


peakcount3$peakset <- gsub("_", "", peakcount3$peakset)

neworder <- c("BTBRmedia<C57media" ,"BTBRmedia>C57media",
              "C57media<LPS1only","C57media<LPS2only", "C57media<LPS1&2",
              "BTBRmedia<LPS1only","BTBRmedia<LPS2only", "BTBRmedia<LPS1&2",
                 "C57media>LPS1only","C57media>LPS2only", "C57media>LPS1&2",
              "BTBRmedia>LPS1only","BTBRmedia>LPS2only", "BTBRmedia>LPS1&2",
              "BTBRLPS1<C57LPS1"  , "BTBRLPS1>C57LPS1" ,  "BTBRLPS2<C57LPS2" ,  "BTBRLPS2>C57LPS2" )

peakcount3$peakset <- factor(peakcount3$peakset, levels = neworder)


#boxplot
#cbPalette <- c("#999999", "#D55E00", "#0072B2")

# Function to produce summary statistics (mean and +/- sd)
# data_summary <- function(x) {
#   m <- mean(x)
#   ymin <- m-sd(x)
#   ymax <- m+sd(x)
#   return(c(y=m,ymin=ymin,ymax=ymax))
# }

pdf("SummitHeight_byPeakSet.pdf", height = 11, width =8.5)    # create PNG for the heat map       

ggplot(peakcount3, aes(x=condition, y=AdjSummitHeight,group=condition)) + 
  facet_wrap(~peakset, scales = "free")+
  #geom_violin(trim=FALSE)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
# stat_summary(fun.data=data_summary)+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=condition))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  #geom_point(position = position_dodge(width = 0.90),aes(group=condition)) + 
  #scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 12)+
  theme(strip.background=element_rect(fill="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Peak Set") + #,limits = order2)+
  scale_y_continuous(name = "Read Depth Adjusted Summit Height")+
  ggtitle(" ") +
  theme(legend.position="none") # Remove legend

dev.off()



########################################################################################
#loading data, compiling, cleaning
########################################################################################
#repeated measures ANOVA
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(lsmeans)
library(nlme)

########################################################################################
#ANOVA for peak sets (not averaged by animal)
########################################################################################


#peakcount3 <- peakcount2 %>% group_by(peakset) %>%
#  gather(Replicate,AdjSummitHeight,1:(ncol(peakcount2)-1)) %>%
#  separate(Replicate,into=c("condition","animal"),sep="_") %>%
#  group_by(peakset,condition,animal) %>%
#  summarize(mean = mean(AdjSummitHeight))

#set factors
#define 
peakcount3$strain <- substr(peakcount3$condition, start=1, stop=3)
peakcount3$strain <- gsub("BTB","BTBR",peakcount3$strain)
peakcount3$strain <- factor(peakcount3$strain, levels=c("C57","BTBR"))

peakcount3$treatment <- gsub("BTBR","",peakcount3$condition)
peakcount3$treatment <- gsub("C57","",peakcount3$treatment)
peakcount3$treatment <- factor(peakcount3$treatment, levels = c("media","LPS1","LPS2"))

#model


# models = lapply(setNames(vars, vars), function(var) {
#   form = paste(var," ~ treatment + strain, ~1|animal")
#   print(form)
#   #lme(form, data=peakcount3)
# })
vars = unique(peakcount3$peakset)

anova <- NULL
posthoc <- NULL
for (i in vars) {
  dat <- peakcount3 %>% filter(peakset == i)
  tmp <- lme(AdjSummitHeight~ treatment * strain, ~1|animal, data = dat)
  
  #model indcludes random effects subject ID 
  anova1 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
  anova1$peakset <- unique(paste(i))
  anova <- rbind(anova,anova1)
  
  #posthoc tests for group at each timepoint
  refgrid <- ref.grid(tmp)
  tmp2 <- lsmeans(refgrid, ~strain*treatment)
  tmp3 <- summary(pairs(tmp2, adjust = "none"))
  tmp3 <- as.data.frame(tmp3)
 # tmp3$contrast <- character(tmp3$contrast)
  #select: C57,LPS1 - BTBR,LPS1 and C57,LPS2 - BTBR,LPS2
  tmp4 <- dplyr::filter(tmp3, grepl('C57,LPS1 - BTBR,LPS1|C57,LPS2 - BTBR,LPS2', contrast))
  tmp4$peakset <- unique(paste(i))
  posthoc <- rbind(posthoc,tmp4)
}

posthoc$FDR <- p.adjust(posthoc$p.value, method="fdr")

anova2 <- anova

anova2$condition <- rownames(anova2)
anova2 <- anova2[!(grepl("(Intercept)",anova2$condition)),]
anova2$condition[(grepl("^strain",anova2$condition))] <- c("strain")


posthoc2 <- posthoc
posthoc2$comparison <- rownames(posthoc2) 

posthoc2$contrast <- gsub("LPS2","LPSx2",posthoc2$contrast)
posthoc2$contrast <- gsub("LPS1","LPSx1",posthoc2$contrast)

#write anova to output file
write.xlsx(anova2,file="RMANOVA_MeanAdjSummitHeight.xlsx")
write.xlsx(posthoc2,file="PostHocs_MeanAdjSummitHeight.xlsx")





