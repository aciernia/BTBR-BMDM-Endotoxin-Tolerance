##########################################################################
#AVC a.ciernia@gmail.com
#8/2018
#ATACseq BMDM +/- 1 or 2 doses LPS in vitro
##########################################################################
#Analysis of overlapping peaks from DiffBind analysis output and SNP/Indels from VCF
##########################################################################
library(tools)
library(valr)
library(regioneR)
library(ChIPseeker)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("regioneR", version = "3.8")
##########################################################################

#make GRanges object for each bed and save to list
#https://stackoverflow.com/questions/27911604/loop-through-all-files-in-a-directory-read-and-save-them-in-an-object-r

#######################################################################################################################
setwd("/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/consensus\ peaks")
path <- getwd()
#get lists
files <- list.files(path=".", pattern="*.mm10.bed", all.files=T, full.names=T)
filelist <- lapply(files, readPeakFile)

#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files)))
names_list <- gsub("_DEpeaks.mm10", "", names_list)
names_list <- gsub("_L_", "<", names_list)
names_list <- gsub("_G_", ">", names_list)

names(filelist) <- names_list
names(filelist)

#######################################################################################################################


# Overlap testing regioneR ---------------------------------------

setwd("/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/SNP_peakoverlaps")

#load BTBR SNPs and Indels
load("BTBR_snps_indels_data.Rdata")

#convert to GRanges
library(GenomicRanges)
library(dplyr)

BTBRvcf_snps$starts <- as.numeric(BTBRvcf_snps$starts)
BTBRvcf_snps$ends <- as.numeric(BTBRvcf_snps$ends)

#error: rror in .Call2("solve_user_SEW0", start, end, width, PACKAGE = "IRanges") : solving row 1: negative widths are not allowed
#solution:https://support.bioconductor.org/p/47332/
BTBRvcf_snps2 <- BTBRvcf_snps[,c(3,4,5,2,6,7,15,16,17)]
BTBRvcf_snps2 <- BTBRvcf_snps2 %>% distinct() %>%
  arrange(seqnames,starts,ends) %>%
  na.omit()
  #all.equal(BTBRvcf_snps, reduce(BTBRvcf_snps))

SNPs <- GRanges(seqnames = BTBRvcf_snps2$seqnames, ranges = IRanges(start = BTBRvcf_snps2$starts, end = BTBRvcf_snps2$ends), mcols = BTBRvcf_snps[,c(2,6,7,15,16,17)])

Indels <- GRanges(seqnames = BTBRvcf_indels$seqnames, ranges = IRanges(start = BTBRvcf_indels$starts, end = BTBRvcf_indels$ends), mcols = BTBRvcf_indels[,c(1,6:ncol(BTBRvcf_indels))])

#GRangesList

SNPList <- list(ALLSNPs = ALLSNPs, NonSynSNPs = NonSynSNPs, SynSNPs =SynSNPs, NonsenseSNPs = NonsenseSNPs)

#save overlap data
setwd("/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/SNP_peakoverlaps")
#save(Bfiles,SNPList,file="BTBRSNP_dataforoverlap.Rdata")
load("BTBRSNP_dataforoverlap.Rdata")

#replace old Bfiles with new ones
Peakfiles <- unlist(filelist[])

#define universe as all consesnsu peaks
universe <- readBed(file = "/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/DiffBind/Allconsensuspeaks.bed")

save(Peakfiles,SNPList,universe,file="BTBRSNP_dataforoverlap_1_5_20.Rdata")
load("/Users/aciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/SNP_peakoverlaps/BTBRSNP_dataforoverlap_1_5_20.Rdata")
#run on server:
#####################BTBR vs C57 comparisons##########################################
# Overlap testing regioneR ---------------------------------------
#overlap SNPs with 
# overlapping A with B > EVobs
# resample from universe (same size as list A)
# compare resampled list with list B > overlap of random sample with list B
# repeat 1000 times > distribution of EVperm
# gives permutation p value and Z-score

#get mouse genome
mouse.genome <- getGenomeAndMask("mm10", mask=NA)$genome
mouse.canonical <- filterChromosomes(mouse.genome, organism="mm")


#number of regions in each list
rl <- rapply(filelist, length, how="list")
rl <- as.data.frame(do.call(rbind, rl))
rl$comparison <- rownames(rl)
colnames(rl)[1] <- c("number of regions")
openxlsx::write.xlsx(rl,file="Number_of_regions.xlsx")

#load("BTBRSNP_dataforoverlap_4_30_19.Rdata")

#######################################################################################################################
#run on server
#######################################################################################################################

#get mouse genome
mouse.genome <- getGenomeAndMask("mm10", mask=NA)$genome
mouse.canonical <- filterChromosomes(mouse.genome, organism="mm")


#compare lists
set.seed(1)

require(regioneR)
output <- NULL
for (i in 1:length(Peakfiles)){
  sum <-NULL
  #loop through each file for comparison  
  for (B in 1:length(SNPList)) {
    pt <- permTest(A=Peakfiles[[i]],B=SNPList[[B]], ntimes=1000,
                   alternative="greater",
                   randomize.function=circularRandomizeRegions,
                   universe=universe,
                   evaluate.function=numOverlaps, 
                   count.once=TRUE,
                   genome="mm10", mc.set.seed=FALSE, mc.cores=10)
    summary <- summary(pt)
    summary$condition1 <- paste(names(Peakfiles[i]))
    summary$condition2 <- paste(names(SNPList[B]))
    summary$length_condition1 <- length(Peakfiles[[i]])
    summary$length_condition2 <- length(SNPList[[B]])
    summary$overlap <- pt$numOverlaps$observed
    
    sum <- rbind(sum,summary)
    #plot(pt)
    write.csv(summary,paste("output",paste(names(Peakfiles[i])),paste(names(SNPList[B])),"csv",sep="."))
  }
  
  output <- rbind(sum,output)
}

write.csv(output,"output_SNPpeakoverlap_4_22_2020.csv")

#######################################################################################################################

#load from server:
#######################################################################################################################
setwd("/Users/aciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/SNP_peakoverlaps/1_5_2020output")
output <- read.csv("output_SNPpeakoverlap_1_5_2020.csv")

#FDR correct comparisons
output$FDR <- p.adjust(output$pvalue, method="fdr")

openxlsx::write.xlsx(output, file = "SNPIndelPeakOverlap_RegioneR.xlsx")
#number of SNPs/Indels per peak
#load("BTBRSNP_dataforoverlap.Rdata")
library(GenomicRanges)
library(rlist)

x <- list()
output2 <- NULL
for (i in 1:length(SNPList)){
  sum <-NULL
  #loop through each file for comparison  
  for (B in 1:length(Peakfiles)) {
    gr1=SNPList[[i]]
    gr2=Peakfiles[[B]]
    
    condition1 = paste(names(SNPList[i]))
    condition2 = paste(names(Peakfiles[B]))
    #peaks with snps
    #ranges <- subsetByOverlaps(gr2,gr1)
    
    #snps inside peaks
    #ranges2 <- subsetByOverlaps(gr1,gr2)
    
    #find the overlaps
    #overlaps <- findOverlaps(gr2, gr1)
    
    count_peaks_with_snps <-  countOverlaps(gr2,gr1)
    
    allpeakscount <- length(count_peaks_with_snps)
    peakswithsnpscount <- length(count_peaks_with_snps[count_peaks_with_snps>0])
    peakswithoutsnpscount <- length(count_peaks_with_snps[count_peaks_with_snps==0])
    
    #save peaks with snps counts to a list
    count_peaks_with_snps <- count_peaks_with_snps[count_peaks_with_snps>0]
    count_peaks_with_snps2 <- list(count_peaks_with_snps)
    names(count_peaks_with_snps2) <- paste(condition1,condition2,sep="_")
    
    #append to list
    x <- list.append(x,count_peaks_with_snps2)
    
    summary <- data.frame(allpeakscount = allpeakscount,
                          peakswithVarcount = peakswithsnpscount,
                          peakswithoutVarcount = peakswithoutsnpscount,
                          condition1 = paste(names(SNPList[i])),
                          condition2 = paste(names(Peakfiles[B])))
    
    sum <- rbind(sum,summary)
    
    #write.csv(summary,paste("output",paste(names(SNPList[i])),paste(names(Bfiles[B])),"csv",sep="."))
  }
  
  output2 <- rbind(sum,output2)
}



#plot
library(ggplot2)
library(reshape2)
library(cowplot)
df <- melt(x)
df2 <- df %>% tidyr::separate(L2, into= c("Variant","RegionsList"), sep="_", remove=F)
df2$RegionsList <- factor(df2$RegionsList,levels= unique(df2$RegionsList))
df2$Variant <- factor(df2$Variant,levels= unique(df2$Variant))
#ggplot(df, aes(x = factor(L2), y = value)) + geom_point()
df2$value <- as.numeric(df2$value)
library(cowplot)
myplot <- ggplot(df2, aes(x = factor(RegionsList), y = value, fill=Variant)) +
  facet_grid(~Variant) + 
  #wiskers = 5th and 95th percentile.https://stackoverflow.com/questions/28961781/changing-whisker-definition-in-faceted-geom-boxplot
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position=position_dodge())+
  #geom_point(size = 1, position=position_dodge(width=0.9)) +
  scale_y_continuous(name="Variant Count within Regions") +
  scale_x_discrete(name = "Region List")+
  theme(legend.title = element_blank())+  
  #scale_fill_manual(values = cbPalette) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(axis.text.x=element_text(angle=65,hjust=1)) +
  theme(legend.key = element_blank(), 
        strip.background = element_blank()
       #panel.border = element_rect(colour = "black")
        ) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20),legend.text=element_text(size=20))

myplot

ggsave(myplot, file="VariantCount_perRegion.pdf",height = 10, width = 10)
```

```{r}
#combine outputs
colnames(output2)[4:5] <- c("condition2","condition1")
outputmerge <- merge(output,output2, by=c("condition1","condition2"))

openxlsx::write.xlsx(outputmerge, file = "SNPIndelPeakOverlap_RegioneR.xlsx")


outputmerge$Percent_regions_with_variants <- (outputmerge$peakswithVarcount/outputmerge$allpeakscount)*100


neworder <- c("BTBRmedia<C57media" ,"BTBRmedia>C57media",
              "C57LPS1>C57media","BTBRLPS1>BTBRmedia", "C57LPS2>C57media","BTBRLPS2>BTBRmedia",
              "C57LPS1<C57media","BTBRLPS1<BTBRmedia", "C57LPS2<C57media","BTBRLPS2<BTBRmedia",
              "BTBRLPS1<C57LPS1","BTBRLPS2<C57LPS2",
              "BTBRLPS1>C57LPS1","BTBRLPS2>C57LPS2")

outputmerge$condition1 <- factor(outputmerge$condition1 ,levels= neworder)
outputmerge$condition2 <- factor(outputmerge$condition2,levels= unique(df2$Variant))



myplot <- ggplot(outputmerge, aes(x = condition1, y = Percent_regions_with_variants, fill=condition2)) +
  facet_grid(~condition2) +
  geom_bar(stat="identity")+
  scale_y_continuous(name="Percent of Regions \n with Variants") +
  scale_x_discrete(name = "Region List")+
  theme(legend.title = element_blank())+  
  theme(strip.text.x = element_text(size = 20)) +
  theme(axis.text.x=element_text(angle=65,hjust=1)) +
  #theme(legend.key = element_blank(), strip.background = element_blank(),panel.border = element_rect(colour = "black"))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))+
  theme(legend.title = element_blank()) 

myplot

ggsave(myplot, file="PercentRegionswithVariants.pdf",height = 6, width = 12)

#######dotplot Percent of Regions with Variants#######

dp <- outputmerge
unique(dp$condition1)

dp$condition1 <- gsub("LPS1","LPSx1",dp$condition1)
dp$condition1 <- gsub("LPS2","LPSx2",dp$condition1)
dp$condition1 <- gsub("7LPSx2","7 LPSx2",dp$condition1)
dp$condition1 <- gsub("7LPSx1","7 LPSx1",dp$condition1)
dp$condition1 <- gsub("RLPSx2","R LPSx2",dp$condition1)
dp$condition1 <- gsub("RLPSx1","R LPSx1",dp$condition1)
dp$condition1 <- gsub("Rmedia","R media",dp$condition1)
dp$condition1 <- gsub("7media","7 media",dp$condition1)

dp$strainDAR <- dp$condition1 
dp$strainDAR[grepl("C57", dp$strainDAR) & grepl("BTBR", dp$strainDAR)] <- c("Strain DAR")
dp$strainDAR[!grepl("Strain DAR",dp$strainDAR)] <- c("Non-Strain DAR")
dp$Percent_regions_with_variants

dp <- dp %>% arrange(desc(strainDAR), condition1)
dp$label <- paste(dp$strainDAR, dp$condition1, sep=" ")
dp$label <- factor(dp$label)
  
dp$neglogFDR <- -log(dp$FDR)

cutoff= -log(0.05)

pdf("PercentRegionswithVariants.dotPlot2.pdf", width=8, height4)

ggplot(dp, aes(y = label, x =condition2,group=condition2)) +
  #facet_grid(~strainDAR, scales="free")+
  theme_bw() +
  geom_point( alpha=1, aes(color = neglogFDR,size=Percent_regions_with_variants),position=position_dodge(width = .5)) + 
  scale_color_gradient2(midpoint=cutoff, low="blue", mid="white",
                        high="red", space ="Lab" )+
 # scale_color_manual(values = c("% BG promoters with motif" = "darkgrey", "% DEG promoters with motif" = "blue")) +
  xlab("BTBR Genetic Variants") +
  ylab("DARs") +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

write.csv(dp, "PercentOverlapOutput_4_22_2020.csv")
