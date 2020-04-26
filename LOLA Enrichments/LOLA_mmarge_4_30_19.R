#author: Annie Vogel Ciernia
#a.ciernia@gmail.com
#10/9/2018
##############################################################################################################
library(dplyr)
library(tidyr)
library(cowplot)
library(gplots)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("LOLA")

library(LOLA)


#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)

#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)

options("scipen"=100, "digits"=4) #prevent exponents

##############################################################################################################
#read in peak bed files to GRanges List


library(ChIPseeker)
library("zonator")
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


#load background regions: all possible peaks called in all samples
Background <- readBed(file = "/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/DiffBind/Allconsensuspeaks.bed")

#load DB:
#regionDB <- loadRegionDB(dbLocation = "/Users/annieciernia/Desktop/regionDB/mm10/",limit = NULL, collections = c("collection1","collection3","collection4","collection5","collection6","collection7","collection8"))

#save(regionDB,file="/Users/annieciernia/Desktop/regionDB/RegionDBmm10_9_3_18.Rdata")
#support is the overlap, and b, c, and d complete the 2x2 table

regionDB_TF <- loadRegionDB(dbLocation = "/Users/annieciernia/Desktop/mm10_LOLA_DB/mm10/",limit = NULL, collections = c("collection5","collection6","collection7","collection8"))

setwd("/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/peak_overlaps1_5_2020")
#two tailed fisher exact test:
Results <- runLOLA(userSets = filelist, userUniverse = Background, regionDB = regionDB_TF, minOverlap = 1, cores=2, redefineUserSets = FALSE,direction = "enrichment")

#locResult = Results[2,]
#extractEnrichmentOverlaps(locResult, filelist, regionDB_TF)

writeCombinedEnrichment(combinedResults = Results, outFolder = "DEpeak_RegionOverlaps", includeSplits=F)

########################################################################################################################
merge <- Results

#pValueLog:=-log10(pValueLog + 10^-322)
#make pvalue
#merge$pvalue <- 10^-(merge$pValueLog)

#undo pseudo count:
# merge$pvalue <- merge$pvalue - 10^-322
# merge$pvalue <- abs(merge$pvalue)
# merge$FDR <- p.adjust(merge$pvalue,method = "fdr")

#merge <- enrichments
names(merge)[names(merge) == 'support'] <- 'userSet.in.target.list'
names(merge)[names(merge) == 'b'] <- 'NonuserSet.in.target.list'
names(merge)[names(merge) == 'c'] <- 'userSet.not.in.target.list'
names(merge)[names(merge) == 'd'] <- 'NonuserSet.not.in.target.list'

#% enrichment
merge$percent_userSet_in_Target <- (merge$userSet.in.target.list/(merge$userSet.in.target.list + merge$userSet.not.in.target.list)*100)
merge$percent_BG_in_Target <- (merge$NonuserSet.in.target.list/(merge$NonuserSet.in.target.list + merge$NonuserSet.not.in.target.list)*100)

#fold enrichment relative to background
#merge$FC <- (merge$percent_userSet_in_Target - merge$percent_BG_in_Target)/merge$percent_BG_in_Target  

#read in list descriptions
listdescript <- read.csv("descriptionorder_fixed.csv")

merge2 <- merge(merge,listdescript,by.x="description",by.y="old")

write.csv(merge2,file="FisherExact_TF_Enrichements_1_5_20.csv")


########################################################################################################################
#clean up names
enrichments <- merge2
unique(enrichments$userSet)

neworder <- c("BTBRmedia<C57media" ,"BTBRmedia>C57media",
              "C57LPS1>C57media","BTBRLPS1>BTBRmedia", "C57LPS2>C57media","BTBRLPS2>BTBRmedia",
              "C57LPS1<C57media","BTBRLPS1<BTBRmedia", "C57LPS2<C57media","BTBRLPS2<BTBRmedia",
              "BTBRLPS1<C57LPS1","BTBRLPS2<C57LPS2",
              "BTBRLPS1>C57LPS1","BTBRLPS2>C57LPS2")

enrichments$userSet <- factor(enrichments$userSet, levels = neworder)


#enrichments_sig <- filter(enrichments,enrichments$qValue<0.05)

########################################################################################################################
#plots
########################################################################################################################
#significant enrichments only:
Collection5678 <- enrichments %>% filter(collection == "collection5"|
                                           collection == "collection6"|
                                           collection == "collection7"|
                                           collection == "collection8") %>%
  filter(qValue < 0.05)

#plot of odds ratios as dot size and pvalues as heatmap color for all lists sig across samples
ggplot(Collection5678, aes(y = newnames, x = oddsRatio)) +
  facet_grid(~userSet)+
  geom_point( alpha=0.75, aes(size = userSet.in.target.list,color=qValue)) + 
  scale_size(name   = "Number of \n Overlapping Regions",
             breaks = signif(fivenum(Collection5678$userSet.in.target.list),2), #returns rounded values for 5 sets
             labels = signif(fivenum(Collection5678$userSet.in.target.list),2))+
  theme_bw() + 
  xlab("Odds Ratio") +
  ylab("Comparison List") + 
  scale_color_gradient(low="blue",high="red")+
  theme(strip.text.x = element_text(size = 14)) +
  #theme(axis.text.x=element_text(angle=65,hjust=1)) +
  theme(legend.key = element_blank(), strip.background = element_blank(),panel.border = element_rect(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

ggsave(filename="Collection5678_BMDMmarks_1_5_20.pdf",width = 16, height = 6, dpi = 300,  useDingbats = FALSE)

########################################################################################################################
#significant enrichments only: media<LPS
Collection5678 <- enrichments %>% filter(collection == "collection5"|
                                           collection == "collection6"|
                                           collection == "collection7"|
                                           collection == "collection8") %>%
  filter(qValue < 0.05) %>%
  filter(grepl("<", userSet))

#plot of odds ratios as dot size and pvalues as heatmap color for all lists sig across samples
ggplot(Collection5678, aes(y = newnames, x = oddsRatio)) +
  facet_grid(~userSet)+
  geom_point( alpha=0.75, aes(size = userSet.in.target.list,color=qValue)) + 
  scale_size(name   = "Number of \n Overlapping Regions",
             breaks = signif(fivenum(Collection5678$userSet.in.target.list),2), #returns rounded values for 5 sets
             labels = signif(fivenum(Collection5678$userSet.in.target.list),2))+
  theme_bw() + 
  xlab("Odds Ratio") +
  ylab("Comparison List") + 
  scale_color_gradient(low="red",high="blue")+
  theme(legend.key = element_blank(), strip.background = element_blank(),panel.border = element_rect(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))


ggsave(filename="LessThanPeaks_BMDMmarks_5_2_19.pdf",width = 12, height = 5, dpi = 300,  useDingbats = FALSE)


########################################################################################################################
#significant enrichments only: media>LPS
Collection5678 <- enrichments %>% filter(collection == "collection5"|
                                           collection == "collection6"|
                                           collection == "collection7"|
                                           collection == "collection8") %>%
  filter(qValue < 0.05) %>%
  filter(grepl(">", userSet))

#plot of odds ratios as dot size and pvalues as heatmap color for all lists sig across samples
ggplot(Collection5678, aes(y = newnames, x = oddsRatio)) +
  facet_grid(~userSet)+
  geom_point( alpha=0.75, aes(size = userSet.in.target.list,color=qValue)) + 
  scale_size(name   = "Number of \n Overlapping Regions",
             breaks = signif(fivenum(Collection5678$userSet.in.target.list),2), #returns rounded values for 5 sets
             labels = signif(fivenum(Collection5678$userSet.in.target.list),2))+
  theme_bw() + 
  xlab("Odds Ratio") +
  ylab("Comparison List") + 
  scale_color_gradient(low="red",high="blue")+
  theme(strip.text.x = element_text(size = 14)) +
  #theme(axis.text.x=element_text(angle=65,hjust=1)) +
  theme(legend.key = element_blank(), strip.background = element_blank(),panel.border = element_rect(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

ggsave(filename="GreaterThanPeaks_BMDMmarks_5_2_19.pdf",width = 12, height = 6, dpi = 300,  useDingbats = FALSE)









