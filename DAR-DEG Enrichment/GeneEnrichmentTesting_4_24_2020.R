#Graphs from overlap permutation analysis
#author: Annie Vogel Ciernia
#a.ciernia@gmail.com
#11/12/2018
##############################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(lsmeans)
library(nlme)
library(cowplot)
library(xlsx)
library("GeneOverlap")
library(qdapRegex)
##############################################################################################################
#repeat analysis using parsed DARs by 1 vs 2 vs 1&2
##############################################################################################################

#get gene lists associated with each peakset from ChIPseeker
setwd("/Users/aciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/geneoverlaps_4_24_2020")

#read in consensus peak genes for BG
BGgenes <- read.table("AllConsensusPeaks.ensembl.txt")

##############################################################################################################
#for Tolerized and NonTolerized Gene lists from Fluidigm Array
##############################################################################################################
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
#get gene lists from Fluidigm:
FluidigmDEG <- read.csv("MasterHeatmapdata.csv")
FluidigmDEG$gene <- as.character(FluidigmDEG$gene)
FluidigmDEG$mgi_symbol <- firstup(FluidigmDEG$gene)

#get ensembl IDs
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

DEGgenes = getBM(attributes = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position","end_position","strand"), 
                    filter= "mgi_symbol",
                    values = FluidigmDEG$mgi_symbol, 
                    mart = mouse)

mergeGene <- merge(FluidigmDEG, DEGgenes,by="mgi_symbol")
mergeGene <- distinct(mergeGene)
write.csv(mergeGene,"fluidigm_DEG.csv")


peakgenelist_FluidigmDEG <- merge(peakgenelist,mergeGene, by.x = "ENSEMBL",by.y = "ensembl_gene_id")


#simplify annotations

peakgenelist_FluidigmDEG$annotation <- rm_between_multiple(peakgenelist_FluidigmDEG$annotation, left="(EN", right=")")
peakgenelist_FluidigmDEG$annotation <- rm_between_multiple(peakgenelist_FluidigmDEG$annotation, left=" (", right=")")


##############################################################################################################
#gene counts

FluidigmCount <- peakgenelist_FluidigmDEG %>% 
  filter(annotation == "Promoter") %>%
  group_by(peaklist, C57.Catagory) %>%
  summarize(count= n())


##############################################################################################################
#DARs in promoters
#####################################################
#read in DARs
peakgenelist <- read.csv("Allgenelists.csv")
peakgenelist <- peakgenelist %>% 
  dplyr::select(V4,annotation,distanceToTSS,ENSEMBL,SYMBOL,GENENAME)
colnames(peakgenelist)[1] <- c("peaklist")


#promoters only list for genes
peakgenelist$annotation <- rm_between_multiple(peakgenelist$annotation, left=" (", right=")")
peakgenelist$annotation <- rm_between_multiple(peakgenelist$annotation, left=" (", right=")")

promoterlist <- peakgenelist %>% filter(annotation == "Promoter") %>%
  dplyr::select(peaklist,ENSEMBL) %>% distinct() 


promoterlist$ENSEMBL <- as.character(promoterlist$ENSEMBL)
promoterlist <- as.data.frame(promoterlist)

#fix names
promoterlist$peaklist <- as.character(promoterlist$peaklist)
unique(promoterlist$peaklist)

promoterlist$peaklist <- gsub("C57","C57 ",promoterlist$peaklist)
promoterlist$peaklist <- gsub("BTBR","BTBR ",promoterlist$peaklist)
promoterlist$peaklist <- gsub("LPS1","LPSx1",promoterlist$peaklist)
promoterlist$peaklist <- gsub("LPS2","LPSx2",promoterlist$peaklist)

promoterlist$peaklist <-factor(promoterlist$peaklist, levels=c(
  "BTBR media<C57 media","BTBR media>C57 media",
  "C57 media<LPSx1only","BTBR media<LPSx1only",
  "C57 media<LPSx2only","BTBR media<LPSx2only",
  "C57 media<LPSx1&2","BTBR media<LPSx1&2",
  "C57 media>LPSx1only","BTBR media>LPSx1only",
  "C57 media>LPSx2only","BTBR media>LPSx2only",
  "C57 media>LPSx1&2","BTBR media>LPSx1&2",
  "BTBR LPSx1<C57 LPSx1",
  "BTBR LPSx2<C57 LPSx2",
  "BTBR LPSx1>C57 LPSx1",
  "BTBR LPSx2>C57 LPSx2"
))


#list split
List <- split(promoterlist$ENSEMBL, promoterlist$peaklist)

#peak counts

write.csv(promoterlist, "PromoterDARs.csv")
write.csv(peakgenelist, "allDARs_subgrouped.csv")

DARpromotercounts <- table(promoterlist$peaklist)
DARsubcounts <- table(peakgenelist$peaklist)

write.csv(DARpromotercounts, "PromoterDARs_counts.csv")
write.csv(DARsubcounts, "allDARs_subgrouped_counts.csv")


##############################################################################################################
#BTBR SNPs and INdels
#####################################################
#coding changes:
coding <- read.csv("/Users/aciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/BTBRgenome/BTBRsnps_10_9_18/BTBR_CodingVariants.csv")
coding$listname <- paste(coding$Type, coding$variant, sep=" ")
codingvar <- coding %>% dplyr::select(listname,ENSEMBL)

#SNPs and indels in any region
#dfvcfannoindel_genes
#dfvcfannosnp_genes
load("/Users/aciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/BTBRgenome/BTBRsnps_10_9_18/snp_data/BTBR_SNP_annotation.RData")
load("/Users/aciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/BTBRgenome/BTBRsnps_10_9_18/snp_data/BTBR_indel_annotation.RData")

#filter for promoters
#pull out promoter SNPs and Indels
SNPpromoters <- dfvcfannosnp_genes %>% filter(LOCATION =="promoter") %>% 
  mutate(listname =  c("Promoter SNP")) %>%
  dplyr::select(listname,ENSEMBL) %>% na.omit() %>% distinct() #9231

IndelPpromoters <- dfvcfannoindel_genes %>% filter(LOCATION =="promoter") %>% 
  mutate(listname =  c("Promoter InDel")) %>%
  dplyr::select(listname,ENSEMBL) %>% na.omit() %>%
  distinct() #8719

#Combine
BTBRvar <- rbind(codingvar,IndelPpromoters, SNPpromoters)
counts <- table(BTBRvar$listname)
write.csv(counts,"CountsBTBRgenevariants.csv")
write.csv(BTBRvar , "BTBRvariants_genes.csv")

VariantList <- split(BTBRvar$ENSEMBL, BTBRvar$listname)


##############################################################################################################
#DEG list splits
##############################################################################################################

#fluidigm
FluidigmDEG_list <- mergeGene %>% 
  filter(C57.Catagory != 'C57 not induced LPS1') %>%
  filter(C57.Catagory != 'C57 sensitized')  %>%
  dplyr::select(ensembl_gene_id,mediaDifferences,C57.Catagory,BTBR.Catagory,HyperLPSx1,HyperLPSx2) %>% distinct()


FluidigmDEG_list <- as.data.frame(FluidigmDEG_list)


  
FluidigmDEG_list1 <- FluidigmDEG_list[,c(1,2)]
colnames(FluidigmDEG_list1) <- c("ensembl_gene_id","catagory")
FluidigmDEG_list2 <- FluidigmDEG_list[,c(1,3)]
colnames(FluidigmDEG_list2) <- c("ensembl_gene_id","catagory")
FluidigmDEG_list3 <- FluidigmDEG_list[,c(1,4)]
colnames(FluidigmDEG_list3) <- c("ensembl_gene_id","catagory")
FluidigmDEG_list4 <- FluidigmDEG_list[,c(1,5)]
colnames(FluidigmDEG_list4) <- c("ensembl_gene_id","catagory")
FluidigmDEG_list5 <- FluidigmDEG_list[,c(1,6)]
colnames(FluidigmDEG_list5) <- c("ensembl_gene_id","catagory")

FluidigmDEG_list <- rbind(
  FluidigmDEG_list1,
  FluidigmDEG_list2,
  FluidigmDEG_list3,
  FluidigmDEG_list4,
  FluidigmDEG_list5)

FluidigmDEG_list$ensembl_gene_id <- as.character(FluidigmDEG_list$ensembl_gene_id)
FluidigmDEG_list$catagory <- as.character(FluidigmDEG_list$catagory )
unique(FluidigmDEG_list$catagory)

FluidigmDEG_list <- FluidigmDEG_list %>% 
filter(catagory != 'BTBR sensitized')  %>%
  filter(catagory != 'Normal')

FluidigmDEG_list$catagory <- factor(FluidigmDEG_list$catagory, levels= c(
  "no baseline differences","baseline differences","C57 tolerized","BTBR tolerized","C57 nontolerized","BTBR nontolerized",
  "BTBR Hyper LPSx1"     ,   "BTBR Hyper LPSx2" ) )

#split list
FluidigmDEG_list <- split(FluidigmDEG_list$ensembl_gene_id, FluidigmDEG_list$catagory)

BG <- length(BGgenes$V1)


##############################################################################################################
#overlaps
##############################################################################################################

gom.obj <- newGOM(DEGs_Var,FluidigmDEG_list,genome.size=BG)

go.obj.test <-  getMatrix(gom.obj, name="pval")

fdr.test <- matrix(p.adjust(as.vector(as.matrix(go.obj.test)), method='fdr'),ncol=ncol(go.obj.test))
rownames(fdr.test) <- rownames(go.obj.test)
colnames(fdr.test) <- colnames(go.obj.test)
fdr.test <- as.matrix(fdr.test)
write.csv(fdr.test,"FDR_PeakGenes_overlap_TolNonTolDEG.csv")

OR <- getMatrix(gom.obj, "odds.ratio")
write.csv(OR,"OR_PeakGenes_overlap_TolNonTolDEG.csv")

inter.nl <- getNestedList(gom.obj, name="intersection")

inter.count <- do.call(rbind, lapply(inter.nl, lengths))
inter.count <- t(inter.count)
write.csv(inter.count,"OverlapCount_PeakGenes_overlap_TolNonTolDEG.csv")


#plot colored by pvalue and labeled with overlap #
#set colors
col3 <- colorRampPalette(c("red","white","blue")) 

#logOverlap <-log(OverlapPercent+1)
neglogFDR <- log(fdr.test)*-1

fdr.test2 <- fdr.test
fdr.test2[fdr.test2 >0.05 ] <- 1


pdf("DEG_FishersExactEnrichmentPeakPromoterGenes.pdf", width=12, height=8)

corrplot(neglogFDR,p.mat = fdr.test2,method ="color", 
         title="Overlap Enrichments DEG and Promoter Peaks", 
         tl.col="black", 
         tl.cex=1.5, is.corr = FALSE,
         #diag=FALSE, 
         addrect=1, 
         mar=c(0,0,2,2), 
         rect.col = "black", 
         addgrid.col = "darkgray",
         cl.lim=c(0, 20),
         #tl.srt = 60,
        # cl.lim = c(0, 25), 
         #type = "upper",
         insig = "label_sig", 
        col = col3(60),
         pch.col = "black",pch.cex = 3
         )

invisible(dev.off())



##############################################################################################################
#overlaps w/genetic variants
##############################################################################################################
#combine lists
List2 <- append(List, FluidigmDEG_list)

gom.obj <- newGOM(FluidigmDEG_list,VariantList,genome.size=BG)

go.obj.test <-  getMatrix(gom.obj, name="pval")

fdr.test <- matrix(p.adjust(as.vector(as.matrix(go.obj.test)), method='fdr'),ncol=ncol(go.obj.test))
rownames(fdr.test) <- rownames(go.obj.test)
colnames(fdr.test) <- colnames(go.obj.test)
fdr.test <- as.matrix(fdr.test)
write.csv(fdr.test,"FDR_PeakGenes_overlap_GeneticVariants.csv")

OR <- getMatrix(gom.obj, "odds.ratio")
write.csv(OR,"OR_PeakGenes_overlap_GeneticVariants.csv")

inter.nl <- getNestedList(gom.obj, name="intersection")

inter.count <- do.call(rbind, lapply(inter.nl, lengths))
inter.count <- t(inter.count)
write.csv(inter.count,"OverlapCount_PeakGenes_overlap_GeneticVariants.csv")


#plot colored by pvalue and labeled with overlap #
#set colors
col3 <- colorRampPalette(c("red","white","blue")) 

#logOverlap <-log(OverlapPercent+1)
neglogFDR <- log(fdr.test)*-1

fdr.test2 <- fdr.test
fdr.test2[fdr.test2 >0.05 ] <- 1


pdf("DEG_FishersExactEnrichmentPeakPromoterGenes.pdf", width=12, height=8)

corrplot(neglogFDR,p.mat = fdr.test2,method ="color", 
         title="Overlap Enrichments DEG and Promoter Peaks", 
         tl.col="black", 
         tl.cex=1.5, is.corr = FALSE,
         #diag=FALSE, 
         addrect=1, 
         mar=c(0,0,2,2), 
         rect.col = "black", 
         addgrid.col = "darkgray",
         cl.lim=c(0, 20),
         #tl.srt = 60,
         # cl.lim = c(0, 25), 
         #type = "upper",
         insig = "label_sig", 
         col = col3(60),
         pch.col = "black",pch.cex = 3
)

invisible(dev.off())


