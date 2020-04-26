##########################################################################
#AVC a.ciernia@gmail.com
#8/2018
#ATACseq BMDM +/- 1 or 2 doses LPS in vitro
##########################################################################
#Annotation of peaks with ChIPseeker and enrichment testing
##########################################################################
library(tools)
library(ChIPseeker)
library(rlist)
##########################################################################

#mouse:
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
#biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
annoDb <- org.Mm.eg.db
#######################################################################################################################
# Annotation -------------------------------------------------------------
#######################################################################################################################
setwd("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/peakfiles_5_5_19")

#get lists
files <- list.files(path=".", pattern="*.mm10.bed", all.files=T, full.names=T)
filelist <- lapply(files, readPeakFile)

#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files)))
names_list <- gsub("_DEpeaks.mm10", "", names_list)
names_list <- gsub(".mm10", "", names_list)
names_list <- gsub("_", " ", names_list)
names_list <- gsub(" ", "", names_list)

names(filelist) <- names_list
names(filelist)


# Region overlaps ---------------------------------------------------------
#BiocManager::install("ChIPpeakAnno", version = "3.8")
library(ChIPpeakAnno)

# peaklistnames <- names(filelist[c(2:14)])
# 
# BTBRmedialessC57media <- filelist$`BTBRmedia<C57media`
# 
# BTBRmediagreaterC57media <- filelist$`BTBRmedia>C57media`
# 
# for (i in peaklistnames) {
#   
#   ol <- findOverlapsOfPeaks(BTBRmedialessC57media,BTBRmediagreaterC57media,filelist[[i]])
#   
#   pdf(file = paste(i,"_Venn.pdf",sep=""), wi = 5, he = 5, useDingbats = F)
#   makeVennDiagram(ol)
#   dev.off()      
#   
# }



# Annotation with genes ---------------------------------------------------


#if overlap="all", then gene overlap with peak will be reported as nearest gene, 
#no matter the overlap is at TSS region or not.

peakAnnoList <- lapply(filelist, annotatePeak, TxDb=txdb, annoDb = "org.Mm.eg.db", overlap = "all")

#annotate background > all called peaks from all replicates
#
BG <- readPeakFile("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/DiffBind/Allconsensuspeaks.bed")

BGanno <- annotatePeak(BG, TxDb=txdb, annoDb = "org.Mm.eg.db", overlap = "all")

#add background to list
peakAnnoList2 <- list.append(peakAnnoList,BGanno)
names(peakAnnoList2)[19] <-c("AllConsensusPeaks")
names(peakAnnoList2) 

#write to excel file
setwd("/Users/annieciernia/Box\ Sync/Ashwood_labcollaboration/Experiments/atac_2019/genes_pathways")
save(peakAnnoList2,file="peakAnnoList2.Rdata")
#load("peakAnnoList2.Rdata")
names(peakAnnoList2) 

#make output file
for (i in 1:length(peakAnnoList2)) {
  print(i)
  
  name <- names(peakAnnoList2[i])
  tmp <- peakAnnoList2[[i]]
  
  tmpdf <- as.data.frame(tmp)
  openxlsx::write.xlsx(tmpdf, file = paste("Annotation_ChIPseeker",name,"peaks.xlsx",sep=""), sep= "")
}
openxlsx::write.xlsx(tmpdf, file = paste("Annotation_ChIPseeker",name,"peaks.xlsx",sep=""), sep= "")

# make output file -------------------------------------------

#make output file > gene names
Genelist <- NULL
lengthout <- NULL
for (i in 1:(length(peakAnnoList2)-1)) {
  print(i)
  
  name <- as.character(names(peakAnnoList2[i]))
  tmp <- peakAnnoList2[[i]]
  
  tmpdf <- as.data.frame(tmp)
  tmpdf$V4 <- name
  
  Genelist <- rbind(tmpdf,Genelist)
  
  tmpdf2 <- tmpdf %>% dplyr::select(ENSEMBL) %>%
    unique()
  #write out gene lists: 
  #write.table(tmpdf2,file = paste(name,".ensembl.txt",sep=""),quote=F, col.names = F,row.names = F,sep="\t" )
  
  #save list length to file
  listlength <- length(tmpdf2$ENSEMBL)
  lengthdf <- cbind(name,listlength)
  lengthout <- rbind(lengthdf,lengthout)
}

write.csv(Genelist,"Allgenelists.csv")

#write out BG regions
tmp <- peakAnnoList2[[19]]
tmpdf <- as.data.frame(tmp)
name <- c("AllConsensusPeaks")
#write out gene lists: 
tmpdf2 <- tmpdf %>% dplyr::select(ENSEMBL) %>%
  unique()
#write out gene lists: 
write.table(tmpdf2,file = paste(name,".ensembl.txt",sep=""),quote=F, col.names = F,row.names = F,sep="\t" )

#length
#save list length to file
listlength <- length(tmpdf2$ENSEMBL)
lengthdf <- cbind(name,listlength)
lengthout <- rbind(lengthdf,lengthout)

library(xlsx)
write.xlsx(lengthout,file="GeneCounts_perPeaklist.xlsx")


# Plots genomic region enrichment -------------------------------------------
names(peakAnnoList2)

#media baseline differences between genotypes
p_all <- plotDistToTSS(peakAnnoList2[c(19,5,6)],title="Distribution of transcription factor-binding loci\nrelative to TSS")
pdf("MediaDiff_Annotation_plotDisToTSS.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

p_all <- plotAnnoBar(peakAnnoList2[c(19,5,6)],title="Distribution of Regions in Genomic Elements")
pdf("MediaDiff_Annotation_barplot.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

#C57 LPS enrichments
p_all <- plotDistToTSS(peakAnnoList2[c(19,14,15,13,17,18,16)],title="Distribution of transcription factor-binding loci\nrelative to TSS")
pdf("C57LPS_Annotation_plotDisToTSS.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

p_all <- plotAnnoBar(peakAnnoList2[c(19,14,15,13,17,18,16)],title="Distribution of Regions in Genomic Elements")
pdf("C57LPS_Annotation_barplot.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

#BTBR LPS enrichments
p_all <- plotDistToTSS(peakAnnoList2[c(19,8,9,7,11,12,10)],title="Distribution of transcription factor-binding loci\nrelative to TSS")
pdf("BTBRLPS_Annotation_plotDisToTSS.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

p_all <- plotAnnoBar(peakAnnoList2[c(19,8,9,7,11,12,10)],title="Distribution of Regions in Genomic Elements")
pdf("BTBRLPS_Annotation_barplot.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

#BTBR vs C57 LPS enrichments
p_all <- plotDistToTSS(peakAnnoList2[c(19,1:4)],title="Distribution of transcription factor-binding loci\nrelative to TSS")
pdf("BTBRvsC57LPS_Annotation_plotDisToTSS.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

p_all <- plotAnnoBar(peakAnnoList2[c(19,1:4)],title="Distribution of Regions in Genomic Elements")
pdf("BTBRvsC57LPS_Annotation_barplot.pdf", height = 6, width =6,useDingbats=FALSE)    # create PNG for the heat map       
print(p_all)
dev.off()

########################################################################################################################
# Fishers exact test genomic region enrichment -------------------------------------------
########################################################################################################################
#Make counts table from annotation file outputs

#function to count the number of occurances of peaks in each genomic feature 
#from annotation output, calculates total peaks and percentage (count/total *100)
library(qdapRegex)
library(dplyr)
sum_count <- function(excelfilename){
  
  #open dataframe
  datframe <- openxlsx::read.xlsx(excelfilename)
  
  #remove intron and exon designations
  datframe$annotation <- rm_between_multiple(datframe$annotation, left="(uc", right=")")
  
  #make df
  count_df <- datframe %>% 
    dplyr::group_by(annotation) %>%
    dplyr::summarise(count = n()) %>%
    mutate(total = sum(count)) %>%
    mutate(percent = count/total*100)
  
  count_df <- as.data.frame(count_df)
  count_df
}

#sum_count(excelfilename = files[2])

files2 <- list.files(pattern = "^Annotation_ChIPseeker")

#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files2)))
names_list <- gsub("peaks", "", names_list)
names_list <- gsub("Annotation_ChIPseeker", "", names_list)

names(files2) <- names_list
names(files2)

#clean files for write to single excel doc
cleanfiles <- function(excelfilename){
  
  #open dataframe
  datframe <- openxlsx::read.xlsx(excelfilename)
  
  #remove intron and exon designations
  datframe$annotation <- rm_between_multiple(datframe$annotation, left="(uc", right=")")
  
  datframe <- as.data.frame(datframe)
  datframe
}

filewrite <- lapply(files2, cleanfiles)
names(filewrite)

#write to one file:
library(writexl) #if gives DLL error, restart R with fewer packages loaded
write_xlsx(filewrite, path="AllFiles_Annotation_ChIPseeker.xlsx")


#get count data
sumcountList <- lapply(files2, sum_count)

#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files2)))
names_list <- gsub("peaks", "", names_list)
names_list <- gsub("Annotation_ChIPseeker", "", names_list)

names(sumcountList) <- names_list
names(sumcountList) 

#two sided fisher's exact test
call <- names(sumcountList)
masterfile <- NULL
for (i in call) {
  print(i)
  
  #get background data
  BG <- sumcountList[]$`AllConsensusPeaks`
  BG <- as.data.frame(BG)
  
  #all annotations possible
  annotations <- unique(BG$annotation)
  
  #get dataset
  data <- sumcountList[[i]]
  data <- as.data.frame(data)
  
  data2 <- merge(BG,data,by ="annotation",all.x=T)
  
  #if missing annotation replace with zero
  data2[is.na(data2)] <- 0
  
  #for each annotation:
  out <- NULL
  for(an in annotations) {
    
    tmp <- filter(data2,annotation == an)
    #for DEpeak list:
    roiplus <- tmp$count.y
    roineg <- tmp$total.y - tmp$count.y
    #for BG peak list:
    BGplus <- tmp$count.x
    BGneg <- tmp$total.x - tmp$count.x
    mytable <- rbind(c(roiplus,BGplus),c(roineg,BGneg))
    test <- fisher.test(mytable,alternative='greater')
    
    pvalue <- as.numeric(as.character(test$p.value))
    HigherCI <- as.numeric(as.character(test$conf.int[1]))
    LowerCI <- as.numeric(as.character(test$conf.int[2]))
    OddsRatio <- as.numeric(as.character(test$estimate))
    
    results <- data.frame(annotation=paste(an),pvalue,HigherCI,LowerCI,OddsRatio)
    out <- rbind(out,results) 
  }
  
  #make masteroutfile    
  out <- as.data.frame(out)
  out$comparison <- paste(names(sumcountList[i]))
  
  masterfile <- rbind(out,masterfile)  
  
}

masterfile <- filter(masterfile, comparison != "AllConsensusPeaks")

masterfile$FDR <- p.adjust(masterfile$pvalue, method="fdr") #FDR

openxlsx::write.xlsx(masterfile, file = "OneTailFishersTest_Annotations_ChIPseeker.xlsx")



# gene ids from annotation ---------------------------------------
# 
genes = lapply(peakAnnoList2, function(i) as.data.frame(i)$geneId)
# 
# #get lists of ensembl and symbols
# ensembl <- lapply(genes, function(i) as.data.frame(i)$ENSEMBL)
# symbol <- lapply(genes, function(i) as.data.frame(i)$SYMBOL)
# 
# #remove duplicates:
# ensembl <-lapply(ensembl, function(x) x[!duplicated(x[])])
# symbol <- lapply(symbol, function(x) x[!duplicated(x[])])

# #write to files
# for(i in 1:length(ensembl)){
#   write.table(data.frame(ensembl[[i]]),file=paste(names(ensembl[i]),'ensemble.txt', sep="."),append=F, row.names=FALSE,col.names=F,quote=F)
# }
# 
# for(i in 1:length(symbol)){
#   write.table(data.frame(symbol[[i]]),file=paste(names(symbol[i]),'symbol.txt', sep="."),append=F, row.names=FALSE,col.names=F,quote=F)
# }

#subset lists
#all conditions:

#media differences
medialist <- genes[c(5,6)]
names(medialist)

#C57 LPS enrichments
C57LPSlist <- genes[c(13:18)]
names(C57LPSlist)

#BTBR LPS enrichments
BTBRLPSlist <- genes[c(7:12)]
names(BTBRLPSlist)

#allconsensuspeaks = univers
BG <- genes[c(19)]
names(BG)

masterlist <- list(medialist,C57LPSlist,BTBRLPSlist)
names(masterlist) <- c("MediaDiff","C57_LPS","BTBR_LPS")

# Functional and Pathway annotation ---------------------------------------
library(clusterProfiler)
library(ggplot2)
#http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#heatmap-of-chip-binding-to-tss-regions
#http://bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

# loop through enrichements for each list --------------------------------------------------------------------

for (i in names(masterlist)) {
  
  print(i)
  tmplist <- masterlist[[i]]
  names(tmplist)
  
#KEGG enrichments 
compKEGG <- compareCluster(geneCluster   = tmplist,
                            fun           = "enrichKEGG",
                            organism = "mmu", #see http://www.genome.jp/kegg/catalog/org_list.html
                            pvalueCutoff  = 0.05,
                            universe = BG$AllConsensusPeaks,#all consesnsus peaks as backgroudn
                            pAdjustMethod = "BH")

p <- dotplot(compKEGG, showCategory = 10, title = paste("KEGG Pathway Enrichment Analysis ",i,sep=""))
p + theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(filename =paste(i,"_KEGG.pdf",sep=""), height = 10, width =20,useDingbats=FALSE) 
#pdf(file =paste(i,"_KEGG.pdf",sep=""), height = 10, width =20,useDingbats=FALSE)    # create PNG for the heat map       


#dev.off();

compKEGG_list <- as.data.frame(compKEGG)
openxlsx::write.xlsx(compKEGG_list, file = paste(i,"_KEGG_enrichements.xlsx",sep=""))

# #source("https://bioconductor.org/biocLite.R")
# #biocLite("pathview")
# library("pathview")
# 
# #IL-6
# mmu04657 <- pathview(gene.data  = list1[[3]],
#                      pathway.id = "mmu04657",
#                      species    = "mmu")

# GO over-representation test ---------------------------------------------

#GO CC
compGO_CC <- compareCluster(geneCluster   = tmplist,
                            fun           = "enrichGO",
                            OrgDb    = org.Mm.eg.db,
                            ont      = "CC",
                            # level    = 3,
                            keyType       = 'ENTREZID',
                            readable = TRUE,
                            pvalueCutoff  = 0.05,
                            universe = BG$AllConsensusPeaks,#all consensus peaks
                            qvalueCutoff  = 0.05
)

compGO_CC_list1 <- as.data.frame(compGO_CC)
openxlsx::write.xlsx(compGO_CC_list1, file= paste(i,"_CCenrichment.xlsx",sep=""))

cc2 <- simplify(compGO_CC, cutoff=0.7, by="p.adjust", select_fun=min)
cc <- as.data.frame(cc2)
openxlsx::write.xlsx(cc2, file=paste(i,"_TrimmedCCenrichment.xlsx",sep=""))

p <- dotplot(cc2, showCategory = 10, title = paste("Cellular Component Gene Ontology Enrichment ",i,sep=""))
p + theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(filename =paste(i,"_TrimCC.pdf",sep=""), height = 10, width =20,useDingbats=FALSE)    # create PNG for the heat map       


#GO bp
compGO_bp <- compareCluster(geneCluster   = tmplist,
                            fun           = "enrichGO",
                            OrgDb    = org.Mm.eg.db,
                            ont      = "BP",
                            # level    = 3,
                            keyType       = 'ENTREZID',
                            readable = TRUE,
                            pvalueCutoff  = 0.05,
                            universe = BG$AllConsensusPeaks,#all consensus peaks
                            qvalueCutoff  = 0.05
)

compGO_bp_list1 <- as.data.frame(compGO_bp)
openxlsx::write.xlsx(compGO_bp_list1, file= paste(i,"_BPenrichment.xlsx",sep=""))

bp2 <- simplify(compGO_bp, cutoff=0.7, by="p.adjust", select_fun=min)
bp <- as.data.frame(bp2)
openxlsx::write.xlsx(bp2, file= paste(i,"_TrimmedBPenrichment.xlsx",sep=""))

p <- dotplot(bp2, showCategory = 10, title = paste("Biological Process Gene Ontology Enrichment ",i,sep=""))
p + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(filename=paste(i,"_TrimBP.pdf",sep=""), height = 10, width =14,useDingbats=FALSE)    # create PNG for the heat map       


#GO mf
compGO_mf <- compareCluster(geneCluster   = tmplist,
                            fun           = "enrichGO",
                            OrgDb    = org.Mm.eg.db,
                            ont      = "MF",
                            # level    = 3,
                            keyType       = 'ENTREZID',
                            readable = TRUE,
                            pvalueCutoff  = 0.05,
                            universe = BG$AllConsensusPeaks,#all consensus peaks
                            qvalueCutoff  = 0.05
)

compGO_mf_list1 <- as.data.frame(compGO_mf)
openxlsx::write.xlsx(compGO_mf_list1, file= paste(i,"_MFenrichment.xlsx",sep=""))

mf2 <- simplify(compGO_mf, cutoff=0.7, by="p.adjust", select_fun=min)
mf <- as.data.frame(mf2)
openxlsx::write.xlsx(mf2, file= paste(i,"_TrimmedMFenrichment.xlsx",sep=""))

p <- dotplot(mf2, showCategory = 10, title = paste("Molecular Function Gene Ontology Enrichment ",i,sep=""))
p + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(filename=paste(i,"_TrimMF.pdf",sep=""), height = 10, width =16,useDingbats=FALSE)    # create PNG for the heat map       



}


# enrichments for media<LPS and media>LPS for both geneotypes--------------------------------------
#media<LPS
media_less_LPS <- genes[c(14,8,15,9,13,7)]
names(media_less_LPS)

#media>LPS
media_greater_LPS <- genes[c(17,11,18,12,16,10)]
names(media_greater_LPS)

#allconsensuspeaks = univers
BG <- genes[c(19)]
names(BG)

masterlist <- list(media_less_LPS,media_greater_LPS)
names(masterlist) <- c("media<LPS","media>LPS")

# Functional and Pathway annotation ---------------------------------------
library(clusterProfiler)
library(ggplot2)
#http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#heatmap-of-chip-binding-to-tss-regions
#http://bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

# loop through enrichements for each list --------------------------------------------------------------------

for (i in names(masterlist)) {
  
  print(i)
  tmplist <- masterlist[[i]]
  names(tmplist)
  
  #KEGG enrichments 
  compKEGG <- compareCluster(geneCluster   = tmplist,
                             fun           = "enrichKEGG",
                             organism = "mmu", #see http://www.genome.jp/kegg/catalog/org_list.html
                             pvalueCutoff  = 0.05,
                             universe = BG$AllConsensusPeaks,#all consesnsus peaks as backgroudn
                             pAdjustMethod = "BH")
  
  p <- dotplot(compKEGG, showCategory = 10, title = paste("KEGG Pathway Enrichment Analysis ",i,sep=""))
  p + theme(axis.text.x=element_text(angle=45, hjust=1))
  
  ggsave(filename =paste(i,"_KEGG.pdf",sep=""), height = 10, width =14,useDingbats=FALSE) 
  #pdf(file =paste(i,"_KEGG.pdf",sep=""), height = 10, width =20,useDingbats=FALSE)    # create PNG for the heat map       
  
  
  #dev.off();
  
  compKEGG_list <- as.data.frame(compKEGG)
  openxlsx::write.xlsx(compKEGG_list, file = paste(i,"_KEGG_enrichements.xlsx",sep=""))
  
  # #source("https://bioconductor.org/biocLite.R")
  # #biocLite("pathview")
  # library("pathview")
  # 
  # #IL-6
  # mmu04657 <- pathview(gene.data  = list1[[3]],
  #                      pathway.id = "mmu04657",
  #                      species    = "mmu")
  
  # GO over-representation test ---------------------------------------------
  
  #GO CC
  compGO_CC <- compareCluster(geneCluster   = tmplist,
                              fun           = "enrichGO",
                              OrgDb    = org.Mm.eg.db,
                              ont      = "CC",
                              # level    = 3,
                              keyType       = 'ENTREZID',
                              readable = TRUE,
                              pvalueCutoff  = 0.05,
                              universe = BG$AllConsensusPeaks,#all consensus peaks
                              qvalueCutoff  = 0.05
  )
  
  compGO_CC_list1 <- as.data.frame(compGO_CC)
  openxlsx::write.xlsx(compGO_CC_list1, file= paste(i,"_CCenrichment.xlsx",sep=""))
  
  cc2 <- simplify(compGO_CC, cutoff=0.7, by="p.adjust", select_fun=min)
  cc <- as.data.frame(cc2)
  openxlsx::write.xlsx(cc2, file=paste(i,"_TrimmedCCenrichment.xlsx",sep=""))
  
  p <- dotplot(cc2, showCategory = 10, title = paste("Cellular Component Gene Ontology Enrichment ",i,sep=""))
  p + theme(axis.text.x=element_text(angle=45, hjust=1))
  
  ggsave(filename =paste(i,"_TrimCC.pdf",sep=""), height = 10, width =14,useDingbats=FALSE)    # create PNG for the heat map       
  
  
  #GO bp
  compGO_bp <- compareCluster(geneCluster   = tmplist,
                              fun           = "enrichGO",
                              OrgDb    = org.Mm.eg.db,
                              ont      = "BP",
                              # level    = 3,
                              keyType       = 'ENTREZID',
                              readable = TRUE,
                              pvalueCutoff  = 0.05,
                              universe = BG$AllConsensusPeaks,#all consensus peaks
                              qvalueCutoff  = 0.05
  )
  
  compGO_bp_list1 <- as.data.frame(compGO_bp)
  openxlsx::write.xlsx(compGO_bp_list1, file= paste(i,"_BPenrichment.xlsx",sep=""))
  
  bp2 <- simplify(compGO_bp, cutoff=0.7, by="p.adjust", select_fun=min)
  bp <- as.data.frame(bp2)
  openxlsx::write.xlsx(bp2, file= paste(i,"_TrimmedBPenrichment.xlsx",sep=""))
  
  p <- dotplot(bp2, showCategory = 10, title = paste("Biological Process Gene Ontology Enrichment ",i,sep=""))
  p + theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(filename=paste(i,"_TrimBP.pdf",sep=""), height = 10, width =14,useDingbats=FALSE)    # create PNG for the heat map       
  
  
  #GO mf
  compGO_mf <- compareCluster(geneCluster   = tmplist,
                              fun           = "enrichGO",
                              OrgDb    = org.Mm.eg.db,
                              ont      = "MF",
                              # level    = 3,
                              keyType       = 'ENTREZID',
                              readable = TRUE,
                              pvalueCutoff  = 0.05,
                              universe = BG$AllConsensusPeaks,#all consensus peaks
                              qvalueCutoff  = 0.05
  )
  
  compGO_mf_list1 <- as.data.frame(compGO_mf)
  openxlsx::write.xlsx(compGO_mf_list1, file= paste(i,"_MFenrichment.xlsx",sep=""))
  
  mf2 <- simplify(compGO_mf, cutoff=0.7, by="p.adjust", select_fun=min)
  mf <- as.data.frame(mf2)
  openxlsx::write.xlsx(mf2, file= paste(i,"_TrimmedMFenrichment.xlsx",sep=""))
  
  p <- dotplot(mf2, showCategory = 10, title = paste("Molecular Function Gene Ontology Enrichment ",i,sep=""))
  p + theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(filename=paste(i,"_TrimMF.pdf",sep=""), height = 10, width =14,useDingbats=FALSE)    # create PNG for the heat map       
  
  
  
}

# gene overlap upset plots ----------------------------------------------------

library(UpSetR)

save(genes,filelist,file="BTBR_C57regionslist.RData")

genes_upset <- genes[c(5,6,14,15,13,8,9,7,17,18,16,11,12,10,1:4)]
names(genes_upset)

pdf("UpSetplot_GenesOverlap.pdf", width=20, height=20,useDingbats = F)

upset(fromList(genes_upset), sets = rev(names(genes_upset)),
      nsets=length(genes_upset), nintersects=50,  keep.order = TRUE,
      show.numbers="yes", 
      main.bar.color="#ea5d4e", sets.bar.color="#317eab", 
      empty.intersections=NULL,order.by ="freq",
      number.angles = 0, mainbar.y.label ="No. of Intersections", 
      sets.x.label ="Set size",text.scale=2)

dev.off()

