##########################################################################
#AVC a.ciernia@gmail.com
#8/2018
#ATACseq BMDM +/- 1 or 2 doses LPS in vitro
##########################################################################
#Analysis of overlapping gene lists from DiffBind analysis output
##########################################################################
library(tools)
library(valr)
library(regioneR)
library(ChIPseeker)
##########################################################################

#make GRanges object for each bed and save to list
setwd("/Users/annieciernia/Box Sync/Ashwood_labcollaboration/Experiments/atac_2019/DiffBind")

path <- getwd()

#https://stackoverflow.com/questions/27911604/loop-through-all-files-in-a-directory-read-and-save-them-in-an-object-r

#get lists
files <- list.files(path=".", pattern="*.mm10.bed", all.files=T, full.names=T)
filelist <- lapply(files, read_bed,n_fields =4)

#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files)))
names_list <- gsub("_DEpeaks.mm10", "", names_list)
names_list <- gsub(".mm10", "", names_list)
names_list <- gsub(" ", "", names_list)
names_list <- gsub("_", " ", names_list)

names(filelist) <- names_list
names(filelist)

# Define main groups with Valr---------------------------------------

#C57 media vs LPS1 shared with C57 media vs LPS2 => genes commonly upregulated by LPS
C57common_up <- valr::bed_intersect(x= filelist[["C57LPS1 > C57media"]],
                                 y= filelist[["C57LPS2 > C57media"]])

#regions open with only 1 LPS but not 2
C57LPSx1only_up <- valr::bed_intersect(x= filelist[["C57LPS1 > C57media"]],
                                    y= filelist[["C57LPS2 > C57media"]],invert=T)

#regions open with only 2 LPS but not 1
C57LPSx2only_up <- valr::bed_intersect(x= filelist[["C57LPS2 > C57media"]],
                                    y= filelist[["C57LPS1 > C57media"]],invert=T)


#C57 media vs LPS1 shared with C57 media vs LPS2 => genes commonly upregulated by LPS
C57common_down <- valr::bed_intersect(x= filelist[["C57LPS1 < C57media"]],
                                 y= filelist[["C57LPS2 < C57media"]])

#regions closed with only 1 LPS but not 2
C57LPSx1only_down <- valr::bed_intersect(x= filelist[["C57LPS1 < C57media"]],
                                    y= filelist[["C57LPS2 < C57media"]],invert=T)

#regions closed with only 1 LPS but not 2
C57LPSx2only_down <- valr::bed_intersect(x= filelist[["C57LPS2 < C57media"]],
                                         y= filelist[["C57LPS1 < C57media"]],invert=T)


#validate <- valr::bed_intersect(x= C57LPSx2only, y= C57LPSx1only,invert=T)

#BTBR media vs LPS1 shared with BTBR media vs LPS2 => genes commonly upregulated by LPS
BTBRcommon_up <- valr::bed_intersect(x= filelist[["BTBRLPS1 > BTBRmedia"]],
                                     y= filelist[["BTBRLPS2 > BTBRmedia"]])

#regions open with only 1 LPS but not 2
BTBRLPSx1only_up <- valr::bed_intersect(x= filelist[["BTBRLPS1 > BTBRmedia"]],
                                        y= filelist[["BTBRLPS2 > BTBRmedia"]],invert=T)

#regions open with only 2 LPS but not 1
BTBRLPSx2only_up <- valr::bed_intersect(x= filelist[["BTBRLPS2 > BTBRmedia"]],
                                        y= filelist[["BTBRLPS1 > BTBRmedia"]],invert=T)


#BTBR media vs LPS1 shared with BTBR media vs LPS2 => genes commonly upregulated by LPS
BTBRcommon_down <- valr::bed_intersect(x= filelist[["BTBRLPS1 < BTBRmedia"]],
                                       y= filelist[["BTBRLPS2 < BTBRmedia"]])

#regions closed with only 1 LPS but not 2
BTBRLPSx1only_down <- valr::bed_intersect(x= filelist[["BTBRLPS1 < BTBRmedia"]],
                                          y= filelist[["BTBRLPS2 < BTBRmedia"]],invert=T)

#regions closed with only 1 LPS but not 2
BTBRLPSx2only_down <- valr::bed_intersect(x= filelist[["BTBRLPS2 < BTBRmedia"]],
                                          y= filelist[["BTBRLPS1 < BTBRmedia"]],invert=T)

# Define sub groups with Valr---------------------------------------

#make 3 col bed files up
C57common_up_bed <- cbind(C57common_up[,1:3],c("C57media<LPS1&2"))
colnames(C57common_up_bed) <- c("chrom","start","end","condition")

BTBRcommon_up_bed <- cbind(BTBRcommon_up[,1:3],c("BTBRmedia<LPS1&2"))
colnames(BTBRcommon_up_bed) <- c("chrom","start","end","condition")

C57LPSx1only_up_bed <- cbind(C57LPSx1only_up[,1:3],c("C57media<LPS1only"))
colnames(C57LPSx1only_up_bed) <- c("chrom","start","end","condition")

BTBRLPSx1only_up_bed <- cbind(BTBRLPSx1only_up[,1:3],c("BTBRmedia<LPS1only"))
colnames(BTBRLPSx1only_up_bed) <- c("chrom","start","end","condition")

C57LPSx2only_up_bed <- cbind(C57LPSx2only_up[,1:3],c("C57media<LPS2only"))
colnames(C57LPSx2only_up_bed) <- c("chrom","start","end","condition")

BTBRLPSx2only_up_bed <- cbind(BTBRLPSx2only_up[,1:3],c("BTBRmedia<LPS2only"))
colnames(BTBRLPSx2only_up_bed) <- c("chrom","start","end","condition")

#make 3 col bed files down
C57common_down_bed <- cbind(C57common_down[,1:3],c("C57media>LPS1&2"))
colnames(C57common_down_bed) <- c("chrom","start","end","condition")

BTBRcommon_down_bed <- cbind(BTBRcommon_down[,1:3],c("BTBRmedia>LPS1&2"))
colnames(BTBRcommon_down_bed) <- c("chrom","start","end","condition")

C57LPSx1only_down_bed <- cbind(C57LPSx1only_down[,1:3],c("C57media>LPS1only"))
colnames(C57LPSx1only_down_bed) <- c("chrom","start","end","condition")

BTBRLPSx1only_down_bed <- cbind(BTBRLPSx1only_down[,1:3],c("BTBRmedia>LPS1only"))
colnames(BTBRLPSx1only_down_bed) <- c("chrom","start","end","condition")

C57LPSx2only_down_bed <- cbind(C57LPSx2only_down[,1:3],c("C57media>LPS2only"))
colnames(C57LPSx2only_down_bed) <- c("chrom","start","end","condition")

BTBRLPSx2only_down_bed <- cbind(BTBRLPSx2only_down[,1:3],c("BTBRmedia>LPS2only"))
colnames(BTBRLPSx2only_down_bed) <- c("chrom","start","end","condition")


#put all _bed dfs in a list
# Put them in a list
ComparisonsBeds <- lapply(ls(pattern = "*_bed"), get)

#set list names:
for (i in 1:length(ComparisonsBeds)) {
  print(i)
  names(ComparisonsBeds)[i] <- as.character(unique(ComparisonsBeds[[i]]$condition))
}

names(ComparisonsBeds) 

#write to bedfiles:
setwd("/Users/annieciernia/Box Sync/Ashwood_labcollaboration/Experiments/atac_2019/valR")

#make filenames
filenames <- paste(names(ComparisonsBeds),".mm10.bed", sep="")
filenames2 <- paste(names(ComparisonsBeds),".bed", sep="")
filenames2 <-gsub(">","greater",filenames2)
filenames2 <-gsub("<","less",filenames2)
filenames2 <-gsub("&","and",filenames2)

#write to bed
lapply(seq_along(ComparisonsBeds), function(i){
  write.table(ComparisonsBeds[[i]], file = filenames[i], col.names = F, row.names = F, quote=F, sep="\t")
})

lapply(seq_along(ComparisonsBeds), function(i){
  write.table(ComparisonsBeds[[i]], file = filenames2[i], col.names = F, row.names = F, quote=F, sep="\t")
})


#####################BTBR vs C57 comparisons##########################################
# Overlap testing regioneR ---------------------------------------

#get mouse genome
mouse.genome <- getGenomeAndMask("mm10", mask=NA)$genome
mouse.canonical <- filterChromosomes(mouse.genome, organism="mm")

#get lists
setwd("/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/consensus\ peaks")
path <- getwd()
files <- list.files(path=".", pattern="*.mm10.bed", all.files=T, full.names=T)
filelist <- lapply(files, readPeakFile)

#get and fix names of files
names_list <- paste0(basename(file_path_sans_ext(files)))
#names_list <- gsub("_DEpeaks.mm10", "", names_list)
names_list <- gsub("_L_", "<", names_list)
names_list <- gsub("_G_", ">", names_list)
names_list <- gsub("_DEpeaks.mm10", "", names_list)

names(filelist) <- names_list
names(filelist)

#number of regions in each list

#setwd("/Users/annieciernia/Box Sync/Ashwood_labcollaboration/Experiments/atac_2019/valR")

rl <- rapply(filelist, length, how="list")
rl <- as.data.frame(do.call(rbind, rl))
rl$comparison <- rownames(rl)
colnames(rl)[1] <- c("number of regions")
openxlsx::write.xlsx(rl,file="Number_of_regions.xlsx")

#define universe as all consesnsu peaks
universe <- readPeakFile("/Users/annieciernia/Sync/collaborations/Ashwood/BTBR_BMDM/atac_2019/DiffBind/Allconsensuspeaks.bed") 
#compare lists
set.seed(1)

#overlapping A with B > EVobs
#resample from universe (same size as list A)
#compare resampled list with list B > overlap of random sample with list B
#repeat 1000 times > distribution of EVperm
#gives permutation p value and Z-score
names(filelist)
C57files <- filelist[]
BTBRfiles <- filelist[]
Bfiles <- unlist(BTBRfiles[])

names(C57files)
names(Bfiles)

output <- NULL
for (i in 1:length(filelist)){
  sum <-NULL
  #loop through each file for comparison  
  for (B in 1:length(Bfiles)) {
    pt <- permTest(A=C57files[[i]],B=Bfiles[[B]], ntimes=1000,
                   alternative="greater",
                   randomize.function=resampleRegions,
                   universe=universe,
                   evaluate.function=numOverlaps, 
                   count.once=TRUE,
                   genome="mm10", mc.set.seed=FALSE, mc.cores=4)
    summary <- summary(pt)
    summary$condition1 <- paste(names(C57files[i]))
    summary$condition2 <- paste(names(Bfiles[B]))
    summary$length_condition1 <- length(C57files[[i]])
    summary$length_condition2 <- length(Bfiles[[B]])
    summary$overlap <- pt$numOverlaps$observed
    
    sum <- rbind(sum,summary)
    plot(pt)
  }
  
  output <- rbind(sum,output)
}

#remove overlaps with self
output2 <- output[which(output[,4] != output[,5]), ]

#unique combos: https://stackoverflow.com/questions/25297812/pair-wise-duplicate-removal-from-dataframe?noredirect=1&lq=1
cols = c(4,5)
newdf = output2[,cols]
for (i in 1:nrow(output2)){
  newdf[i, ] = sort(output2[i,cols])
}

output3 <- output2[!duplicated(newdf),]

#FDR correct comparisons
output3$FDR <- p.adjust(output3$pvalue, method="fdr")

write.csv(output3, file = "EnrichPeakOverlap_RegioneR_1_5_2020.csv")

