# Fluidigm HD 48x48 IFC data analysis
#10/30/18 AVC
#normalizes to average of wt media condition for each gene instead of within animal
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(lsmeans)
library(nlme)
library(xlsx)

################################################################
#8/17/16 data
################################################################
path = ("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/Fluidigm\ chip/8_17_16_BMDM_MG/")
setwd(path)

genelistinfo <- read.delim("genelistinfo.txt", header =T)

file <- read.csv("Ctvalues_BTBRBMDM_10_5_16_mod.csv",header =T)


#fix name of gene column
names(file)[names(file) == 'Name.1'] <- 'gene'


file.replace <- file[,1:9]

#replace all 999 low ct values with 30 (max Ct recorded)
file.replace$Ct.Value[file.replace$Ct.Value == 999.000000] <- 30


#######################################################
#BMDM gene expression
#######################################################
#filter out MG samples
dat <- file.replace %>%
  separate(col = Name, into = c("celltype","genotype","animal", "condition"), sep = " ", remove = F) %>%
  filter(celltype == "BMDM")


#calculate deltaCt between Ct gene-Ct hprt
samples <- unique(dat$Name) #unique samples
bind <- NULL

for(i in samples){
  
  print(i)
  #subset data for matching each sample/condition combo
  tmp <- dat %>%
    filter(Name == i) 

  #extract row/colmun position with hprt Ct value (col #9)
  hprt.Ct <- tmp[tmp$Type.1 == "Reference",12]
  
  #caluclate the relative ratio: Ct value - hprt Ct value
  tmp2 <- tmp %>%
  mutate(deltaCt = Ct.Value - hprt.Ct)

  tmp2 <- as.data.frame(tmp2)
  tmp2$deltaCt <-as.numeric(tmp2$deltaCt)
  bind <- rbind(bind, tmp2)
}

######################################

#calculate deltadeltaCt between media and LPSA, LPSA+P
genes <- unique(bind$gene) #unique samples

out <- NULL

for(i in genes){
  
  print(i)
  #return wt dct media value for that gene
  
  mediaddct <- bind %>%
    filter(gene == i) %>%
    filter(condition == "media") %>%
    filter(genotype == "C57")
  

  media_dct_ave <- mean(mediaddct$deltaCt) #mean of wt dct for that gene
  
  #subset data for matching each sample/condition combo
  tmp <- bind %>%
    filter(gene == i) %>%
    mutate(ddCt = deltaCt - media_dct_ave) %>%
    mutate(fold_change = 2^-ddCt)
    
    
    #select(Name,condition,deltaCt) %>%
    #spread(condition,deltaCt) %>%
  #mutate(ddCt.LPSA = LPS1 - media_dct_ave) %>%
  #mutate(ddCt.LPSA_P = LPS2 - media_dct_ave) %>%
  #mutate(ddCt.media = media - media_dct_ave) %>%
  #mutate(fold_change_media = 2^-ddCt.media) %>%
  #mutate(fold_change_LPSA = 2^-ddCt.LPSA)%>%
  #mutate(fold_change_LPSA_P = 2^-ddCt.LPSA_P) %>%
  #mutate(gene = paste(i))
  
  out <- rbind(out, tmp)
}  
#merge in gene list info
genelistinfo$Gene = tolower(genelistinfo$Gene)
data <- merge(out,genelistinfo, by.x = "gene", by.y ="Gene")

write.table(data,file="DeltaDeltaCt_8_17_16.txt",sep="\t", quote=FALSE)



################################################################
#5/9/16 data
################################################################
path = ("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/Fluidigm\ chip/5_9_16_chiprun")
setwd(path)

expinfo <- read.delim("expinfo.txt", header =T)
colnames(expinfo)[1] <- c("Name")

expinfo2 <- expinfo %>%
  mutate(Sample = paste(substr(Name,1,6))) %>%
  dplyr::select(Sample, genotype)

expinfo2 <- as.data.frame(expinfo2)




genelistinfo <- read.delim("genelistinfo.txt", header =T)

file <- read.csv("5_9_16_fludigm_resultsmod2.csv",header =T)

#file$Name.1 <- as.character(file$Name.1)
#file$Ct.Call <- as.character(file$Ct.Call)

#fix name of gene column
names(file)[names(file) == 'Name.1'] <- 'gene'

file.replace <- file[,1:9]

#replace all 999 low ct values with 30 (max Ct recorded)
file.replace$Ct.Value[file.replace$Ct.Value == 999.000000] <- 30

#average Ct values across technical replicates
file.tech <- file.replace %>%
  filter(Name != "NTC") %>%
  mutate(Sample = paste(substr(Name,1,6))) %>%
  mutate(condition = paste(substr(Name,8,8))) %>%
  group_by(gene,condition,Name,Type.1) %>%
  dplyr::summarize(tech.repl.count= n(), Ct.value.mean = mean(Ct.Value)) 

file.tech$condition[file.tech$condition == 1] <- "media"
file.tech$condition[file.tech$condition == 2] <- "LPS1"
file.tech$condition[file.tech$condition == 3] <- "LPS2"

#calculate deltaCt between Ct gene-Ct hprt
samples <- unique(file.tech$Name) #unique samples
bind <- NULL
file.tech$Ct.value.mean <- as.numeric(file.tech$Ct.value.mean)

for(i in samples){
  
  print(i)
  #subset data for matching each sample/condition combo
  tmp <- file.tech %>%
    filter(Name == i) 
  
  #extract row/colmun position with hprt Ct value (col #6)
  hprt.Ct <- tmp[tmp$Type.1 == "Reference",6]
  
  hprt.Ct <- as.numeric(hprt.Ct)
  #caluclate the relative ratio: Ct value - hprt Ct value
  tmp2 <- tmp %>%
    mutate(deltaCt = Ct.value.mean - hprt.Ct)
  
  tmp2 <- as.data.frame(tmp2)
  tmp2$deltaCt <-as.numeric(tmp2$deltaCt)
  bind <- rbind(bind, tmp2)
}

bind2 <- merge(bind, expinfo, by  = "Name")
#calculate deltadeltaCt between media and LPSA, LPSA+P
genes <- unique(bind$gene) #unique samples

out2 <- NULL

for(i in genes){
  
  print(i)
  #return wt dct media value for that gene
  
  mediaddct <- bind2 %>%
    filter(gene == i) %>%
    filter(condition == "media") %>%
    filter(genotype == "C57") 
  
  media_dct_ave <- mean(mediaddct$deltaCt) #mean of wt dct for that gene
  
  #subset data for matching each sample/condition combo
  tmp <- bind2 %>%
    filter(gene == i) %>%
    mutate(ddCt = deltaCt - media_dct_ave) %>%
    mutate(fold_change = 2^-ddCt)
  
  
  out2 <- rbind(out2, tmp)
}  
#merge in experiment info

out2 <- as.data.frame(out2)


genelistinfo$Gene = tolower(genelistinfo$Gene)

repl3 <- merge(out2, genelistinfo, by.x = "gene", by.y ="Gene")

write.table(repl3,file="DeltaDeltaCt_5_9_16.txt",sep="\t", quote=FALSE)

################################################################
################################################################

#add sample column
data$Sample <- paste(data$genotype, data$animal)

batch1 <- data %>%
  dplyr::select(gene,Name,genotype,condition,fold_change,deltaCt,class)


batch2 <- repl3 %>%
  dplyr::select(gene,Name,genotype,condition,fold_change,deltaCt,class)

batch1$batch <- c("batch1")
batch2$batch <- c("batch2")

alldat <- rbind(batch1,batch2)

setwd("/Users/annieciernia/Box\ Sync/Ashwood\ lab\ collaboration/Experiments/Fluidigm\ chip/Analysis_10_30_18")
write.xlsx(alldat,file="FluidigData_TwoBatches_BMDM_ET.xlsx")
########################################################################################################################
#add DDCT value and FC for within animal math
########################################################################################################################
genes <- unique(as.character(alldat$gene)) #unique samples

#sample names
alldat$SampleName <- alldat$Name
alldat$SampleName <- gsub("BMDM ","",alldat$SampleName)
alldat$SampleName <- gsub(" media","",alldat$SampleName)
alldat$SampleName <- gsub(" LPS1","",alldat$SampleName)
alldat$SampleName <- gsub(" LPS2","",alldat$SampleName)
alldat$SampleName <- gsub("C57 ","",alldat$SampleName)
alldat$SampleName <- gsub("\\.1","",alldat$SampleName)
alldat$SampleName <- gsub("\\.2","",alldat$SampleName)
alldat$SampleName <- gsub("\\.3","",alldat$SampleName)
unique(alldat$SampleName)


names <- unique(as.character(alldat$SampleName))

output <- NULL
for(i in genes){
  
  print(i)

  #subset data for matching each sample/condition combo
  tmp <- alldat %>%
    filter(gene == i) 
  
        gene1 <- NULL
        for (n in names){
          tmp2 <- tmp %>% filter(SampleName == n)
          
          mediadct <- tmp2 %>%
            filter(condition == "media") 
          
      
          #subset data for matching each sample/condition combo
          tmp2$ddCt <- tmp2$deltaCt - mediadct$deltaCt
          tmp2$fold_change <- 2^-tmp2$ddCt
          
          #output
          gene1 <- rbind(gene1,tmp2)
            
        }
        
  output <- rbind(output,gene1)
}
  
  

########################################################################################################################
#loop through each gene and calculate 2-way repeated measures ANOVA
#Save to table with animal count 
#make and save a graph of each gene 

data2 <- output %>%
  filter(gene != "hprt") #remove hprt (ref gene)
  
data2$gene <- as.factor(data2$gene)
data2$genotype <- as.factor(data2$genotype)
data2$condition <- as.factor(data2$condition)
data2$batch <- as.factor(data2$batch)

#relevel variables
data2$genotype <- relevel(data2$genotype, ref = "C57")
data2$condition <- relevel(data2$condition, ref = "media")


genes <- unique(data2$gene) #unique genes
######################################
setwd("/Users/aciernia/Sync/collaborations/Ashwood/BTBR_BMDM/Fluidigm/Analysis2020")
#complete model: 2 genotypes x 3 treatments with animal as repeated measure + batch
#stats on delta Ct but graphs of fold change
anova_fullmodel <- NULL
posthoc_fullmodel <- NULL
count_out <- NULL

#complete treatment and genotypes:
for(i in genes){

  print(i)
  
  #subset data for matching gene
  tmp <- data2 %>%
    filter(gene == i) #%>%
    #group_by(animal,genotype,treatment) %>%
    #summarize(FC=mean(FC))
  

  
  #count
  count <- tmp %>% group_by(genotype,condition) %>% summarize(count = n(), meanFC = mean(fold_change),meanFC = mean(deltaCt))
  count <- as.data.frame(count)
  count$gene <- paste(i)
  
  count_out <- rbind(count_out, count)
  
  
  #full model ANOVA:
  #batch as a covariate
      model <- lme(deltaCt ~ genotype * condition + batch, ~1|SampleName, data = tmp) #genotype included
      #sum <- summary(model)
     # plot(fitted(model),residuals(model))
      #qqnorm(residuals(model))
      
      anova = anova(model, type = "marginal")
      anova = as.data.frame(anova)
      anova$gene = paste(i)
      
      
      #posthocs for C57 vs BTBR within each treatment:
      refgrid <- ref.grid(model)
      posthocM <- lsmeans(refgrid, ~genotype|condition)
      C57vsBTBR_bytreat <- summary(pairs(posthocM, adjust = "none"))
      C57vsBTBR_bytreat <- as.data.frame(C57vsBTBR_bytreat)
      C57vsBTBR_bytreat$gene <- paste(i)
      colnames(C57vsBTBR_bytreat)[2] <- c("condition")
      
      #posthocs for C57 only media,LPS1, LPS2:
      refgrid <- ref.grid(model)
      posthocM <- lsmeans(refgrid, ~condition|genotype)
      C57vsBTBR <- summary(pairs(posthocM, adjust = "none"))
      C57vsBTBR <- as.data.frame(C57vsBTBR)
      C57vsBTBR$gene <- paste(i)
      colnames(C57vsBTBR)[2] <- c("condition")
      
      posthocout <- rbind(C57vsBTBR_bytreat,C57vsBTBR)
      #adjust p value for all comparisons:
      posthocout$p.adjust <- p.adjust(posthocout$p.value,method = "BH") #BH adjusted p values
      
      anova_fullmodel <- rbind(anova_fullmodel,anova)
      posthoc_fullmodel <- rbind(posthoc_fullmodel,posthocout)
      
      
      #plot
      cbPalette <- c("lightblue2","lightgoldenrod","palegreen3") #blue and orange http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
      library(cowplot)
      myplot <- ggplot(tmp, aes(x = genotype, y = log2(fold_change), fill=condition)) +
        #geom_boxplot(aes(fill=condition)) + 
        #wiskers = 5th and 95th percentile.https://stackoverflow.com/questions/28961781/changing-whisker-definition-in-faceted-geom-boxplot
        stat_summary(geom = "boxplot", 
                     fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
                     position=position_dodge())+
        geom_point(position = position_dodge(width = 0.90),aes(group=condition),size=3)+
        scale_y_continuous(name=paste("Log2 Fold Change Relative C57 Media",i, sep = " ")) +
        scale_x_discrete(name = " ")+
        theme(legend.title = element_blank())+  
        scale_fill_manual(values = cbPalette) +
        theme(strip.text.x = element_text(size = 20)) +
        #theme(legend.key = element_blank(), strip.background = element_blank()) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        theme(axis.text=element_text(size=20),axis.title=element_text(size=20),legend.text=element_text(size=20))
      
      
      #save
      ggsave(myplot,height = 5, width = 5, filename = paste(i,"Log2foldchange_Fluidigm.pdf",sep="."))
      
      
}



#write output files
write.xlsx(anova_fullmodel,file="Fluidigm_ANOVA.xlsx")
write.xlsx(posthoc_fullmodel,file="Fluidigm_posthocs.xlsx")
write.xlsx(count_out,file="Fluidigm_CountperCondition.xlsx")

  

#########################################################################################
#catagorize as Tolerized or Non-tolerized
#########################################################################################
#means by condition
data2$treatment <- data2$condition
data2$strain <- data2$genotype

outdf <- data2 %>% dplyr::select(gene,strain,treatment,fold_change,class) %>%
  group_by(strain,treatment,gene) %>%
  summarise(meanFC = mean(fold_change)) 
  

#calculate mean differences
tmp <- outdf %>%
  group_by(gene) %>%
  spread(treatment,meanFC)

tmp <- as.data.frame(tmp)

#media_vs_LPS1
tmp$media_vs_LPS1 <- tmp$LPS1 - tmp$media

#media_vs_LPS2
tmp$media_vs_LPS2 <- tmp$LPS2 - tmp$media

#LPS1vsLPS2
tmp$LPS1_vs_LPS2 <- tmp$LPS1 - tmp$LPS2

#add in stats data from posthocs

#comparsion
posthoc_fullmodel$comparison <- paste(as.character(posthoc_fullmodel$condition),as.character(posthoc_fullmodel$contrast),sep= " ")

keepcomparisons <-c("media C57 - BTBR",
                    "C57 media - LPS1",
                    "C57 media - LPS2",
                    "C57 LPS1 - LPS2",
                    "BTBR media - LPS1",
                    "BTBR media - LPS2",
                    "BTBR LPS1 - LPS2",
                    "LPS1 C57 - BTBR",
                    "LPS2 C57 - BTBR")

posthoc_fullmodel2 <- subset(posthoc_fullmodel, comparison %in% keepcomparisons)
posthoc_fullmodel2 <- as.data.frame(posthoc_fullmodel2)

#filter
anova_fullmodel3 <- posthoc_fullmodel2 %>% 
  dplyr::select(comparison,gene,p.value) %>%
  group_by(gene) %>%
  spread(comparison,p.value)


MasterFilter <- merge(tmp, anova_fullmodel3, by = "gene")


#for C57
MasterFilterC57 <- MasterFilter %>% filter(strain == "C57")

#set up catagories
#induced = media_vs_LPS1 > 0 and "C57,media - C57,LPS1" <0.05
MasterFilterC57$induced <- ifelse(MasterFilterC57$`C57 media - LPS1` <0.05 & MasterFilterC57$media_vs_LPS1 > 0,"induced","not")

#tolerized: induced, LPS1_vs_LPS2 > 0, and C57,LPS1 - C57,LPS2 <0.05 
MasterFilterC57$tolerized <- ifelse(MasterFilterC57$`C57 LPS1 - LPS2` <0.05 & MasterFilterC57$induced == "induced" & MasterFilterC57$LPS1_vs_LPS2 > 0,"tolerized","not")

#sensitized
#LPS1vsLPS2 is < 0 (LPS1- LPS2 = negative)
MasterFilterC57$sensitized <- ifelse(MasterFilterC57$`C57 LPS1 - LPS2` <0.05 & MasterFilterC57$induced == "induced" & MasterFilterC57$LPS1_vs_LPS2 < 0,"sensitized","not")

#non-tolerized:C57,LPS1 - C57,LPS2` not significant
MasterFilterC57$nontolerized <- ifelse(MasterFilterC57$`C57 LPS1 - LPS2` >0.05 & MasterFilterC57$induced == "induced","nontolerized","not")


write.xlsx(MasterFilterC57,"MasterFilterC57_catagorization.xlsx")


#for BTBR
MasterFilterBTBR <- MasterFilter %>% filter(strain == "BTBR")

#set up catagories
#induced = media_vs_LPS1 > 0 and "BTBR,media - BTBR,LPS1" <0.05
MasterFilterBTBR$induced <- ifelse(MasterFilterBTBR$`BTBR media - LPS1` <0.05 & MasterFilterBTBR$media_vs_LPS1 > 0,"induced","not")

#tolerized: induced, LPS1_vs_LPS2 > 0, and BTBR,LPS1 - BTBR,LPS2 <0.05 
MasterFilterBTBR$tolerized <- ifelse(MasterFilterBTBR$`BTBR LPS1 - LPS2` <0.05 & MasterFilterBTBR$induced == "induced" & MasterFilterBTBR$LPS1_vs_LPS2 > 0,"tolerized","not")

#sensitized
#LPS1vsLPS2 is < 0 (LPS1- LPS2 = negative)
MasterFilterBTBR$sensitized <- ifelse(MasterFilterBTBR$`BTBR LPS1 - LPS2` <0.05 & MasterFilterBTBR$induced == "induced" & MasterFilterBTBR$LPS1_vs_LPS2 < 0,"sensitized","not")

#non-tolerized:BTBR,LPS1 - BTBR,LPS2` not significant
MasterFilterBTBR$nontolerized <- ifelse(MasterFilterBTBR$`BTBR LPS1 - LPS2` >0.05 & MasterFilterBTBR$induced == "induced","nontolerized","not")

write.xlsx(MasterFilterBTBR,"MasterFilterBTBR_catagorization.xlsx")

#combinedoutput
combinedoutput <- rbind(MasterFilterC57,MasterFilterBTBR)
combinedoutput$mediaDifferences <- ifelse(combinedoutput$`media C57 - BTBR` < 0.05,"baseline differences","no baseline differences")

#master catagory for toler, nontol, sens, none
combinedoutput$Catagory <- combinedoutput$nontolerized 

combinedoutput$Catagory[combinedoutput$Catagory == "not" & combinedoutput$tolerized == "tolerized"] <- c("tolerized")

combinedoutput$Catagory[combinedoutput$Catagory == "not" & combinedoutput$sensitized == "sensitized"] <- c("sensitized")

combinedoutput$Catagory[combinedoutput$induced == "not"] <- c("not induced LPS1")




write.xlsx(combinedoutput,"MasterFilterAllData_catagorization.xlsx")

#make catagories for plotting from media differences and C57 
Allcatagories <-combinedoutput %>%
  dplyr::select(gene,strain,Catagory, mediaDifferences) %>%
  mutate(StrainCatagory =  paste(strain,Catagory, sep =" ")) %>%
  dplyr::select(gene,StrainCatagory,mediaDifferences)


# plots -------------------------------------------------------------------

#average FC by condition
outFC <- outdf %>% ungroup() %>% # fixes Adding missing grouping variables: `strain`, `treatment` error
  mutate(Condition =  paste(strain,treatment)) %>%
  dplyr::select(-strain,-treatment) %>% 
  group_by(gene) %>%
  spread(Condition,meanFC)

#add catagories 
outFCmerge <- merge(outFC, Allcatagories, by = c("gene"))

# outFCmerge$Catagory <- factor(outFCmerge$Catagory, levels = c("not induced LPS1","tolerized","nontolerized","sensitized"))
# 
# outFCmerge$mediaDifferences <- factor(outFCmerge$mediaDifferences, levels = unique(heatDF$mediaDifferences))
# 
# heatDF3 <- outFCmerge %>% arrange(gene,Sex,Strain,mediaDifferences,Catagory) %>%
#   mutate(StrainSex = paste(Strain,Sex)) %>%
#   mutate(Condition =  paste(Strain, Sex, treatment)) %>%
#   arrange(gene,StrainSex,Catagory)


#for catagory
heatmatrixFC2 <- outFCmerge %>%
  dplyr::select(gene,StrainCatagory) %>%
  distinct() %>%
  arrange(gene,StrainCatagory) %>%
  group_by(gene) %>%
  mutate(Order = seq_along(gene)) %>%
  spread(key = Order, value = StrainCatagory)

colnames(heatmatrixFC2) <- c("gene","BTBR Catagory","C57 Catagory")

heatmatrixFC <- distinct(outFCmerge[,c(1:7,9)])
tmp2 <- distinct(combinedoutput[,c(1,9:17)])

MasterHeat <- merge(heatmatrixFC,heatmatrixFC2, by="gene")

MasterHeat2 <- merge(MasterHeat,tmp2, by="gene")

#hyper respond to LPS1
MasterHeat2$HyperLPSx1 <- ifelse(MasterHeat2$`LPS1 C57 - BTBR` < 0.05 & MasterHeat2$`BTBR LPS1`-MasterHeat2$`C57 LPS1`>0,"BTBR Hyper LPSx1","Normal")
MasterHeat2$HyperLPSx2 <- ifelse(MasterHeat2$`LPS2 C57 - BTBR` < 0.05 & MasterHeat2$`BTBR LPS2`-MasterHeat2$`C57 LPS2`>0,"BTBR Hyper LPSx2","Normal")


write.xlsx(MasterHeat2,"MasterHeatmapdata.xlsx")

MasterHeat2 <- MasterHeat2 %>% arrange(desc(`C57 Catagory`),desc(`BTBR Catagory`),
                                       desc(HyperLPSx1),desc(HyperLPSx2),`mediaDifferences`)

#matrix for heatmap
heatDF2 <- as.data.frame(MasterHeat2[,2:7])

#reorder columns
colorder <- c("C57 media","C57 LPS1", "C57 LPS2","BTBR media","BTBR LPS1" , "BTBR LPS2")

rownames(heatDF2) <- MasterHeat2$gene

heatDF2 <- heatDF2[,colorder]

heatDF2 <- data.matrix(heatDF2)

##=====================================================================================
#all samples
##=====================================================================================

#heatmap:
library(pheatmap)

#define row groups:
annotation_row <- MasterHeat2[,c(8,10,9,20,21)]
colnames(annotation_row) <- c("Media Baseline","C57 Category","BTBR Category","BTBR Hyper LPSx1","BTBR Hyper LPSx2")

rownames(annotation_row) = rownames(heatDF2)

head(annotation_row)
#col


annotation_col <- as.data.frame(colnames(heatDF2))
names(annotation_col) <- c("col")
annotation_col$col <- as.character(annotation_col$col)

annotation_col <- tidyr::separate(annotation_col,col,into=c("Strain","LPS Treatment"),sep="\\.")

rownames(annotation_col) = colnames(heatDF2)

head(annotation_col)

# Specify colors
# #https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
# library(RColorBrewer)
# n <- 27
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# #pie(rep(1,n), col=sample(col_vector, n))

#gradient
#colfunc_blue <- colorRampPalette(c("lightblue", "darkblue"))
#colfunc_green <- colorRampPalette(c("lightgreen", "darkgreen"))

#media,LPS1,LPS2
#cbPalette <- c("royalblue3","darkorange1","brown4") #blue and orange http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# treatment
Var1        <- c("royalblue3","coral1","darkorchid4")
names(Var1) <- unique(annotation_col$`LPS Treatment`)

#strain
Var3        <- c("darkseagreen4","blueviolet")
names(Var3) <- unique(annotation_col$Strain)

Var4 <- c("darkred","gold","grey","darkorange")
names(Var4) <- unique(annotation_row$`C57 Category`)

Var5 <- c("darkred","darkorange","gold")
names(Var5) <- unique(annotation_row$`BTBR Category`)

Var8 <- c("cyan3","darkslategray4")
names(Var8) <- unique(annotation_row$`Media Baseline`)

Var9 <- c("gray","firebrick3")
names(Var9) <- unique(annotation_row$`BTBR Hyper LPSx1`)

Var10 <- c("gray","deeppink3")
names(Var10) <- unique(annotation_row$`BTBR Hyper LPSx2`)


#combined for heatmap labels
anno_colors <- list(`LPS Treatment`=Var1,Strain=Var3,
                    `C57 Category`=Var4,
                    `BTBR Category`=Var5,
                    `Media Baseline`=Var8,
                    `BTBR Hyper LPSx1`=Var9,
                    `BTBR Hyper LPSx2`=Var10)


#my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

pdf(file = "Heatmap_Fluidigm_FCexpression_cluster.pdf", wi = 8, he = 12)
pheatmap::pheatmap(heatDF2, 
                   cluster_row = T,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = T,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=16, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Fold Change Gene Expression")

dev.off()

##=====================================================================================
#complex heatmap
#mat <- read.xlsx("MasterHeatmapdata.xlsx",sheetIndex = 1)
mat <- MasterHeat2

mat <- mat %>% filter(`C57 Catagory` != "C57 not induced LPS1") %>%
  filter(`C57 Catagory` != "C57 sensitized") %>%
  arrange(`C57 Catagory`,`BTBR Catagory`,mediaDifferences,HyperLPSx1,HyperLPSx2)

mat$`C57 Catagory` <- as.factor(as.character(mat$`C57 Catagory`))
rownames(mat) <- mat$gene


pval <- mat[,c(11:21)]
pval <- pval[,grepl("media",colnames(pval))]
pval <- pval[,c(3,4,1,2)]


rownames(pval) <- mat$gene
name <- rownames(pval)
rownames(pval) <- paste(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name)), sep="")


library(ComplexHeatmap)

library(dendsort)


#reorder columns
mat2 <-  mat[,c(5,6,2,3)]

rownames(mat2) <- mat$gene
name <- rownames(mat2)
rownames(mat2) <- paste(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name)), sep="")

#dataRowNorm <- t(apply(mat2, 1, function(x)  (x - mean(x)) / sd(x)))
dataRowNorm <- log2(as.matrix(mat2))
#mat[is.na(mat)] <- 0

dend = dendsort(hclust(dist(dataRowNorm )))


#row order for clustering
rownames(dataRowNorm )
mat$col.cat <- paste(mat$`C57 Catagory`,mat$`BTBR Catagory`, sep=".")

max <- max(dataRowNorm [,c(1:4)], na.rm=T)
min <- min(dataRowNorm [,c(1:4)], na.rm=T)

library(circlize)
col_fun = colorRamp2(c(min,max), c("white",  "red"))

#row annotation by TF binding site
row_ha = rowAnnotation(df = mat[,c(9,8,20,21)])

#ha <- HeatmapAnnotation(region = collab,
 #                       col = list(region= c("HC" = "blue","PFC" = "orange",  "STR" = "purple")),
 #                       gp = gpar(col = "black"))

ht0 = Heatmap(dataRowNorm, name = "Log2 Fold Change", #cluster_rows = dend, 
              width = unit(6, "cm"), height = unit(16, "cm"),
              column_order = colnames(dataRowNorm), 
              row_order = rownames(dataRowNorm),
              row_split = mat$`C57 Catagory`, #number of clusters
              col = col_fun, 
              right_annotation = row_ha,
              #add pvalues
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   if(pval[i, j] <0.05)
              #     #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
              #     grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              # }
             )
draw(ht0) 
pdf("Log2FCcomplexheat.pdf", width=8, height=16)

draw(ht0)

dev.off()

#sc
##=====================================================================================
#split heatmap by col for strain
##=====================================================================================
library(circlize)

C57data <- dataRowNorm[,c(1:3)]
BTBRdata <- dataRowNorm[,c(4:6)]

col_rnorm = colorRamp2(c(min(BTBRdata), 0, max(BTBRdata)), c("blue", "white", "red"))

ht1 = Heatmap(C57data, name = "C57 Log2 Fold Change", #cluster_rows = dend, 
              width = unit(4, "cm"),
              column_order = colnames(C57data), row_split = mat$C57.Catagory, #number of clusters
              rect_gp = gpar(col = "grey", lwd = .5),col = col_rnorm
              # top_annotation = ha,
              #add pvalues
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   if(pval[i, j] <0.05)
              #     #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
              #     grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              # }
)



ht2 = Heatmap(BTBRdata, name = "BTBR Log2 Fold Change", #cluster_rows = dend, 
              width = unit(4, "cm"),
              column_order = colnames(BTBRdata), row_split = mat$BTBR.Catagory, #number of clusters
              rect_gp = gpar(col = "grey", lwd = .5), col = col_rnorm
              # top_annotation = ha,
              #add pvalues
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   if(pval[i, j] <0.05)
              #     #grid.text(sprintf("%.1f", adjpval[i, j]), x, y, gp = gpar(fontsize = 4))
              #     grid.points(x, y, pch = 8, size = unit(2, "mm")) #pch =8  for *
              # }
)


ht_list = ht1 + ht2

pdf("complexheatmap_FluidigmLog2FC.pdf", width=10, height=8)

draw(ht_list, row_title = " ", row_title_gp = gpar(col = "black"),
     column_title = " ", column_title_gp = gpar(fontsize = 16))

dev.off()

#counts
df = mat[,c(10,9,8,20,21)]
k <- as.data.frame(table(df))

write.csv(k,"Counts_catagories.csv")

