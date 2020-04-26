#!/bin/bash

module load deeptools/3.1.0

#########################################################################################

#Make bigwig files with no duplicates, extend reads between mates, and normalized for read depth
# Report read coverage normalized to 1x sequencing depth (also known as Reads Per Genomic Content (RPGC)). 
#Sequencing depth is defined as: (total number of mapped reads * fragment length) / effective genome size. 
# The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. 
#To use this option, the effective genome size has to be indicated after the option. The effective genome size is the portion of the genome that is mappable. 
#Large fractions of the genome are stretches of NNNN that should be discarded. Also, if repetitive regions were not included in the mapping of reads, the effective genome size needs to be adjusted accordingly. 


############################################################################################bamCoverage --bam /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003A.NFRchrMfilter.sort.bam -o /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003A.NFRSeqDepthNorm_mm10.bw --ignoreDuplicates --binSize 10 --effectiveGenomeSize 2652783500 --ignoreForNormalization chrM chrX chrY --extendReads -p 30 --centerReads --minFragmentLength 0 --maxFragmentLength 150 --blackListFileName /share/lasallelab/Annie/BMDM/ATACanalysis2019/deeptools/mm10.blacklist.bed

# for sample in $(<AllSamples.txt) 
# 	do
#    	#do stuff
# 	echo ${sample}
# 	
# 	bamCoverage --bam /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/${sample}.NFRchrMfilter.sort.bam -o /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/${sample}.NFRSeqDepthNorm.bw --ignoreDuplicates --binSize 10 --effectiveGenomeSize 2652783500 --ignoreForNormalization chrM chrX chrY --extendReads -p 20 --centerReads --minFragmentLength 0 --maxFragmentLength 150 --blackListFileName /share/lasallelab/Annie/BMDM/ATACanalysis2019/deeptools/mm10.blacklist.bed
# 
# 
# done

##########################################################################################
# multibigwig summary
# 
# multiBigwigSummary bins -b /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006C.NFRSeqDepthNorm_mm10.bw 
# 
# /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006D.NFRSeqDepthNorm_mm10.bw 
# 
# /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006E.NFRSeqDepthNorm_mm10.bw 
# 
# /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006F.NFRSeqDepthNorm_mm10.bw --labels C57media1 C57media2 C57media3 C57media4 C571xLPS1 C571xLPS2 C571xLPS3 C571xLPS4 C572xLPS1 C572xLPS2 C572xLPS3 C572xLPS4 BTBRmedia1 BTBRmedia2 BTBRmedia3 BTBRmedia4 BTBR1xLPS1 BTBR1xLPS2 BTBR1xLPS3 BTBR1xLPS4 BTBR2xLPS1 BTBR2xLPS2 BTBR2xLPS3 BTBR2xLPS4 -out Allscores_per_bin.npz --outRawCounts Allscores_per_bin.tab --numberOfProcessors 20
# 

#plot PCA
#plotPCA -in Allscores_per_bin.npz -o PCA_AllreadCounts.pdf -T "PCA of Sequencing Depth Normalized Read Counts" --labels C57media1 C57media2 C57media3 C57media4 C571xLPS1 C571xLPS2 C571xLPS3 C571xLPS4 C572xLPS1 C572xLPS2 C572xLPS3 C572xLPS4 BTBRmedia1 BTBRmedia2 BTBRmedia3 BTBRmedia4 BTBR1xLPS1 BTBR1xLPS2 BTBR1xLPS3 BTBR1xLPS4 BTBR2xLPS1 BTBR2xLPS2 BTBR2xLPS3 BTBR2xLPS4 

##########################################################################################
#individual conditions
##########################################################################################

#C57media
echo "C57media"
#multibigwig summary
multiBigwigSummary bins -b /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006A.NFRSeqDepthNorm_mm10.bw -out C57media_multibigwigsum.npz --outRawCounts C57media_multibigwigsum.tab --numberOfProcessors 20

#plot PCA
plotPCA -in C57media_multibigwigsum.npz -o PCA_C57media_multibigwigsum.pdf -T "PCA of Sequencing Depth Normalized Read Counts" --plotHeight 7 --plotWidth 9

plotCorrelation -in C57media_multibigwigsum.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Sequencing Depth Normalized Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o SpearmanCorr_C57media_multibigwigsum.pdf

##########################################################################################
#C57LPSx1
echo "C57LPS1"
multiBigwigSummary bins -b /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006B.NFRSeqDepthNorm_mm10.bw -out C57x1LPS_multibigwigsum.npz --outRawCounts C57x1LPS_multibigwigsum.tab --numberOfProcessors 20


#plot PCA
plotPCA -in C57x1LPS_multibigwigsum.npz -o PCA_C57x1LPS_multibigwigsum.pdf -T "PCA of Sequencing Depth Normalized Read Counts" --plotHeight 7 --plotWidth 9

plotCorrelation -in C57x1LPS_multibigwigsum.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Sequencing Depth Normalized Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o SpearmanCorr_C57x1LPS_multibigwigsum.pdf

##########################################################################################
#C57LPSx2
echo "C57LPS2"
multiBigwigSummary bins -b /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006C.NFRSeqDepthNorm_mm10.bw -out C57x2LPS_multibigwigsum.npz --outRawCounts C57x2LPS_multibigwigsum.tab --numberOfProcessors 20

#plot PCA
plotPCA -in C57x2LPS_multibigwigsum.npz -o PCA_C57x2LPS_multibigwigsum.pdf -T "PCA of Sequencing Depth Normalized Read Counts" --plotHeight 7 --plotWidth 9

plotCorrelation -in C57x2LPS_multibigwigsum.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Sequencing Depth Normalized Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o SpearmanCorr_C57x2LPS_multibigwigsum.pdf

##########################################################################################
#BTBR media
echo "BTBRmedia"
multiBigwigSummary bins -b /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006D.NFRSeqDepthNorm_mm10.bw  -out BTBRmedia_multibigwigsum.npz --outRawCounts BTBRmedia_multibigwigsum.tab --numberOfProcessors 20

#plot PCA
plotPCA -in BTBRmedia_multibigwigsum.npz -o PCA_BTBRmedia_multibigwigsum.pdf -T "PCA of Sequencing Depth Normalized Read Counts" --plotHeight 7 --plotWidth 9


plotCorrelation -in BTBRmedia_multibigwigsum.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Sequencing Depth Normalized Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o SpearmanCorr_BTBRmedia_multibigwigsum.pdf

##########################################################################################
#BTBR 1xLPS
echo "BTBRLPS1"
multiBigwigSummary bins -b /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006E.NFRSeqDepthNorm_mm10.bw  -out BTBRx1LPS_multibigwigsum.npz --outRawCounts BTBRx1LPS_multibigwigsum.tab --numberOfProcessors 20

#plot PCA
plotPCA -in BTBRx1LPS_multibigwigsum.npz -o PCA_BTBRx1LPS_multibigwigsum.pdf -T "PCA of Sequencing Depth Normalized Read Counts" --plotHeight 7 --plotWidth 9

plotCorrelation -in BTBRx1LPS_multibigwigsum.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Sequencing Depth Normalized Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o SpearmanCorr_BTBRx1LPS_multibigwigsum.pdf
##########################################################################################
#BTBR 2xLPS
echo "BTBRLPS2"
multiBigwigSummary bins -b /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006F.NFRSeqDepthNorm_mm10.bw -out BTBRx2LPS_multibigwigsum.npz --outRawCounts BTBRx2LPS_multibigwigsum.tab --numberOfProcessors 20

#plot PCA
plotPCA -in BTBRx2LPS_multibigwigsum.npz -o PCA_BTBRx2LPS_multibigwigsum.pdf -T "PCA of Sequencing Depth Normalized Read Counts" --plotHeight 7 --plotWidth 9

plotCorrelation -in BTBRx2LPS_multibigwigsum.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Sequencing Depth Normalized Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o SpearmanCorr_BTBRx2LPS_multibigwigsum.pdf
#######################################################################################