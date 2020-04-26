#!/bin/bash
#

##########################################################################################
module load  kentutils/302.0.0 
module load deeptools/3.1.0
#merge by condition into one bigwig
##########################################################################################
#http://wresch.github.io/2014/01/31/merge-bigwig-files.html
#merge big wigs with normalization so that each replicate contributes equally to the final average
	
./merge_bigwig_server.sh C57media_mm10.bw mm10.chrom.sizes /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006A.NFRSeqDepthNorm_mm10.bw 
# 
./merge_bigwig_server.sh C57LPS1_mm10.bw mm10.chrom.sizes /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006B.NFRSeqDepthNorm_mm10.bw
# 
./merge_bigwig_server.sh C57LPS2_mm10.bw mm10.chrom.sizes /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006C.NFRSeqDepthNorm_mm10.bw
# 
./merge_bigwig_server.sh BTBRmedia_mm10.bw mm10.chrom.sizes /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003D.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005A.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006D.NFRSeqDepthNorm_mm10.bw
# 
./merge_bigwig_server.sh BTBRLPS1_mm10.bw mm10.chrom.sizes /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003E.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005B.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006E.NFRSeqDepthNorm_mm10.bw
# 
./merge_bigwig_server.sh BTBRLPS2_mm10.bw mm10.chrom.sizes /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC003F.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAC004C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC005C.NFRSeqDepthNorm_mm10.bw /share/lasallelab/Annie/BMDM/ATACanalysis2019/alignment/JLAVC006F.NFRSeqDepthNorm_mm10.bw

##########################################################################################

#computeMatrix scale-regions -S Lbiwig file(s)G -R Lbed fileG -b 1000
computeMatrix reference-point --referencePoint center -S C57media_mm10.bw C57LPS1_mm10.bw C57LPS2_mm10.bw BTBRmedia_mm10.bw BTBRLPS1_mm10.bw BTBRLPS2_mm10.bw --samplesLabel C57media C57LPSx1 C57LPSx2 BTBRmedia BTBRLPSx1 BTBRLPSx2 -b 500 -a 500 --skipZeros -p 10 --outFileSortedRegions Mediamatrix_sigpeaks.bed -o Mediamatrix_sigpeaks.gz -R BTBRmedia_L_C57media.mm10.bed BTBRmedia_G_C57media.mm10.bed

# with legend media ymax chosen based on initial graphs
plotProfile -m Mediamatrix_sigpeaks.gz --numPlotsPerRow 1 --regionsLabel "BTBR media < C57 media" "BTBR media > C57 media" --plotFileFormat "pdf" --samplesLabel "C57 media" "C57 LPSx1" "C57 LPSx2" "BTBR media" "BTBR LPSx1" "BTBR LPSx2" -out profile_AverageBW_mediasigpeaks_G.pdf --averageType mean --perGroup --yMax 40
# 

#same but heatmap
plotHeatmap -m Mediamatrix_sigpeaks.gz --regionsLabel "BTBR media < C57 media" "BTBR media > C57 media" --plotFileFormat "pdf" --samplesLabel "C57 media" "C57 LPSx1" "C57 LPSx2" "BTBR media" "BTBR LPSx1" "BTBR LPSx2" -out heatmap_AverageBW_mediasigpeaks_G.pdf --outFileNameMatrix heatmapt_AverageBW_mediasigpeaks --refPointLabel "peak center"  --averageTypeSummaryPlot mean --perGroup --colorMap RdBu --yMax 40

  
# media>LPS ymax chosen based on initial graphs
   
computeMatrix reference-point --referencePoint center -S C57media_mm10.bw C57LPS1_mm10.bw C57LPS2_mm10.bw BTBRmedia_mm10.bw BTBRLPS1_mm10.bw BTBRLPS2_mm10.bw --samplesLabel C57media C57LPSx1 C57LPSx2 BTBRmedia BTBRLPSx1 BTBRLPSx2 -b 500 -a 500 --skipZeros -p 20 --outFileSortedRegions Greater_matrix_sigpeaks.bed -o Greater_matrix_sigpeaks.gz -R C57media_G_LPS1only.mm10.bed BTBRmedia_G_LPS1only.mm10.bed C57media_G_LPS2only.mm10.bed BTBRmedia_G_LPS2only.mm10.bed C57media_G_LPS1and2.mm10.bed BTBRmedia_G_LPS1and2.mm10.bed BTBRLPS1_G_C57LPS1.mm10.bed BTBRLPS2_G_C57LPS2.mm10.bed BTBRLPS1_L_C57LPS1.mm10.bed BTBRLPS2_L_C57LPS2.mm10.bed


plotProfile -m Greater_matrix_sigpeaks.gz --numPlotsPerRow 2 --regionsLabel "C57 LPSx1 < media" "BTBR LPSx1 < media" "C57 LPSx2 < media" "BTBR LPSx2 < media" "C57 media < LPSx1&2" "BTBR media < LPSx1&2" "BTBR LPSx1 > C57 LPSx1" "BTBR LPSx2 > C57 LPSx2" "BTBR LPSx1 < C57 LPSx1" "BTBR LPSx2 < C57 LPSx2" --plotFileFormat "pdf" --samplesLabel "C57 media" "C57 LPSx1" "C57 LPSx2" "BTBR media" "BTBR LPSx1" "BTBR LPSx2" -out profile_AverageBW_mediaGreaterLPSsigpeaks_G.pdf --averageType mean --perGroup --yMax 40
# 

plotHeatmap -m Greater_matrix_sigpeaks.gz --regionsLabel "C57 LPSx1 < media" "BTBR LPSx1 < media" "C57 LPSx2 < media" "BTBR LPSx2 < media" "C57 media < LPSx1&2" "BTBR media < LPSx1&2" "BTBR LPSx1 > C57 LPSx1" "BTBR LPSx2 > C57 LPSx2" "BTBR LPSx1 < C57 LPSx1" "BTBR LPSx2 < C57 LPSx2" --plotFileFormat "pdf" --samplesLabel "C57 media" "C57 LPSx1" "C57 LPSx2" "BTBR media" "BTBR LPSx1" "BTBR LPSx2" -out heatmap_AverageBW_mediaGreaterLPS.pdf --outFileNameMatrix heatmap_AverageBW_mediaGreaterLPSsigpeaks.gz --refPointLabel "peak center"  --averageTypeSummaryPlot mean --perGroup --colorMap RdBu --yMax 40


# media<LPS ymax chosen based on initial graphs
computeMatrix reference-point --referencePoint center -S C57media_mm10.bw C57LPS1_mm10.bw C57LPS2_mm10.bw BTBRmedia_mm10.bw BTBRLPS1_mm10.bw BTBRLPS2_mm10.bw --samplesLabel C57media C57LPSx1 C57LPSx2 BTBRmedia BTBRLPSx1 BTBRLPSx2 -b 500 -a 500 --skipZeros -p 20 --outFileSortedRegions Less_matrix_sigpeaks.bed -o Less_matrix_sigpeaks.gz -R C57media_L_LPS1only.mm10.bed BTBRmedia_L_LPS1only.mm10.bed C57media_L_LPS2only.mm10.bed BTBRmedia_L_LPS2only.mm10.bed C57media_L_LPS1and2.mm10.bed BTBRmedia_L_LPS1and2.mm10.bed BTBRLPS1_G_C57LPS1.mm10.bed BTBRLPS2_G_C57LPS2.mm10.bed BTBRLPS1_L_C57LPS1.mm10.bed BTBRLPS2_L_C57LPS2.mm10.bed


plotProfile -m Less_matrix_sigpeaks.gz --numPlotsPerRow 2 --regionsLabel "C57 media < LPSx1" "BTBR media < LPSx1" "C57 media < LPSx2" "BTBR media < LPSx2" "C57 media < LPSx1&2" "BTBR media < LPSx1&2" "BTBR LPSx1 > C57 LPSx1" "BTBR LPSx2 > C57 LPSx2" "BTBR LPSx1 < C57 LPSx1" "BTBR LPSx2 < C57 LPSx2" --plotFileFormat "pdf" --samplesLabel "C57 media" "C57 LPSx1" "C57 LPSx2" "BTBR media" "BTBR LPSx1" "BTBR LPSx2" -out profile_AverageBW_mediaLessLPSsigpeaks_G.pdf --averageType mean --perGroup --yMax 40
# 

plotHeatmap -m Less_matrix_sigpeaks.gz --regionsLabel "C57 media < LPSx1" "BTBR media < LPSx1" "C57 media < LPSx2" "BTBR media < LPSx2" "C57 media < LPSx1&2" "BTBR media < LPSx1&2" "BTBR LPSx1 > C57 LPSx1" "BTBR LPSx2 > C57 LPSx2" "BTBR LPSx1 < C57 LPSx1" "BTBR LPSx2 < C57 LPSx2" --plotFileFormat "pdf" --samplesLabel "C57 media" "C57 LPSx1" "C57 LPSx2" "BTBR media" "BTBR LPSx1" "BTBR LPSx2" -out heatmap_AverageBW_mediaLessLPSsigpeaks.pdf --outFileNameMatrix heatmap_AverageBW_mediaLessLPSsigpeaks.gz --refPointLabel "peak center"  --averageTypeSummaryPlot mean --perGroup --colorMap RdBu --yMax 40
