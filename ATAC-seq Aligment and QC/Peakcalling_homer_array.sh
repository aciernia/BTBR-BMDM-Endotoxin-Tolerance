#!/bin/bash

#SBATCH --array=1-24  # NEED TO CHANGE THIS!
#SBATCH --job-name=Tagshomer # Job name
#SBATCH --workdir /share/lasallelab/Annie/BMDM/ATACanalysis2019/
#SBATCH -p production              
#SBATCH -c 2    
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=0-02:00:00
#SBATCH --output=arrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=arrayJob_%A_%a.err # File to which STDERR will be written
#######################################################################################

begin=`date +%s`

echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" AllSamples.txt`

#######################################################################################
module load homer/4.9
module load kentutils/302.0.0
module load bedtools2/2.27.0 
#######################################################################################
#make tag directories for each bam file
#makeTagDirectory TagDirectory/tag_${sample} alignment/${sample}.NFRchrMfilter.sort.bam 
 

#By default, "-style histone" evokes a peak size of 500. 
#By default this is 2x the peak size, by default "-style histone" sets this to 1000

#findPeaks <tag directory> -style histone -o auto -i <input tag directory>

findPeaks TagDirectory/tag_${sample} -style histone -i TagDirectory/ > HomerPeaks/${sample}.peaks.txt

#######################################################################################
#convert to bed
	#remove header: 31lines
	sed -i '1,31d' HomerPeaksFDR01/${sample}.peaksFDR01.txt 
		
	#modified homer pos2bed.pl to  make findPeaks score in the forth column
	#print $filePtr "$chr\t$start\t$end\t$line[0]\t$v\t$dir\n"; to
	# print $filePtr "$chr\t$start\t$end\t$line[7]\t$v\t$dir\n"; 
	./pos2bedmod.pl ${sample}.peaksFDR01.txt > HomerPeaksFDR01/${sample}.peaksFDR01.bed

	sort -k1,1 -k2,2n HomerPeaksFDR01/${sample}.peaksFDR01.bed> HomerPeaksFDR01/${sample}.peaksFDR01.sort.bed
	
	#remove header: 30lines
	#sed -i '1,30d' HomerPeaksFDR01/${sample}.peaksFDR01.sort.bed
	
	#remove blacklist regions
	#bedtools	intersect	–v	–a	<peaks>	-b	<black>	>	<output>
	
	bedtools intersect -a HomerPeaksFDR01/${sample}.peaksFDR01.sort.bed -b deeptools/mm10.blacklist_sort.bed -v > HomerPeaksFDR01/${sample}.peaksFDR01.nobl.bed
	
	sort -k1,1 -k2,2n HomerPeaksFDR01/${sample}.peaksFDR01.nobl.bed > HomerPeaksFDR01/${sample}.peaksFDR01.nobl_sort.bed

	#bedToBigBed input.bed chrom.sizes myBigBed.bb
	bedToBigBed HomerPeaksFDR01/${sample}.peaksFDR01.nobl_sort.bed mm10.chrom.sizes HomerPeaksFDR01/${sample}.peaksFDR01.bb

#######################################################################################

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
