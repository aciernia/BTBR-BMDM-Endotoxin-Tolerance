#!/bin/bash

#SBATCH --array=1-12  # NEED TO CHANGE THIS!
#SBATCH --job-name=BTBRalignment # Job name
#SBATCH --workdir /share/lasallelab/Annie/BMDM/ATACanalysis2019
#SBATCH -p production              
#SBATCH -c 12     
#SBATCH --mem-per-cpu=8000 
#SBATCH --time=3-00:00:00
#SBATCH --output=arrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=arrayJob_%A_%a.err # File to which STDERR will be written
#######################################################################################

begin=`date +%s`

echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" BTBRsamples.txt`

#alignment for all BTBR samples
#######################################################################################
#module load fastqc/0.11.7
#module load fastq_screen/0.11.4
module load bowtie2/2.3.4.1
module load samtools/1.8
module load bedtools2/2.27.0
module load picard-tools/2.18.4
module load java
module load mmarge/1.0
module load homer/4.9
#######################################################################################
#trim and fastqc
#######################################################################################

#fastqc
#fastqc raw_sequences/${sample}*R1*.fastq.gz --outdir fastqc_pretrim
#fastqc raw_sequences/${sample}*R2*.fastq.gz --outdir fastqc_pretrim

#Trimmomatic for PE
#trim adapters (TruSeq3-PE.fa)
# TruSeq3-PE.fa must be in the same folder as the fastq files


#PE trimming for adapters and quality
#java -jar /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/mmarge/Trimmomatic-0.38/trimmomatic-0.38.jar PE raw_sequences/${sample}*R1*.fastq.gz raw_sequences/${sample}*R2*.fastq.gz trimmed/${sample}_1.paired.fastq.gz trimmed/${sample}_1.unpaired.fastq.gz trimmed/${sample}_2.paired.fastq.gz trimmed/${sample}_2.unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:T LEADING:3 TRAILING:3 MINLEN:15 &> trimlogs/trim_log_${sample}

	#trimmomatic will look for seed matches of 16 bases with 2 mismatches allowed
	#will then extend and clip if a score of 30 for PE or 10 for SE is reached (~17 base match)
	#minimum adapter length is 8 bases
	#T = keeps both reads even if only one passes criteria
	#trims low quality bases at leading and trailing end if quality score < 15
	#sliding window: scans in a 4 base window, cuts when the average quality drops below 15
	#log outputs number of input reads, trimmed, and surviving

#fastqc
#fastqc trimmed/${sample}_1.paired.fastq.gz --outdir fastqc_posttrim
#fastqc trimmed/${sample}_2.paired.fastq.gz --outdir fastqc_posttrim

###############################################################################
#align trimmed reads using bowtie2 to BTBR genome made with mmarge
# --no-unal          suppress SAM records for unaligned reads

echo "starting sample ${sample}"

bowtie2 -x /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/bowtie2/index/btbr -1 /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/mmarge/trimmed/${sample}_1.paired.fastq.gz -2 /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/mmarge/trimmed/${sample}_2.paired.fastq.gz -S alignment/${sample}.sam -p 12 --no-unal --time 2>alignment/bowtielogs/${sample}_bowtie2log.log  

#######################################################################################
#Shift data to reference coordinates
#shift to mm10 C57 coordinates
#######################################################################################
#######################################################################################
echo "shifting sample ${sample}"


MMARGE.pl shift -paired -files alignment/${sample}.sam -ind btbr -data_dir /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/

# Will produce name sorted BAM
#samtools view -bS alignment/${sample}_shifted_from_BTBR.sam  | samtools sort -n -o alignment/${sample}_shifted_from_BTBR_sort.bam 


#######################################################################################
#sort with samtools and filter
#######################################################################################

#properly mapped and paired reads:
#-u Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command.
#only output properly paired reads -f 0x2
echo "samtools filtering sample ${sample}"

samtools view -bS -f 0x2 -u alignment/${sample}_shifted_from_BTBR.sam  | samtools sort -n -o alignment/${sample}_paired.bam # Will produce name sorted BAM

#fix mates and remove 2ndary and unmapped reads:
#needs name sorted bam as input
#need -m for markdup
#-r removes 2ndary and unmapped reads
samtools fixmate -m -r alignment/${sample}_paired.bam alignment/${sample}_fixmate.bam

#coordinate sort:
samtools sort alignment/${sample}_fixmate.bam -o alignment/${sample}_fixmate_sort.bam

#build bam index file
picard BuildBamIndex INPUT=alignment/${sample}_fixmate_sort.bam OUTPUT=alignment/${sample}_fixmate_sort.bam.bai

#picard alignment summary
picard CollectAlignmentSummaryMetrics INPUT=alignment/${sample}_fixmate_sort.bam O=alignment/${sample}_fixmate_alignsum.tsv

#flagstat
samtools flagstat alignment/${sample}_fixmate_sort.bam  > alignment/${sample}_fixmateFlagstat.txt

echo "samtools removing duplicates ${sample}"
#mark duplicates and remove
samtools markdup -r -s alignment/${sample}_fixmate_sort.bam alignment/${sample}_dedup.bam

#fix mate information
picard FixMateInformation I=alignment/${sample}_dedup.bam

#build bam index file
picard BuildBamIndex INPUT=alignment/${sample}_dedup.bam OUTPUT=alignment/${sample}_dedup.bam.bai

#picard alignment summary
picard CollectAlignmentSummaryMetrics INPUT=alignment/${sample}_dedup.bam O=alignment/${sample}_dedup_alignsum.tsv

#find insert size
picard CollectInsertSizeMetrics I=alignment/${sample}_dedup.bam H=alignment/${sample}_histogram.pdf O=alignment/${sample}_dedupinsertmetric.tsv

#flagstat
samtools flagstat alignment/${sample}_dedup.bam > alignment/${sample}_dedupFlagstat.txt


#######################################################################################
#Shift filtered bam file for Tn5 insertions
#Arguments required: 1. sam or bam alignment file (must end in .bam or .sam) 2. Filehandle for the output
echo "shifting sample ${sample}"
perl ATAC_BAM_shifter_gappedAlign.pl alignment/${sample}_dedup.bam alignment/${sample}_dedupTn5shift

#build bam index file
picard BuildBamIndex INPUT=alignment/${sample}_dedupTn5shift.bam OUTPUT=alignment/${sample}_dedupTn5shift.bam.bai

#######################################################################################
#filter bam files for regions of open chromatin < 147 base pairs
#use bash script FilterFragmentSizeBam.sh made from https://www.biostars.org/p/310670/ 
#filters by fragment size using samtools field TLEN and awk
#requires samtools
#Arguments required: bam file and bp filter
#outputs filtered bam file with this name: input.bam_147filter.bam

echo "filtering sample ${sample}"
./FilterFragmentSizeBam.sh alignment/${sample}_dedupTn5shift.bam 147

#build bam index file
picard BuildBamIndex INPUT=alignment/${sample}_dedupTn5shift.bam_147filter.bam OUTPUT=alignment/${sample}_dedupTn5shift.bam_147filter.bam.bai

# removed chrM reads and sort by reads
samtools view -h alignment/${sample}_dedupTn5shift.bam_147filter.bam | awk '{if($3 != "chrM"){print $0}}' | samtools view -Sb | samtools sort -n > alignment/${sample}.NFRchrMfilter.bam

#######################################################################################
#sort and build bam index file
samtools view -bS alignment/${sample}.NFRchrMfilter.bam | samtools sort -o alignment/${sample}.NFRchrMfilter.sort.bam


picard BuildBamIndex INPUT=alignment/${sample}.NFRchrMfilter.sort.bam OUTPUT=alignment/${sample}.NFRchrMfilter.sort.bam.bai

#find insert size
picard CollectInsertSizeMetrics I=alignment/${sample}.NFRchrMfilter.sort.bam H=alignment/${sample}_NFRhistogram.pdf O=alignment/${sample}_NFRinsertmetric.tsv

#flagstat
samtools flagstat alignment/${sample}.NFRchrMfilter.sort.bam > alignment/${sample}_NFRflagstat.txt


#######################################################################################
echo "making tags ${sample}"
#make tag directories for each bam file
makeTagDirectory TagDirectory/tag_${sample} alignment/${sample}.NFRchrMfilter.sort.bam 
 
#call peaks for ATAC data (from Verena settings Link 2018 paper)
#findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
#If "-o auto" is specified, the peaks will be written to:
#"<tag directory>/peaks.txt" (-style factor)
#	$command = "findPeaks " . $output_dir . "/pooled_tag_dirs -L 0 -C 0 -fdr 0.9 -minDist 200 -size 200 > " . $output_dir . "/peaks_pooled_tag_dirs.txt 2>> $output_dir/error.txt";

#run with super liberal and then filter by IDR
#-L <#> (fold enrichment over local tag count, default: 4.0)
# -C <#> (fold enrichment limit of expected unique tag positions, default: 2.0)
#-minDist <#> (minimum distance between peaks, default: peak size x2)
# didn't manually set size to 200, left as auto
echo "finding peaks ${sample}"
findPeaks TagDirectory/tag_${sample} -style factor -center -F 0 -L 0 -C 0 -fdr 0.01 -minDist 200 -size 200 -o auto 

#######################################################################################

#######################################################################################

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
