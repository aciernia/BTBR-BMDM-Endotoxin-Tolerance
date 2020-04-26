#!/bin/bash

#SBATCH --array=1-12  # NEED TO CHANGE THIS!
#SBATCH --job-name=C57alignment # Job name
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

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" C57samples.txt`

#alignment for all C57 samples
#######################################################################################
#module load fastqc/0.11.7
#module load fastq_screen/0.11.4
module load bowtie2/2.3.4.1
module load samtools/1.8
module load bedtools2/2.27.0
module load picard-tools/2.18.4
module load java
module load mmarge/1.0
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
#align trimmed reads using bowtie2
# --no-unal          suppress SAM records for unaligned reads
bowtie2 -x /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/bowtie2/index/c57bl6j -1 /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/mmarge/trimmed/${sample}_1.paired.fastq.gz -2 /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/mmarge/trimmed/${sample}_2.paired.fastq.gz -S alignment/${sample}.sam -p 12 --no-unal --time 2>alignment/bowtielogs/${sample}_bowtie2log.log  

#######################################################################################
#sort with samtools and filter
#######################################################################################

#properly mapped and paired reads:
#-u Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command.
#only output properly paired reads -f 0x2

samtools view -bS -f 0x2 -u alignment/${sample}.sam | samtools sort -n -o alignment/${sample}_paired.bam # Will produce name sorted BAM

#fix mates and remove 2ndary and unmapped reads:
#needs name sorted bam as input
#need -m for markdup
#-r removes 2ndary and unmapped reads
samtools fixmate -m -r alignment/${sample}_paired.bam alignment/${sample}_fixmate.bam

#coordinate sort:
samtools sort alignment/${sample}_fixmate.bam -o alignment/${sample}_fixmate_sort.bam

#mark duplicates and remove
samtools markdup -r -s alignment/${sample}_fixmate_sort.bam alignment/${sample}_dedup.bam

#fix mate information
picard FixMateInformation I=alignment/${sample}_dedup.bam

#build bam index file
picard BuildBamIndex INPUT=alignment/${sample}_dedup.bam OUTPUT=alignment/${sample}_dedup.bam.bai

#picard alignment summary
picard CollectAlignmentSummaryMetrics INPUT=alignment/${sample}_dedup.bam O=alignment/${sample}_dedup_alignsum.tsv

#find insert size
picard CollectInsertSizeMetrics I=alignment/${sample}_dedup.bam H=alignment/${sample}_histogram.pdf O=alignment/${sample}_insertmetric.tsv

#flagstat
samtools flagstat alignment/${sample}_dedup.bam > alignment/${sample}_dedupFlagstat.txt



#######################################################################################

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
