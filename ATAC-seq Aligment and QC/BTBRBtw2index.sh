#!/bin/bash
#
#SBATCH --workdir /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/BTBR
#SBATCH -p production,assembly              
#SBATCH -c 12   
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=1-00:00:00

#from screen
#kinit -l 22d
#klist
#aklog
#sbatch --exclude cycle-0 /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/BTBR_genome/mmarge_genome.sh
#mkdir -p bowtie2/index
module load bowtie2/2.3.4.1

#concatenate all fasta files into one large file
cat /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/BTBR/*.fa >> /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/BTBR/genome-btbr.fa

#generate a bowtie2 index,
bowtie2-build /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/BTBR/genome-btbr.fa /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/bowtie2/index/btbr

#move these files into a general bowtie2 index folder
#mv *bt2 /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/bowtie2/index