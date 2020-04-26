#!/bin/bash
#
#SBATCH --workdir /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/REFERENCE
#SBATCH -p production,assembly              
#SBATCH -c 12   
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=1-00:00:00

#from screen
#kinit -l 22d
#klist
#aklog

module load bowtie2/2.3.4.1

#concatenate all fasta files into one large file
cat /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/REFERENCE/*.fa >> /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/REFERENCE/genome-c57.fa

#generate a bowtie2 index,
bowtie2-build /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/REFERENCE/genome-c57.fa /share/lasallelab/Annie/BMDM/ATACanalysis2019/genomes_annie/bowtie2/index/c57bl6j

#move these files into a general bowtie2 index folder
#mv /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/BTBR_genome/C57BL6J/*bt2 /share/lasallelab/Annie/BMDM/ATACanalysis_4_2018/BTBR_genome/bowtie2/index
