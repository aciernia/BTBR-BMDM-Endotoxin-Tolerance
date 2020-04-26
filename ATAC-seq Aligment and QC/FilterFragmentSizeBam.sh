#!/bin/bash

#https://www.biostars.org/p/310670/
#filters bam file by fragment size
#filters based on sam file field 9 = TLEN 
#number of bases from the leftmost mapped base to the rightmost mapped base.
#module load samtools/1.9

samtools view -h $1 | \
  awk -v LEN=$2 '{if ($9 <= LEN && $9 >= -(LEN) && $9 != 0 || $1 ~ /^@/) print $0}' | \
  samtools view -bh -o $1_147filter.bam

## Example to get fragments of 147bp or smaller:
#./code.sh in.bam 147