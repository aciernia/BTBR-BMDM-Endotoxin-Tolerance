#!/bin/bash
#


#md5sum -c master_MD5_sam.txt > MD5_samchecksum.txt

#loop through folders and check checksum against md5sum.txt

for f in */;
do
echo "$f"
cd $f
md5sum -c md5sum.txt > MD5_checksum.txt
cd ..

done;

cat */MD5_checksum.txt > masterMD5check.txt

