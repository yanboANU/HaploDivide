#!/bin/bash

/home/admin-u6260133/Tools/bwa-0.7.17/bwa index $1

/home/admin-u6260133/Tools/bwa-0.7.17/bwa mem -t 16 $1 ../../reads1_2.fastq >yeast1.mem.sam

python3 /home/admin-u6260133/script/HaploDivide/read_sam.py yeast1.mem.sam 1_filter.sam >1_filter.log

/home/admin-u6260133/Tools/samtools/samtools/samtools view -@ 16 -S -b 1_filter.sam >1_filter.bam

/home/admin-u6260133/Tools/samtools/samtools/samtools sort -@ 16 1_filter.bam -o 1_filter.sorted.bam

/home/admin-u6260133/Tools/samtools/samtools/samtools index -@ 16 1_filter.sorted.bam

mkdir after_filter

cd after_filter

python3 /home/admin-u6260133/script/HaploDivide/filter_bam.py ../1_filter.sorted.bam

/home/admin-u6260133/Tools/samtools/samtools/samtools index -@ 16 2_filter.bam

python3 /home/admin-u6260133/script/HaploDivide/version2/main.py  2_filter.bam ../$1 >log
