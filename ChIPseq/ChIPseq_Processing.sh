#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -N PDX_HCI005_ER_ChIPSeq_Processing
#$ -q short.q
#$ -m bea
#$ -M e.parbery@garvan.org.au

# this will print every command before running it,
# can help with debugging
set -x

# this will log the environment, useful for debugging!
(set -o posix; set; ulimit -a; uname -a; lsb_release -a; hostname -A) 1>&2

#load modules
module load elypar/bowtie2/2.2.9 
module load phuluu/samtools/1.11

#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc> | b <bam>} -S [<sam>]
cd /share/ScratchGeneral/elypar/PDX_HCI005_ER
for i in PDXHCI005_Veh_2_ER_ChIP PDXHCI005_Veh_3_ER_ChIP PDXHCI005_Veh_4_ER_ChIP PDXHCI005_Dec_pooled_input PDXHCI005_Veh_pooled_input
do
bowtie2 -x hg38 -U ./fastq/$i/$i.fq.gz -S ./bowtie2/$i.sam
samtools view -b ./bowtie2/$i.sam > ./bowtie2/$i.bam
samtools sort -o ./bowtie2/$i_sorted.bam -O bam ./bowtie2/$i.bam
samtools index -b ./bowtie2/$i_sorted.bam > ./bowtie2/$i.bai
echo $i Done
done

#run peak calling
cd /share/ScratchGeneral/elypar/PDX_HCI005_ER

for i in PDXHCI005_Dec_1_ER_ChIP PDXHCI005_Dec_2_ER_ChIP PDXHCI005_Dec_3_ER_ChIP PDXHCI005_Dec_pooled_input
do
macs2 callpeak -t ./bowtie2/$i_sorted.bam -c ./bowtie2/PDXHCI005_Dec_pooled_input_sorted.bam --outdir ./macs2 -n $i 2> ./macs2/$i.log 
echo $i Done
done

for i in PDXHCI005_Veh_1_ER_ChIP PDXHCI005_Veh_2_ER_ChIP PDXHCI005_Veh_3_ER_ChIP PDXHCI005_Veh_4_ER_ChIP PDXHCI005_Veh_pooled_input
do
macs2 callpeak -t ./bowtie2/$i_sorted.bam -c ./bowtie2/PDXHCI005_Veh_pooled_input_sorted.bam --outdir ./macs2 -n $i 2> ./macs2/$i.log 
echo $i Done
done

#make bigwig files with deeptools
module load elypar/deeptools/3.5.0
     
bamCoverage -b ./bowtie2/PDXHCI005_Veh_pooled_input_sorted.bam -o ./deeptools/PDXHCI005_Veh_pooled_input.bw

bamCoverage -b ./bowtie2/PDXHCI005_Veh_4_ER_ChIP_sorted.bam -o ./deeptools/PDXHCI005_Veh_4_ER_ChIP.bw

bamCoverage -b ./bowtie2/PDXHCI005_Veh_2_ER_ChIP_sorted.bam -o ./deeptools/PDXHCI005_Veh_2_ER_ChIP.bw

bamCoverage -b ./bowtie2/PDXHCI005_Veh_1_ER_ChIP_sorted.bam -o ./deeptools/PDXHCI005_Veh_1_ER_ChIP.bw

bamCoverage -b ./bowtie2/PDXHCI005_Dec_pooled_input_sorted.bam -o ./deeptools/PDXHCI005_Dec_pooled_input.bw

bamCoverage -b ./bowtie2/PDXHCI005_Dec_3_ER_ChIP_sorted.bam -o ./deeptools/PDXHCI005_Dec_3_ER_ChIP.bw

bamCoverage -b ./bowtie2/PDXHCI005_Dec_2_ER_ChIP_sorted.bam -o ./deeptools/PDXHCI005_Dec_2_ER_ChIP.bw

bamCoverage -b ./bowtie2/PDXHCI005_Dec_1_ER_ChIP_sorted.bam -o ./deeptools/PDXHCI005_Dec_1_ER_ChIP.bw
