#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -N prognosis_ER_ChIPSeq_processing
#$ -q short.q
#$ -m bea
#$ -M e.parbery@garvan.org.au


# the script will now error out if something goes wrong
# instead of trying to continue
# it is not foolproof, subshells and various things do not behave as expected
set -euo pipefail

# this will print every command before running it,
# can help with debugging
set -x

# this will log the environment, useful for debugging!
(set -o posix; set; ulimit -a; uname -a; lsb_release -a; hostname -A) 1>&2

#load modules
module load elypar/bowtie2/2.2.9 
module load phuluu/samtools/1.11

WORKDIR=/share/ScratchGeneral/elypar/prognosis_ER

#read config file containing base sample names
file_list=$WORKDIR/bowtie2.config

# Loop over each sample listed in bowtie.config file
while read SAMPLE; do
   INPUT_1=${SAMPLE}
    echo "Input files:" $INPUT_1
    bowtie2 -x hg38 -U /share/ScratchGeneral/elypar/GSE32222/$INPUT_1.fastq.gz |samtools view -bS > /share/ScratchGeneral/elypar/prognosis_ER/bowtie2/$INPUT_1.bam
	samtools sort -o /share/ScratchGeneral/elypar/prognosis_ER/bowtie2/${INPUT_1}_sorted.bam -O bam /share/ScratchGeneral/elypar/prognosis_ER/bowtie2/$INPUT_1.bam
	samtools index -b /share/ScratchGeneral/elypar/prognosis_ER/bowtie2/${INPUT_1}_sorted.bam 
	echo $INPUT_1 Done
done < ${file_list}

#run peak calling
cd /share/ScratchGeneral/elypar/prognosis_ER

module load phuluu/python/3.9.1


WORKDIR=/share/ScratchGeneral/elypar/prognosis_ER

#read config file containing base sample names
sample_list=$WORKDIR/macs3_test.config

while read SAMPLE; do
   test=${SAMPLE}
    echo "Input files:" $test
    macs3 callpeak -t ./bowtie2/${test}_sorted.bam -c ./bowtie2/${test}_input.bam --outdir ./macs2 -n $test
	echo $test Done
done < ${sample_list}

#Make bigwig files with deeptools
module load elypar/deeptools/3.5.0 

WORKDIR=/share/ScratchGeneral/elypar/prognosis_ER

#read config file containing base sample names
sample_list=$WORKDIR/macs3_test.config

# Loop over each sample listed in bowtie.config file
while read SAMPLE; do
   test=${SAMPLE}
    echo "Input files:" $test
    bamCoverage -b ./bowtie2/${test}_sorted.bam -o ./deeptools/$test.bw
	echo $test Done
done < ${sample_list}
