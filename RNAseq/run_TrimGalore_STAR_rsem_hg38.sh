#!/bin/bash

source /etc/profile.d/modules.sh
export MODULEPATH=/share/ClusterShare/Modules/modulefiles/noarch:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/contrib:$MODULEPATH

chrom_sizes_file="/share/ScratchGeneral/joaach/genomes/hg38/genome.chrom.sizes"

if [ ! -d "$1" ]; then
  echo "Directory does not exist" && exit
fi
directory=$1

if [ ! -f "$2" ]; then
  echo "Genome fasta file does not exist" && exit
fi

genome=$2

if [ ! -d "$3" ]; then
  echo "Genome index directory does not exist" && exit
fi

index=$3

if [ ! -f "$4" ]; then
  echo "GTF file does not exist" && exit
fi

gtf=$4

if [ ! -f ${5}".grp" ]; then
  echo "RSEM index (.grp) does not exist" && exit
fi

rsemindex=$5

if ! [ $6 -eq $6 ] 2> /dev/null
then
    echo "Read length must be an integer" && exit
fi

if [ -z $7 ]
then
  echo "Specify an output directory"
fi

mkdir -p $7

length=$( expr $6 - 1 )

echo $1 $2 $3 $4 $5 $6 $7

echo $1
R1s=$( ls $1/*_R1.fastq.gz | wc -l )
echo $R1s
R2s=$( ls $1/*_R2.fastq.gz | wc -l )
echo $R2s
if [ $R1s -ne $R2s ]; then
  echo "Missing mates for paired-end - please check all fastqs are present"
fi

samplename=$( echo $directory | sed s/.*\\///g )
echo $samplename


module load gi/trim_galore/0.4.0
module load gi/fastqc/0.11.2

mkdir -p $7/trimgalore
CMD="trim_galore --paired --fastqc_args \"--nogroup\" --output_dir ${7}/trimgalore ${1}/${samplename}\"_R1.fastq.gz\" ${1}/${samplename}\"_R2.fastq.gz\""
echo $CMD  && eval $CMD

module load nenbar/star/2.4.0d
module load gi/samtools/1.1
module load gi/novosort/1.02.02
module load gi/bedtools/2.22.0
module load gi/ucsc_utils/283

mkdir -p $7/star_TMP
mkdir -p $7/star

# Build STAR index
# check assets/build_STAR_index.sh for required annotation/index

# Argument --outFilterMatchNmin will have to be tailored to <read length> - 1

CMD="STAR --runMode alignReads --genomeDir $3 --sjdbGTFfile $4 --readFilesIn ${7}/trimgalore/${samplename}\"_R1_val_1.fq.gz\" ${7}/trimgalore/${samplename}\"_R2_val_2.fq.gz\" --runThreadN 6 --genomeLoad LoadAndKeep --readFilesCommand zcat --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMax 1500000 --alignMatesGapMax 1500000 --alignIntronMin 20 --alignSJoverhangMin 6 --alignSJDBoverhangMin 1 --outFilterMatchNmin $length --quantMode TranscriptomeSAM --outFileNamePrefix ${7}/star_TMP/${samplename}_ --outSAMtype BAM Unsorted"
echo $CMD && eval $CMD

# mv STAR alignment details to main STAR folder
mv ${7}/star_TMP/${samplename}_Log.final.out ${7}/star/${samplename}.log

# -f 3 means only output alignments which are read paired and read mapped in proper pair.
CMD="samtools view -@ 6 ${7}/star_TMP/${samplename}_Aligned.toTranscriptome.out.bam -f 3 -b > ${7}/star/${samplename}\"_out.bam\""
echo $CMD && eval $CMD

CMD="novosort -n -m 16G -c 6 ${7}/star/${samplename}\"_out.bam\" > ${7}/star/${samplename}\".sorted.transcriptome.bam\""
echo $CMD && eval $CMD

rm ${7}/star/${samplename}"_out.bam"

mkdir -p $7/star_TMP/tmp/
# sorted genome bam file by coordinate
CMD="samtools sort -T $7/star_TMP/tmp/${samplename}_Aligned.out.sorted -o $7/star_TMP/${samplename}_Aligned.out.sorted.bam $7/star_TMP/${samplename}_Aligned.out.bam"
echo $CMD && eval $CMD

# index sorted bam file
CMD="samtools index $7/star_TMP/${samplename}_Aligned.out.sorted.bam"
echo $CMD && eval $CMD

mkdir -p $7/genome
# move genome files to genome/ folder
mv $7/star_TMP/${samplename}_Aligned.out.sorted.bam $7/genome/${samplename}.genome.sorted.bam
mv $7/star_TMP/${samplename}_Aligned.out.sorted.bam.bai $7/genome/${samplename}.genome.sorted.bam.bai

# generate bigwig file for genome files
CMD="genomeCoverageBed -bg -split -ibam $7/genome/${samplename}.genome.sorted.bam -g ${chrom_sizes_file} > $7/genome/${samplename}.bedGraph"
echo $CMD && eval $CMD
CMD="sed -n '/chr/p' $7/genome/${samplename}.bedGraph > $7/genome/${samplename}.bedGraph.filt"
echo $CMD && eval $CMD
CMD="sort -k1,1 -k2,2n $7/genome/${samplename}.bedGraph.filt > $7/genome/${samplename}.bedGraph.filt.sorted"
echo $CMD && eval $CMD
CMD="bedGraphToBigWig $7/genome/${samplename}.bedGraph.filt.sorted ${chrom_sizes_file} $7/genome/${samplename}.genome.bw"
echo $CMD && eval $CMD

# tidy up and remove intermediary files
rm $7/genome/${samplename}.bedGraph
rm $7/genome/${samplename}.bedGraph.filt
rm $7/genome/${samplename}.bedGraph.filt.sorted

module load gi/rsem/1.2.21

mkdir -p $7/rsem

# Build RSEM index
# check assets/build_RSEM_index.sh for required annotation/index

CMD="rsem-calculate-expression --paired-end --bam --no-bam-output --seed 12345 -p 6 --forward-prob 0 ${7}/star/${samplename}\".sorted.transcriptome.bam\" $5 ${7}/rsem/${samplename}"
echo $CMD && eval $CMD
