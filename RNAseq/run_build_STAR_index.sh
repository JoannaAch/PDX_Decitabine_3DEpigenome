# build STAR index
module load centos7.8/joaach/STAR/2.7.7a

read_length=150
length=$(expr $read_length - 1)

STAR_DIR="STAR_index_hg38/"
mkdir -p ${STAR_DIR}

echo ${STAR_DIR}

STAR_CMD="STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ${STAR_DIR} --genomeFastaFiles /share/ScratchGeneral/joaach/genomes/hg38/genome.fa --sjdbGTFfile /share/ScratchGeneral/joaach/genomes/hg38_GTF/Homo_sapiens.GRCh38.102.withchr_noheader5.gtf --sjdbOverhang ${length}"
echo $STAR_CMD && eval $STAR_CMD

