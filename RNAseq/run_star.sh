##run STAR for TEtranscripts analyses

module load centos7.8/joaach/STAR/2.7.7a

STAR --runMode alignReads --genomeDir /share/ScratchGeneral/joaach/Repeats/Data/star/STAR_index_hg38 --sjdbGTFfile /share/ScratchGeneral/joaach/genomes/hg38_GTF/Homo_sapiens.GRCh38.102.withchr_noheader5.gtf --readFilesCommand zcat --readFilesIn /share/ScratchGeneral/joaach/Repeats/Data/star/fastq/GAR15-13D_Vehicle_1_RNA_R1.fastq.gz /share/ScratchGeneral/joaach/Repeats/Data/star/fastq/GAR15-13D_Vehicle_1_RNA_R2.fastq.gz --runThreadN 16 --outFilterMultimapNmax 100 --outFilterMismatchNmax 100 --outFileNamePrefix GAR15-13D_Vehicle_1_RNA

module load centos6.10/gi/samtools/1.2 
samtools view -S -b  GAR15-13D_Vehicle_1_RNAAligned.out.sam >  GAR15-13D_Vehicle_1_RNAAligned.out.bam
