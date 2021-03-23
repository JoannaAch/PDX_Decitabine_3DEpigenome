##calculate number of "on target" reads per sample

echo "Dec_4"

###Step 1
module load gi/samtools/1.2 
samtools view PCHi-C_PDX_GAR15-13D_Dec4_R1_2.hicup.bam | wc -l



##Step 3
module load joaach/bedtools/2.25.0

bedtools intersect -u -bed -a PCHi-C_PDX_GAR15-13D_Dec4_R1_2.hicup.bam -b hg38_baitmap_final_sort.uniq.baitmap.bed > PCHi-C_PDX_GAR15-13D_Dec4_R1_2.hicup.bam.on-targeted.bed

##Step 4
cut -f4 PCHi-C_PDX_GAR15-13D_Dec4_R1_2.hicup.bam.on-targeted.bed | cut -f1 -d"/" | sort | uniq | wc -l 

echo "DONE"

