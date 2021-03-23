##make Tag directory for homer for A/B compartment calling

module load joaach/homer/4.7
module load qiadu/R/3.2.3
module load gi/samtools/1.2
module load gi/seqlogo/2.8.2


makeTagDirectory Dec4 /share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Decitabine/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Dec4/Hi-C_PDX_GAR15-13D_Dec4_R1_genome.bwt2merged.bam,/share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Decitabine/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Dec4/Hi-C_PDX_GAR15-13D_Dec4_R2_genome.bwt2merged.bam -illuminaPE

makeTagDirectory Dec5 /share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Decitabine/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Dec5/Hi-C_PDX_GAR15-13D_Dec5_R1_genome.bwt2merged.bam,/share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Decitabine/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Dec5/Hi-C_PDX_GAR15-13D_Dec5_R2_genome.bwt2merged.bam -illuminaPE

makeTagDirectory Dec9 /share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Decitabine/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Dec9/Hi-C_PDX_GAR15-13D_Dec9_R1_genome.bwt2merged.bam,/share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Decitabine/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Dec9/Hi-C_PDX_GAR15-13D_Dec9_R2_genome.bwt2merged.bam -illuminaPE

makeTagDirectory Veh1 /share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Vehicle/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Veh1/Hi-C_PDX_GAR15-13D_Veh1_R1_genome.bwt2merged.bam,/share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Vehicle/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Veh1/Hi-C_PDX_GAR15-13D_Veh1_R2_genome.bwt2merged.bam -illuminaPE

makeTagDirectory Veh2 /share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Vehicle/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Veh2/Hi-C_PDX_GAR15-13D_Veh2_R1_genome.bwt2merged.bam,/share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Vehicle/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Veh2/Hi-C_PDX_GAR15-13D_Veh2_R2_genome.bwt2merged.bam -illuminaPE

makeTagDirectory Veh4 /share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Vehicle/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Veh4/Hi-C_PDX_GAR15-13D_Veh4_R1_genome.bwt2merged.bam,/share/ScratchGeneral/joaach/Hi-C_Level2_Arima_PDX_GAR15-13D/hicpro_gar15-13d_Vehicle/bowtie_results/bwt2/Hi-C_PDX_GAR15-13D_Veh4/Hi-C_PDX_GAR15-13D_Veh4_R2_genome.bwt2merged.bam -illuminaPE

echo "done"
