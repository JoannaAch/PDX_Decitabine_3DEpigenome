module load gi/bowtie/2.2.4
module load fabbus/python/2.7.3
module load phuluu/R/3.1.2
module load gi/samtools/1.2
module load gi/gcc/4.8.2

time /home/joaach/HiC-Pro_2/bin/HiC-Pro_2.11.3 -c config_hic-pro_hic_pdx.txt -i fastq_Dec -o hicpro_gar15-13d_Decitabine

time /home/joaach/HiC-Pro_2/bin/HiC-Pro_2.11.3 -c config_hic-pro_hic_pdx.txt -i fastq_Veh -o hicpro_gar15-13d_Vehicle

echo "DONE"
