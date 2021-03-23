module load centos7.8/joaach/python/2.7.8
module load centos7.8/joaach/TEtranscripts/2.2.1

TEtranscripts --format BAM --GTF Homo_sapiens.GRCh38.102.withchr_noheader5.gtf --TE GRCm38_Ensembl_rmsk_TE_withchr.gtf --mode multi --stranded reverse -t GAR15-13D_Dec_4_RNAAligned.out.bam GAR15-13D_Dec_5_RNAAligned.out.bam GAR15-13D_Dec_9_RNAAligned.out.bam GAR15-13D_Dec_10_RNAAligned.out.bam -c GAR15-13D_Vehicle_1_RNAAligned.out.bam GAR15-13D_Vehicle_2_RNAAligned.out.bam GAR15-13D_Vehicle_3_RNAAligned.out.bam GAR15-13D_Vehicle_4_RNAAligned.out.bam --project sample_nosort_Veh.vs.Dec

