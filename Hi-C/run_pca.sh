##Run PCA analyses in Homer for A/B compartment calling

module load joaach/homer/4.7
module load qiadu/R/3.2.3
module load gi/samtools/1.2
module load gi/seqlogo/2.8.2

runHiCpca.pl pcaOut_Dec4 Dec4/ -res 100000 -cpu 16 -genome hg38
runHiCpca.pl pcaOut_Dec5 Dec5/ -res 100000 -cpu 16 -genome hg38
runHiCpca.pl pcaOut_Dec9 Dec9/ -res 100000 -cpu 16 -genome hg38


runHiCpca.pl pcaOut_Veh1 Veh1/ -res 100000 -cpu 16 -genome hg38
runHiCpca.pl pcaOut_Veh2 Veh2/ -res 100000 -cpu 16 -genome hg38
runHiCpca.pl pcaOut_Veh4 Veh4/ -res 100000 -cpu 16 -genome hg38

echo "done"
