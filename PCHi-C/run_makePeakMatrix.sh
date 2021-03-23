##make peak matrix for use with chicdiff (Cairns Bioinformatics 2019)

module load centos7.8/joaach/R/3.6.3


Rscript makePeakMatrix.R --twopass --notrans --rda --var cd --scorecol score names-file --cutoff 5 peakMatrix_PDX_Veh.vs.Dec.txt
