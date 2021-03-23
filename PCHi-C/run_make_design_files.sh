##Create design files from baitmap and rmap for use with CHiCAGO

module load phuluu/python/2.7.8

python makeDesignFiles.py --rmapfile=hg38_GATC_GANTC_fragments_final.rmap --baitmapfile=hg38_baitmap_final_sort.uniq.5col.baitmap --outfilePrefix=hg38_GATC_GANTC_fragments --minFragLen=25 --maxFragLen=1200 --maxLBrownEst=75000 --binsize=1500 --removeb2b=True --removeAdjacent=True
