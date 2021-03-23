##Create .chinput files for CHiCAGO from hicup bam files

#!/bin/sh

source /etc/profile.d/modules.sh
export MODULEPATH=/share/ClusterShare/Modules/modulefiles/noarch:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/contrib:$MODULEPATH

module load joaach/bedtools/2.25.0
module load fabbus/perl/5.14.2

./bam2chicago.sh /share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_PDX_GAR15-13D/hicup/PCHi-C_PDX_GAR15-13D_Dec4_R1_2.hicup.bam hg38_baitmap_final_sort.uniq.5col.baitmap hg38_GATC_GANTC_fragments_final.rmap PCHi-C_PDX_GAR15-13D_Dec4 [nodelete]

