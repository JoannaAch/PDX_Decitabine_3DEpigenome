##call TADs using tadtool (Kruse, Bioinformatics, 2016) with pre-defined parameters (best parameters identifed with tadtool subset to chr1 and tadtool plot for 3-5 different regions)

module load joaach/python/2.7.8

tadtool tads Hi-C_PDX_GAR15-13D_Dec4_50000_iced.matrix Hi-C_PDX_GAR15-13D_Dec4_50000_abs.bed 102353 50 -v > Dec4_values_50kb.txt
tadtool tads Hi-C_PDX_GAR15-13D_Dec4_50000_iced.matrix Hi-C_PDX_GAR15-13D_Dec4_50000_abs.bed 102353 50 > Dec4_TADs_50kb.bed
tadtool tads Hi-C_PDX_GAR15-13D_Dec5_50000_iced.matrix Hi-C_PDX_GAR15-13D_Dec5_50000_abs.bed 102353 50 -v > Dec5_values_50kb.txt
tadtool tads Hi-C_PDX_GAR15-13D_Dec5_50000_iced.matrix Hi-C_PDX_GAR15-13D_Dec5_50000_abs.bed 102353 50 > Dec5_TADs_50kb.bed
tadtool tads Hi-C_PDX_GAR15-13D_Dec9_50000_iced.matrix Hi-C_PDX_GAR15-13D_Dec9_50000_abs.bed 102353 50 -v > Dec9_values_50kb.txt
tadtool tads Hi-C_PDX_GAR15-13D_Dec9_50000_iced.matrix Hi-C_PDX_GAR15-13D_Dec9_50000_abs.bed 102353 50 > Dec9_TADs_50kb.bed

echo "done"
