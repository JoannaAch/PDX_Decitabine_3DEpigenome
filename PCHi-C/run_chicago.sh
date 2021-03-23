##run chicago (2020 updated R script)


module load centos7.8/joaach/R/3.6.3

Rscript runChicago_2020.R --design-dir DESIGN-DIR --en-feat-list feature_file --rda --export-format interBed,washU_text,washU_track /share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_PDX_GAR15-13D/chicago/PCHi-C_PDX_GAR15-13D_Dec4/PCHi-C_PDX_GAR15-13D_Dec4.chinput PCHi-C_PDX_GAR15-13D_Dec4_chicago
