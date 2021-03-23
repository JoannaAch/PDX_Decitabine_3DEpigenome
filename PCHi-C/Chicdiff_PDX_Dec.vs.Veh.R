library(Chicdiff)


setwd("/share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_PDX_GAR15-13D/chicago/chicdiff/")
dataPath <- ("/share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_PDX_GAR15-13D/chicago/chicdiff/input_rda")
testDesignDir <- file.path(dataPath, "DESIGN-DIR/")
dir(testDesignDir)




peakFiles <- file.path(dataPath,"peakMatrix_PDX_Veh.vs.Dec_final.txt.txt")


inputdataPath <- ("/share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_PDX_GAR15-13D/chicago/chicdiff/")
testDataPath_Dec <- file.path(inputdataPath, "Dec_chinputfiles")
testDataPath_Veh <- file.path(inputdataPath, "Veh_chinputfiles")
dir(testDataPath_Dec)
dir(testDataPath_Veh)



countData <- list(
  Veh = c(Veh1 = file.path(testDataPath_Veh, "PCHi-C_PDX_GAR15-13D_Veh1.chinput"),
          Veh2 = file.path(testDataPath_Veh, "PCHi-C_PDX_GAR15-13D_Veh2.chinput"),
          Veh4 = file.path(testDataPath_Veh, "PCHi-C_PDX_GAR15-13D_Veh4.chinput")),
  
  
  Dec = c(Dec4 = file.path(testDataPath_Dec, "PCHi-C_PDX_GAR15-13D_Dec4.chinput"),
          Dec5 = file.path(testDataPath_Dec, "PCHi-C_PDX_GAR15-13D_Dec5.chinput"),
          Dec9 = file.path(testDataPath_Dec, "PCHi-C_PDX_GAR15-13D_Dec9.chinput"))
)





testDataPath_rds <- ("/share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_PDX_GAR15-13D/chicago/chicdiff/input_rda/")
dir(testDataPath_rds)

chicagoData <- list(
  Veh = c(Veh1 = file.path(testDataPath_rds, "PCHi-C_PDX_GAR15-13D_Veh1_chicago.RDa"),
          Veh2 = file.path(testDataPath_rds, "PCHi-C_PDX_GAR15-13D_Veh2_chicago.RDa"),
          Veh4 = file.path(testDataPath_rds, "PCHi-C_PDX_GAR15-13D_Veh4_chicago.RDa")),
  Dec = c(Dec4 = file.path(testDataPath_rds, "PCHi-C_PDX_GAR15-13D_Dec4_chicago.RDa"),
          Dec5 = file.path(testDataPath_rds, "PCHi-C_PDX_GAR15-13D_Dec5_chicago.RDa"),
          Dec9 = file.path(testDataPath_rds, "PCHi-C_PDX_GAR15-13D_Dec9_chicago.RDa"))
)



###RUN


chicdiff.settings <- setChicdiffExperiment(designDir = testDesignDir, chicagoData = chicagoData, countData = countData, peakfiles = peakFiles,  outprefix="test")

output <- chicdiffPipeline(chicdiff.settings)

save.image("/share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_CTCFKD/chicago/bams/chicdiff/chicdiff_Veh.vs.Dec_output.RData")

outCI <- getCandidateInteractions(chicdiff.settings = chicdiff.settings, 
                                  output = output, peakFiles = peakFiles, 
                                  pvcut = 1, minDeltaAsinhScore = 1)


write.table(outCI, "outCI_Veh.vs.Dec_all.txt", quote = FALSE, sep="\t", row.names = TRUE, col.names = TRUE)


outCI <- getCandidateInteractions(chicdiff.settings = chicdiff.settings, 
                                  output = output, peakFiles = peakFiles, 
                                  pvcut = 0.1, minDeltaAsinhScore = 1)


write.table(outCI, "outCI_Veh.vs.Dec_FDR510.txt", quote = FALSE, sep="\t", row.names = TRUE, col.names = TRUE)


save.image("/share/ScratchGeneral/joaach/PCHi-C_Level2_Arima_CTCFKD/chicago/bams/chicdiff/chicdiff_Veh.vs.Dec.RData")


quit()
n
