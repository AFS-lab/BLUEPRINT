library(DESeq2)
library(countsimQC)

load("../data/SupplementaryFigure8.RData")

countsimQCReport(ddsList = countsimExample, 
                 outputFile = "SupplementaryFigure8.html", 
                 outputDir = "../pdfs", 
                 savePlots = TRUE,
                 description = "This is a comparison of three count data sets.")

ggplots <- readRDS("../pdfs/SupplementaryFigure8_ggplots.rds")
if (!dir.exists("../pdfs/figures")) dir.create("../pdfs/figures")
generateIndividualPlots(ggplots, device = "pdf", nDatasets = 3, 
                        outputDir = "../pdfs/figures")