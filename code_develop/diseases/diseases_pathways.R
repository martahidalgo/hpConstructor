
library(hipathia)

setwd("appl/hipathia/private/")

pathwaydis <- read.csv("diseases/CTDfiles/CTD_diseases_pathways.csv", 
                     header = F, sep=",", stringsAsFactors = F, 
                    comment.char = "#")
colnames(pathwaydis) <- c("DiseaseName", "DiseaseID", "PathwayName", 
                       "PathwayID", "InferenceGeneSymbol")

# Filter Reactome pathways
pathwaydis <- pathwaydis[grepl("KEGG:", pathwaydis$PathwayID),]
# Filter modules
pathwaydis <- pathwaydis[!grepl("_M", pathwaydis$PathwayID),]
# Simplify Pathway ID
pathwaydis$PathwayID_simple <- gsub("KEGG:", "", pathwaydis$PathwayID)
pathwaydis$PathwayID_num <- gsub("KEGG:hsa", "", pathwaydis$PathwayID)

save(pathwaydis, file = "diseases/diseases_pathways.RData")


