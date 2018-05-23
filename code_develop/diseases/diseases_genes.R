
library(hipathia)

setwd("appl/hipathia/private/")

# Load annoation
gendis <- read.csv("diseases/CTDfiles/CTD_genes_diseases.csv", 
                    header = F, sep=",", stringsAsFactors = F, 
                    comment.char = "#")
colnames(gendis) <- c("GeneSymbol", "GeneID", "DiseaseName", "DiseaseID", 
                       "DirectEvidence", "InferenceChemicalName", 
                       "InferenceScore", "OmimsIDs", "PubMedIDs")
# Check dimension
length(unique(gendis$GeneSymbol))

# Load metaginfo & annotations
mgi <- hipathia:::load.mgi("hsa")
entrez2hgnc <- hipathia:::load.entrez.hgnc("hsa")
hgnc.path <- entrez2hgnc[entrez2hgnc$V1 %in% mgi$all.genes,]
length(intersect(hgnc.path$V2, unique(gendis$GeneSymbol)))
table(mgi$all.genes %in% entrez2hgnc$V1)

# Check GeneID == Entrez Gene
table(mgi$all.genes %in% gendis$GeneID)
gendis[which(gendis$GeneSymbol == "PPARGC1A")[1],]
gendis[which(gendis$GeneSymbol == "BRCA1")[1],]

# Filter pathway genes
gendis <- gendis[gendis$GeneID %in% mgi$all.genes,]

# Annotate paths
disannot <- unique(gendis[,c(1,3)])
colnames(disannot) <- c("gene", "function")
pathdis <- hipathia:::annotate.paths(metaginfo = mgi, dbannot = disannot)
annotated.dis <- unique(pathdis$funs)

# Save
save(gendis, file = "diseases/diseases_genes.RData")
save(pathdis, file = "diseases/annotations_disease_subpathways.RData")
save(annotated.dis, file = "diseases/list_of_annotated_diseases.RData")
