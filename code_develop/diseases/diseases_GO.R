library(hipathia)

setwd("appl/hipathia/private/")

# Load Go annotations
godis <- read.csv("diseases/CTDfiles/CTD_Disease-GO_biological_process_associations.csv", 
                  header = F, sep=",", stringsAsFactors = F, comment.char = "#")
colnames(godis) <- c("DiseaseName", "DiseaseID", "GOName", "GOID", 
                     "InferenceGeneQty", "InferenceGeneSymbol")

# Load GO annotations
annogos <- hipathia:::load.annofuns("GO", "hsa")

# Check number of intersections
length(unique(intersect(godis$GOName, annogos$funs))) # 1607
length(unique(annogos$funs)) # 1655
length(unique(godis$GOName)) # 12059
length(unique(godis$DiseaseName)) # 5702

# Common annotated GOs
comms <- unique(intersect(godis$GOName, annogos$funs))
# Diseases in the annotated GOs
diss <- unique(godis[godis$GOName %in% comms, "DiseaseName"])

# Filter godis by pathway annotated GOs
godis <- godis[godis$GOName %in% comms, ]

save(godis, file = "diseases/GO_diseases.RData")


