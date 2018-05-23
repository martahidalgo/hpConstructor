
library(hipathia)
library(devtools)
library(SummarizedExperiment)

load("private/example_data/brca_data.RData")
load("private/example_data/brca_design.RData")

brca <- SummarizedExperiment(assays=SimpleList(raw=brca_data), 
                             colData=brca_design)

# load and filter graphs
pathways <- load_pathways(species = "hsa",
                          pathways_list = c("hsa03320", "hsa04012"))

# prepare data
trans_data <- translate_data(brca_data, "hsa")
exp_data <- normalize_data(trans_data)

results <- hipathia(exp_data, pathways, verbose=TRUE)

path_vals <- get_paths_matrix(results)
go_vals <- quantify_terms(results, pathways, "GO")
uni_vals <- quantify_terms(results, pathways, "uniprot")

sample_group <- brca_design[colnames(path_vals),"group"]
comp <- do_wilcoxon(path_vals, sample_group, g1 = "Tumor", g2 = "Normal")
path_names <- get_path_names(pathways, rownames(comp))
comp <- cbind(path_names, comp)


use_data(brca, pkg = ".", internal = F)

# use_data(brca_data, pkg = ".", internal = F)
# use_data(brca_design, pkg = ".", internal = F)
# use_data(exp_data, pkg = ".", internal = F)
# use_data(results, pkg = ".", internal = F)
# use_data(path_vals, pkg = ".", internal = F)
# use_data(go_vals, pkg = ".", internal = F)
# use_data(uni_vals, pkg = ".", internal = F)
# use_data(comp, pkg = ".", internal = F)



