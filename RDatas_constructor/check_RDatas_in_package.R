
## LIBRARIES
library(devtools)
load_all("~/appl/hipathia/")
load_all("~/appl/hpAnnot/")
source("~/appl/hipathia/R/load_hpAnnot.R")


## PARAMETERS
species <- c("hsa", "mmu", "rno")
groupings <- c("uniprot", "GO", "genes")

## CHECKING
for(i in 1:length(species)){
    spe <- species[i]
    
    # Load pathways
    pathways <- load_pathways(spe)

    if(spe == "hsa"){
        data("brca")
        g1 <- "Tumor"
        g2 <- "Normal"
    }else if(spe == "rno"){
        data <- read.delim("~/UBB_projects/2018/2018_12_Vicky_Moreno/saved_data/norm_data.txt", 
                           header = T, sep = "\t")
        data <- as.matrix(data)
        ed <- load("~/UBB_projects/2018/2018_12_Vicky_Moreno/saved_data/covs.RData")
        ed <- get(ed)
        ed <- ed[,c(2,3)]
        colnames(ed) <- c("sample", "group")
        brca <- SummarizedExperiment(assays=SimpleList(raw=data), 
                                     colData=ed)
        g1 <- "t1"
        g2 <- "t0"
    }else if(spe == "mmu"){
        data <- read.delim("/data/projectsData/180_checkSpecies/mouse_norm.txt", 
                           header = T, sep = "\t")
        data <- as.matrix(data)
        ed <- read.delim("/data/projectsData/180_checkSpecies/mouse_ED.txt", 
                         header = F, sep = "\t")
        colnames(ed) <- c("sample", "group")
        brca <- SummarizedExperiment(assays=SimpleList(raw=data), 
                                     colData=ed)
        g1 <- "ko"
        g2 <- "control"
    }else if(spe == "dme"){
        load("/data/datos/hipathia_package_test/genital_disc.Rda")
        brca <- genital_disc
        colData(brca) <- DataFrame(colData(brca), group = colData(brca)[["Stage"]])
        g1 <- "P6"
        g2 <- "P20"
    }
    
    brca_trans <- translate_data(brca, spe)
    boxplot(assay(brca_trans))
    brca_norm <- normalize_data(brca_trans)

    # Hipathia 
    brca_results <- hipathia(brca_norm, pathways, verbose = F)
    brca_paths <- get_paths_data(brca_results)

    # Comparison
    comp <- do_limma(brca_paths, "group", g1, g2, order = T)

    colors_de <- node_color_per_de(brca_results, pathways, "group", g1, g2)
    # Save and serve all results to browser
    temp.folder <- create_report(comp, pathways, node_colors = colors_de)
    visualize_report(temp.folder, port = 4000 + 4*i)

    # Subpathways grouped by functions (This may take a bit longer)
    for( j in 1:length(groupings)){
        group <- groupings[j]
        colors_group <- node_color_per_de(brca_results, pathways, "group", g1,
                                        g2, group_by = group)
        temp.folder <- create_report(comp, pathways, node_colors = colors_group, group_by = group)
        visualize_report(temp.folder, port = 4000 + 4*i + j)
    }
    
    
    # FUNCTIONAL ANALYSIS
    brca_uniprot <- quantify_terms(brca_results, pathways, "uniprot")
    heatmap_plot(brca_uniprot, "group", variable_clust = TRUE)
    comp_uniprot <- do_wilcoxon(brca_uniprot, "group", g1, g2, order = T)
    head(comp_uniprot)
    
    brca_go <- quantify_terms(brca_results, pathways, "GO")
    heatmap_plot(brca_go, "group", variable_clust = TRUE)
    comp_go <- do_wilcoxon(brca_go, "group", g1, g2, order = TRUE)
    head(comp_go)
    
    
}
