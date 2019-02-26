## marta.hidalgo@outlook.es
## Create metaginfo objects 

rm(list = ls())

library(KEGGgraph)
library(igraph)
library(graph)
library(hipathia)
# library(hpAnnot)

path <- "~/appl/hpConstructor/RDatas_constructor/pathways/"
setwd(path)
hipath <- paste0(path, "scripts/")
source(paste0(hipath, "/graphs.R"))
source(paste0(hipath, "/KEGG_net.R"))
source(paste0(hipath, "/layout.R"))
source("~/appl/hipathia/R/utils.R")
# source(paste0(hipath, "/R/load.R"))

ammend.file <- paste0(hipath, "/../sif_amendments.txt")
comp.file <- paste0(hipath, "/../compounds_list.txt")

# FUNCTION
load_compounds <- function(comp.file){
    compounds <- utils::read.table(comp.file, header = FALSE, sep = "\t", 
                                   stringsAsFactors = FALSE, 
                                   colClasses = "character", quote = "")
    colnames(compounds) <- c("ID", "names")
    compounds$simple.ID <- gsub("cpd:", "", compounds$ID)
    compounds$simple.name <- sapply(strsplit(compounds$names, split = ";"), "[[", 1) 
    return(compounds[,c(3,4)])
}

compounds <- load_compounds(comp.file)
# save(compounds, file=paste0(hipath, "/../compounds_list.RData"))

# Parameters
species <- c("hsa", "rno", "mmu")
date <- gsub("-", "_", Sys.Date())

for(spe in species){

    # set folders
    kgml.folder <- paste0(hipath, "/../", spe, "/kgml/")
    sif.folder <- paste0(hipath, "/../", spe, "/sif/")
    tmp.folder <- paste0(hipath, "/../", spe, "/temp/")
    pathway.names <- unique(gsub(".xml", "", list.files(kgml.folder, 
                                                        pattern="xml")))
    
    # Create folders
    if(!dir.exists(sif.folder))
        dir.create(sif.folder)
    if(!dir.exists(tmp.folder))
        dir.create(tmp.folder)
    
    # Load annotations
    dbannot <- hipathia:::load_annots("uniprot", spe)
    entrez2hgnc <- hipathia:::load_entrez_hgnc(spe)

    # Process KGML files
    #-------------------------------------------------
    # create name_pathways.txt file
    create.pathway.names.file(pathway.names, spe, kgml.folder, sif.folder)

    # Transform KGML to SIF files
    transform.XML.to.SIF(pathway.names, kgml.folder, sif.folder, compounds)

    # Load pathways from created SIF files
    pgs <- load.graphs(sif.folder, spe)
    save(pgs, file=paste0(tmp.folder, "/pgs_", date, ".RData"))

    # Ammend pathways
    apgs <- amend.kegg.pathways(ammend.file, pgs, spe)
    save(apgs, file=paste0(tmp.folder, "/apgs_", date, ".RData"))

    # Add final functions to the pathways
    fpgs <- add.functions.to.pathigraphs(apgs, entrez2hgnc, dbannot, 
                                         maxiter = 1000)
    save(fpgs, file=paste0(tmp.folder, "/fpgs_", date, ".RData"))

    # Compute Path Normalization Values
    metaginfo <- create.metaginfo.object(fpgs, spe)
    save(metaginfo, file=paste0(tmp.folder, "/meta_graph_info_", spe,
                                "_", date, ".RData"))

}


