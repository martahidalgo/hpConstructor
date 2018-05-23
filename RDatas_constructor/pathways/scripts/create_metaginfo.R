## marta.hidalgo@outlook.es
## Create metaginfo objects 

rm(list = ls())

library(KEGGgraph)
library(igraph)
library(graph)
library(hipathia)
library(hpAnnot)


hipath <- getwd()
source(paste0(hipath, "/private/pathways/scripts/graphs.R"))
source(paste0(hipath, "/private/pathways/scripts/KEGG_net.R"))
source(paste0(hipath, "/private/pathways/scripts/layout.R"))
source(paste0(hipath, "/../hipathia/R/utils.R"))
# source(paste0(hipath, "/R/load.R"))

# Parameters
species <- c("hsa", "rno", "mmu")

for(spe in species){

    # set folders
    kgml.folder <- paste0(hipath, "/private/pathways/", spe, "/kgml/")
    sif.folder <- paste0(hipath, "/private/pathways/", spe, "/sif/")
    tmp.folder <- paste0(hipath, "/private/pathways/", spe, "/temp/")
    ammend.file <- paste0(hipath, "/private/pathways/sif_amendments.txt")
    pathway.names <- unique(gsub(".xml", "", list.files(kgml.folder, 
                                                        pattern="xml")))
    
    # Create folders
    if(!dir.exists(sif.folder))
        dir.create(sif.folder)
    if(!dir.exists(tmp.folder))
        dir.create(tmp.folder)
    
    # Load annotations
    dbannot <- hipathia:::load.annots("uniprot", spe)
    entrez2hgnc <- hipathia:::load.entrez.hgnc(spe)


    # Process KGML files
    #-------------------------------------------------
    # create name_pathways.txt file
    create.pathway.names.file(pathway.names, spe, kgml.folder, sif.folder)

    # Transform KGML to SIF files
    transform.XML.to.SIF(pathway.names, kgml.folder, sif.folder)

    # Load pathways from created SIF files
    pgs <- load.graphs(sif.folder, spe)
    save(pgs, file=paste0(tmp.folder, "/pgs.RData"))

    # Ammend pathways
    apgs <- amend.kegg.pathways(ammend.file, pgs, spe)
    save(apgs, file=paste0(tmp.folder, "/apgs.RData"))

    # Add final functions to the pathways
    fpgs <- add.functions.to.pathigraphs(apgs, entrez2hgnc, dbannot, 
                                         maxiter = 1000)
    save(fpgs, file=paste0(tmp.folder, "/fpgs.RData"))

    # Compute Path Normalization Values
    metaginfo <- create.metaginfo.object(fpgs, spe)
    save(metaginfo, file=paste0(tmp.folder, "/meta_graph_info_", spe,
                                ".RData"))

}



