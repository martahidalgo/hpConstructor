
library(igraph)
library(devtools)
library(hipathia)
library(hpAnnot)

species <- c("hsa", "mmu", "rno")

for(spe in species){

    print(spe)
    mgi <- load(paste0("private/pathways/", spe, "/temp/meta_graph_info_", 
                      spe, ".RData"))
    metaginfo <- get(mgi)
    
    # Uniprot
    annofuns <- hipathia:::annotate.paths(metaginfo, "uniprot")
    save(annofuns, file = paste0("private/annofuns/annofuns_uniprot_",
                                 spe, ".RData"))

    # GO
    annofuns <- hipathia:::annotate.paths(metaginfo, "GO")
    save(annofuns, file = paste0("private/annofuns/annofuns_GO_",
                                 spe, ".RData"))


}

