
library(igraph)
library(devtools)
library(hipathia)
library(hpAnnot)

species <- c("hsa", "mmu", "rno")
rda_path <- "RDatas_constructor/RDatas/"
version <- "v2"

for(spe in species){

    print(spe)
    mgi <- load(paste0(rda_path, "/meta_graph_info_", spe, "_", version, ".rda"))
    metaginfo <- get(mgi)
    
    # Uniprot
    unidb <- load(paste0(rda_path,  "annot_uniprot_", spe, "_", version, ".rda"))
    uni_bp_annot <- get(unidb)
    annofuns <- hipathia:::annotate_paths(metaginfo, uni_bp_annot)
    save(annofuns, file = paste0("RDatas_constructor/annofuns/annofuns_uniprot_",
                                 spe, ".RData"))

    # GO
    godb <- load(paste0(rda_path, "annot_GO_", spe, "_", version, ".rda"))
    go_bp_annot <- get(godb)
    annofuns <- hipathia:::annotate_paths(metaginfo, go_bp_annot)
    save(annofuns, file = paste0("RDatas_constructor/annofuns/annofuns_GO_",
                                 spe, ".RData"))


}

