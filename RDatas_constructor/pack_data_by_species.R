

# Script for storing data generated in the RDatas_constructor folder into the 
# data folder as Package Data


species <- c("hsa", "mmu", "rno")
dbs <- c("GO", "uniprot")
setwd("RDatas_constructor/")
edpath <- "RDatas/"
version <- "v2"

# XREF
#------------------------
for(spe in species){
    xref_folder <- paste0("xref/", spe, "/")
    xref_file <- paste0(xref_folder, "xref_", spe, ".RData")
    xr <- load(xref_file)
    xref <- get(xr)
    xref_file <- paste0(edpath, "xref_", spe, "_", version, ".rda")
    save(xref, file = xref_file, compress = "bzip2")
}


# ANNOTATIONS
#------------------------
load("annotations/annotations/go_bp_frame.RData")
goframe_file <- paste0(edpath, "go_bp_frame_", version, ".rda")
save(go_bp_frame, file = goframe_file, compress = "bzip2")

load("annotations/annotations/go_bp_net.RData")
gonet_file <- paste0(edpath, "go_bp_net_", version, ".rda")
save(go_bp_net, file = gonet_file, compress = "bzip2")

# GO_bp_annots
for(spe in species){
    gba_folder <- paste0("annotations/annotations/", spe, "/")
    gba_file <- paste0(gba_folder, "/go_bp_", spe, ".annot")
    go_bp_annot <- read.delim(gba_file,
                              sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE,
                              comment.char = "!")
    go_bp_annot <- go_bp_annot[,1:3]
    colnames(go_bp_annot) <- c("gene", "function", "term")
    bpfile <- paste0(edpath, "annot_GO_", spe, "_", version, ".rda")
    save(go_bp_annot, file = bpfile, compress = "bzip2")
}


# Uni_bp_annots
for(spe in species){
    gba_folder <- paste0("annotations/annotations/", spe, "/")
    uba_file <- paste0(gba_folder, "/uniprot_keywords_",
                       spe, "__biological_process.annot")
    uni_bp_annot <- read.delim(uba_file,
                              sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE,
                              comment.char = "!")
    uni_bp_annot <- uni_bp_annot[,1:2]
    colnames(uni_bp_annot) <- c("gene","function")
    ubpfile <- paste0(edpath, "annot_uniprot_", spe, "_", version, ".rda")
    save(uni_bp_annot, file = ubpfile, compress = "bzip2")
}

# Entrez_HGNC
for(spe in species){
    eh_folder <- paste0("annotations/annotations/", spe, "/")
    eh_file <- paste0(eh_folder, "/entrez_hgnc_", spe, ".annot")
    entrez_hgnc <- read.delim(eh_file,
                              sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE,
                              comment.char = "!")
    ehfile <- paste0(edpath, "entrez_hgnc_", spe, "_", version, ".rda")
    save(entrez_hgnc, file = ehfile, compress = "bzip2")
}



# PATHWAYS
#------------------------
# All together in different files
for(spe in species){
    path <- paste0("pathways/", spe, "/temp/meta_graph_info_",
                   spe, ".RData")
    mgi <- load(path)
    meta_graph_info <- get(mgi)
    meta_graph_info$species <- spe
    mgifile <- paste0(edpath, "meta_graph_info_", spe, "_", version, ".rda")
    save(meta_graph_info, file = mgifile, compress = "bzip2")
}



# ANNOFUNS
#------------------------
for(db in dbs){
    for(spe in species){
        load(paste0("annofuns/annofuns_", db, "_", spe, ".RData"))
        af <- paste0(edpath, "annofuns_", db, "_", spe, "_", version, ".rda")
        save(annofuns, file = af, compress = "bzip2")
    }
}



# PSEUDO PATHWAYS
#------------------------
feats <- c("genes", "uniprot", "GO")
for(spe in species){
    for(feat in feats){
        file <- paste0("pathways/pseudo/pmgi_", spe, "_", feat,".RData")
        pseudo <- load(file)
        pmgi <- get(pseudo)
        pf <- paste0(edpath, "pmgi_", spe, "_", feat, "_", version, ".rda")
        save(pmgi, file = pf, compress = "bzip2")
        rm(list = pseudo)
    }
}

