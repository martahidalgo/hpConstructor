

# Script for storing data generated in the private folder into the data folder
# as Package Data


species <- c("hsa", "mmu", "rno")
dbs <- c("GO", "uniprot")
edpath <- "inst/extdata/"

# XREF
#------------------------
for(spe in species){
    xref_folder <- paste0("private/xref/", spe, "/")
    xref_file <- paste0(xref_folder, "xref_", spe, ".RData")
    xr <- load(xref_file)
    xref <- get(xr)
    save(xref, file=paste0(edpath, "xref_", spe, ".rda"), compress = "bzip2")
}


# ANNOTATIONS
#------------------------
load("private/annotations/annotations/go_bp_frame.RData")
save(go_bp_frame, file = paste0(edpath, "go_bp_frame.rda"), compress = "bzip2")
load("private/annotations/annotations/go_bp_net.RData")
save(go_bp_net, file = paste0(edpath, "go_bp_net.rda"), compress = "bzip2")

# GO_bp_annots
for(spe in species){
    gba_folder <- paste0("private/annotations/annotations/", spe, "/")
    gba_file <- paste0(gba_folder, "/go_bp_", spe, ".annot")
    go_bp_annot <- read.delim(gba_file,
                              sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE,
                              comment.char = "!")
    go_bp_annot <- go_bp_annot[,1:2]
    colnames(go_bp_annot) <- c("gene","function")
    bpfile <- paste0(edpath, "annot_GO_", spe, ".rda")
    save(go_bp_annot, file = bpfile, compress = "bzip2")
}


# Uni_bp_annots
for(spe in species){
    gba_folder <- paste0("private/annotations/annotations/", spe, "/")
    uba_file <- paste0(gba_folder, "/uniprot_keywords_",
                       spe, "__biological_process.annot")
    uni_bp_annot <- read.delim(uba_file,
                              sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE,
                              comment.char = "!")
    uni_bp_annot <- uni_bp_annot[,1:2]
    colnames(uni_bp_annot) <- c("gene","function")
    ubpfile <- paste0(edpath, "annot_uniprot_", spe, ".rda")
    save(uni_bp_annot, file = ubpfile, compress = "bzip2")
}

# Entrez_HGNC
for(spe in species){
    eh_folder <- paste0("private/annotations/annotations/", spe, "/")
    eh_file <- paste0(eh_folder, "/entrez_hgnc_", spe, ".annot")
    entrez_hgnc <- read.delim(eh_file,
                              sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE,
                              comment.char = "!")
    ehfile <- paste0(edpath, "entrez_hgnc_", spe, ".rda")
    save(entrez_hgnc, file = ehfile, compress = "bzip2")
}



# PATHWAYS
#------------------------
# All together in different files
for(spe in species){
    path <- paste0("private/pathways/", spe, "/temp/meta_graph_info_",
                   spe, ".RData")
    mgi <- load(path)
    meta_graph_info <- get(mgi)
    meta_graph_info$species <- spe
    mgifile <- paste0(edpath, "meta_graph_info_", spe, ".rda")
    save(meta_graph_info, file = mgifile, compress = "bzip2")
}


# # Separatedly
# # HSA
# mgi <- load("private/pathways/hsa/temp/meta_graph_info_hsa.RData")
# meta_graph_info_hsa <- get(mgi)
# meta_graph_info_hsa$species <- "hsa"
# save(meta_graph_info_hsa, pkg = ".")
# rm(metaginfo, meta_graph_info_hsa, mgi)
# 
# # RNO
# mgi <- load("private/pathways/rno/temp/meta_graph_info_rno.RData")
# meta_graph_info_rno <- get(mgi)
# meta_graph_info_rno$species <- "rno"
# save(meta_graph_info_rno, pkg = ".")
# rm(metaginfo, meta_graph_info_rno, mgi)
# 
# # MMU
# mgi <- load("private/pathways/mmu/temp/meta_graph_info_mmu.RData")
# meta_graph_info_mmu <- get(mgi)
# meta_graph_info_mmu$species <- "mmu"
# save(meta_graph_info_mmu, pkg = ".")
# rm(metaginfo, meta_graph_info_mmu, mgi)




# ANNOFUNS
#------------------------
for(db in dbs){
    for(spe in species){
        load(paste0("private/annofuns/annofuns_", db, "_", spe, ".RData"))
        save(annofuns, file = paste0(edpath, "annofuns_", db, "_", spe, ".rda"), 
             compress = "bzip2")
    }
}



# PSEUDO PATHWAYS
#------------------------
feats <- c("genes", "uniprot", "GO")
for(spe in species){
    for(feat in feats){
        file <- paste0("private/pathways/pseudo/pmgi_", spe, "_", feat,".RData")
        pseudo <- load(file)
        pmgi <- get(pseudo)
        save(pmgi, file = paste0(edpath, "pmgi_", spe, "_", feat, ".rda"), 
             compress = "bzip2")
        rm(list = pseudo)
    }
}

# # HSA
# pmgi <- load("private/pathways/pseudo/pmgi_hsa_genes.RData")
# pmgi_hsa_genes <- get(pmgi)
# save(pmgi_hsa_genes, pkg = ".")
# rm(pmgi_hsa_genes, pmgi)
# 
# pmgi <- load("private/pathways/pseudo/pmgi_hsa_uniprot.RData")
# pmgi_hsa_uniprot <- get(pmgi)
# save(pmgi_hsa_uniprot, pkg = ".")
# rm(pmgi_hsa_uniprot, pmgi)
# 
# pmgi <- load("private/pathways/pseudo/pmgi_hsa_GO.RData")
# pmgi_hsa_GO <- get(pmgi)
# save(pmgi_hsa_GO, pkg = ".")
# rm(pmgi_hsa_GO, pmgi)
# 
# # RNO
# pmgi <- load("private/pathways/pseudo/pmgi_rno_genes.RData")
# pmgi_rno_genes <- get(pmgi)
# save(pmgi_rno_genes, pkg = ".")
# rm(pmgi_rno_genes, pmgi)
# 
# pmgi <- load("private/pathways/pseudo/pmgi_rno_uniprot.RData")
# pmgi_rno_uniprot <- get(pmgi)
# save(pmgi_rno_uniprot, pkg = ".")
# rm(pmgi_rno_uniprot, pmgi)
# 
# pmgi <- load("private/pathways/pseudo/pmgi_rno_GO.RData")
# pmgi_rno_GO <- get(pmgi)
# save(pmgi_rno_GO, pkg = ".")
# rm(pmgi_rno_GO, pmgi)
# 
# # MMU
# pmgi <- load("private/pathways/pseudo/pmgi_mmu_genes.RData")
# pmgi_mmu_genes <- get(pmgi)
# save(pmgi_mmu_genes, pkg = ".")
# rm(pmgi_mmu_genes, pmgi)
# 
# pmgi <- load("private/pathways/pseudo/pmgi_mmu_uniprot.RData")
# pmgi_mmu_uniprot <- get(pmgi)
# save(pmgi_mmu_uniprot, pkg = ".")
# rm(pmgi_mmu_uniprot, pmgi)
# 
# pmgi <- load("private/pathways/pseudo/pmgi_mmu_GO.RData")
# pmgi_mmu_GO <- get(pmgi)
# save(pmgi_mmu_GO, pkg = ".")
# rm(pmgi_mmu_GO, pmgi)
# 
# 
