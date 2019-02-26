## marta.hidalgo@outlook.es
## Create metaginfo objects grouped by functions or genes


rm(list = ls())

library(hipathia)

hipath <- getwd()
hipath <- paste0(hipath, "/RDatas_constructor/pathways/")
source(paste0(hipath, "scripts/pseudo.R"))
date <- gsub("-", "_", Sys.Date())

# Parameters
species <- c("hsa", "rno", "mmu")
groupings <- c("uniprot", "GO", "genes")

pmgis <- NULL
for(spe in species){
    print(spe)
    spe.pmgis <- NULL
    mgi <- load_pathways(spe)
    for(group in groupings){
        print(group)
        pseudo.mgi <- create.pseudo.metaginfo(mgi, group)
        spe.pmgis[[group]] <- pseudo.mgi
    }
    pmgis[[spe]] <- spe.pmgis
}

sapply(species, function(spe) sapply(groupings, function(group){
    file <- paste0(hipath, spe, "/temp/pmgi_", spe, "_", group, "_", date, ".RData")
    tosave <- pmgis[[spe]][[group]]
    save(tosave, file = file)
}))

# 
# pmgi_hsa_uniprot <- pmgis[["hsa"]][["uniprot"]]
# save(pmgi_hsa_uniprot, file = "private/pathways/pseudo/pmgi_hsa_uniprot.RData")
# 
# pmgi_hsa_GO <- pmgis[["hsa"]][["GO"]]
# save(pmgi_hsa_GO, file = "private/pathways/pseudo/pmgi_hsa_GO.RData")
# 
# pmgi_hsa_genes <- pmgis[["hsa"]][["genes"]]
# save(pmgi_hsa_genes, file = "private/pathways/pseudo/pmgi_hsa_genes.RData")
# 
# pmgi_mmu_uniprot <- pmgis[["mmu"]][["uniprot"]]
# save(pmgi_mmu_uniprot, file = "private/pathways/pseudo/pmgi_mmu_uniprot.RData")
# 
# pmgi_mmu_GO <- pmgis[["mmu"]][["GO"]]
# save(pmgi_mmu_GO, file = "private/pathways/pseudo/pmgi_mmu_GO.RData")
# 
# pmgi_mmu_genes <- pmgis[["mmu"]][["genes"]]
# save(pmgi_mmu_genes, file = "private/pathways/pseudo/pmgi_mmu_genes.RData")
# 
# pmgi_rno_uniprot <- pmgis[["rno"]][["uniprot"]]
# save(pmgi_rno_uniprot, file = "private/pathways/pseudo/pmgi_rno_uniprot.RData")
# 
# pmgi_rno_GO <- pmgis[["rno"]][["GO"]]
# save(pmgi_rno_GO, file = "private/pathways/pseudo/pmgi_rno_GO.RData")
# 
# pmgi_rno_genes <- pmgis[["rno"]][["genes"]]
# save(pmgi_rno_genes, file = "private/pathways/pseudo/pmgi_rno_genes.RData")
# 
