


pathway_index <- lapply(c("hsa", "mmu", "rno"), function(spe){
    mgi <- load_pathways(spe)
    pi <- data.frame(unique(mgi$all.labelids[,c(3,4)]))
    rownames(pi) <- NULL
    pi
})
names(pathway_index) <- c("hsa", "mmu", "rno")
pathway_index$hsa

save(pathway_index, file= "~/appl/hpConstructor/code_develop/pathway_categories/pathway_index.RData")
