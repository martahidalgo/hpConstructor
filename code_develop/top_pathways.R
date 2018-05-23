library(hipathia)
library(mdgsa)

data(exp_data)
data("brca_design")
mgi <- load_pathways("hsa")
results <- hipathia(exp_data, mgi)
new_comp <- do_wilcoxon(results$all$path.vals, brca_design$group, "Tumor", "Normal")

rindex <- pval2index(pval = new_comp$p.value, sign = new_comp$statistic, names = rownames(new_comp))
rindex <- indexTransform(rindex)
plot(new_comp$statistic, rindex)
plot(new_comp$p.value, rindex)

pathways <- sapply(strsplit(names(rindex), split = "-"), "[[", 2)
m <- cbind(paths = names(rindex), pathways)
annot <- annotMat2list(m)
# annot <- annotFilter(annot, rindex)

res <- uvGsa(rindex, annot)
res

# Sólo salen 3 pathways enriquecidos con este método, los que son indudablemente
# DOWN. Hipotetizo que debe ser porque los que tienen subpaths UP y DOWN no 
# salen significativos porque sus subpaths están repartidos a lo largo del valor 
# del estadístico


# Probamos si haciendo que todos los estadísticos vayan en la misma dirección 
# (sean positivos) corregimos este problema.

rindex <- pval2index(pval = new_comp$p.value, sign = abs(new_comp$statistic), names = rownames(new_comp))
rindex <- indexTransform(rindex)
plot(abs(new_comp$statistic), rindex)
plot(new_comp$p.value, rindex)

pathways <- sapply(strsplit(names(rindex), split = "-"), "[[", 2)
m <- cbind(paths = names(rindex), pathways)
annot <- annotMat2list(m)
# annot <- annotFilter(annot, rindex)

res <- uvGsa(rindex, annot)
res


# No parece que consigamos corregir este efecto, en este caso sólo salen 6 
# subpaths significativos.


# Probamos a hacer un Fisher Test comparando la proporción de pathways 
# significativos en cada pathway, con la proporción del resto de pathways.
new_comp$pathways <- sapply(strsplit(rownames(new_comp), split = "-"), "[[", 2)

path <- "hsa03320"
tests <- do.call(rbind, lapply(paths[,1], function(path) {
    t1 <- table(new_comp[new_comp$pathways == path,"FDRp.value"] < 0.05)
    t2 <- table(new_comp[new_comp$pathways != path,"FDRp.value"] < 0.05)
    
    ft <- fisher.test(rbind(t1, t2))
    data.frame(path = path, pval = ft$p.value, stat = ft$estimate)
}))

tests$path.name <- paths[,2]
tests <- tests[,c(1,4,2,3)]
tests <- tests[order(tests$pval),]
tests[tests$pval < 0.05 & tests$stat < 1,]



# Ahora comparamos si la proporción de pathways significativos de un pathway
# es significativa (comparando contra no tener ningún pathway significativo)

path <- "hsa03320"
tests2 <- do.call(rbind, lapply(paths[,1], function(path) {
    t1 <- table(new_comp[new_comp$pathways == path,"FDRp.value"] < 0.05)
    t2 <- c(sum(t1), 0)
    
    ft <- fisher.test(rbind(t1, t2))
    data.frame(path.id = path, pval = ft$p.value, stringsAsFactors = F)
}))

all(paths[,1] == tests2$paths)
tests2$path.name <- paths[,2]
rownames(tests2) <- tests2$path.id
tests2 <- tests2[,c(1,3,2)]
tests2 <- tests2[order(tests2$pval),]
table(tests2$pval < 0.05)
tests2[tests2$pval < 0.05,]

# Esto parece que no va mal, asigna un p-valor a cada pathway según la proporción
# de subpaths significativos, y los ordena de mayor proporción a menor. Al menos 
# como medida de qué pathways mirar primero puede ser interesante.


