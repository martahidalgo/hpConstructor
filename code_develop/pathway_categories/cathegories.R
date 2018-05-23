
library(XML)
setwd("appl/hipathia/private/")

html <- readLines("pathway_categories/KEGG PATHWAY Database.html")

# Take lines with middle cathegory
bs <- grep("<b>", html)
# Take lines with high cathegory
h4s <- grep("<h4>", html)
# Take lines with pathways
dts <- grep("<dt>", html)

# # Remove "Metabolism" and "Drug Development" labels from h4s
# h4s <- h4s[2:6]

# Remove middle cathegories and pathways before first h4 class
bs <- bs[bs > h4s[1]]
dts <- dts[dts > h4s[1]]

# Make table
t <- do.call("rbind", lapply(dts, function(dt){
    b <- max(bs[bs < dt])
    middle.cat <- unlist(strsplit(html[b], split = " "))
    middle.cat <- middle.cat[-1]
    middle.cat <- paste(middle.cat, collapse = " ")
    middle.cat <- unlist(strsplit(middle.cat, split = "<"))[1]
    
    h4 <- max(h4s[h4s < dt])
    high.cat.split <- unlist(strsplit(html[h4], split = "\\. "))
    high.cat <- unlist(strsplit(high.cat.split[2], split = "<"))[1]
    
    pathway <- html[dt]
    p <- unlist(strsplit(pathway, split = "<"))
    pathway.name <- unlist(strsplit(p[5], split = ">"))[2]
    pathway.id <- unlist(strsplit(p[2], split = ">"))[2]
    
    data.frame(Pathway.id = pathway.id, Pathway.name = pathway.name, 
               Subcathegory = middle.cat, Cathegory = high.cat, 
               stringsAsFactors = FALSE)
}))

cats <- t

save(cats, file = "pathway_categories/pathway_cathegories.RData")
