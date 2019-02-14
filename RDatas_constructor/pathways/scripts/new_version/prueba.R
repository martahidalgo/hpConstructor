removeInnecessaryBindings <- function(graph, verbose = TRUE){
    es <- graph@edgeData@defaults$KEGGEdge$edges
    types <- sapply(es, function(x) if(!is.null(x@subtype$subtype)){x@subtype$subtype@name}else{x@type})
    nodes1 <- as.character(sapply(es, function(x) x@entry1ID))
    nodes2 <- as.character(sapply(es, function(x) x@entry2ID))
    es.mat <- data.frame(nodes1, nodes2, types, stringsAsFactors = F)
    es.mat$include <- apply(es.mat[,c(1,2)], 1, function(v){
        v1 <- unlist(strsplit(v[1], split = " "))
        v2 <- unlist(strsplit(v[2], split = " "))
        all(v1 %in% v2) | all(v2 %in% v1)
    })
    if(any(es.mat$include)){
        head(es.mat)
        toremove <- which(es.mat$include == TRUE & es.mat$type == "binding/association")
        for(i in toremove){
            if(verbose)
                message("Removing edge: ", es.mat$nodes1[i], " -> ", es.mat$nodes2[i])
            graph <- removeKEGGgraphEdge(es.mat$nodes1[i], es.mat$nodes2[i], graph)
        }
    }
    return(graph)
}



removeKEGGgraphEdge <- function(n1, n2, graph){
    g <- removeEdge(n1, n2, graph)
    eds <- g@edgeData@defaults$KEGGEdge$edges
    index <- which(names(eds) == paste0(n1, "~", n2))
    g@edgeData@defaults$KEGGEdge$edges <- eds[-index]
    return(g)
}
