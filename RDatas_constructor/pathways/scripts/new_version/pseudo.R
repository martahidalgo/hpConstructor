## marta.hidalgo@outlook.es
## Functions to create metaginfo objects grouped by functions or genes


create.pseudo.metaginfo <- function(pathways, group.by, verbose = TRUE){
    
    subgraphs <- unlist(lapply(pathways$pathigraphs, function(pg){
        pg$effector.subgraphs
    }), recursive = FALSE)
    names(subgraphs) <- sapply(names(subgraphs), function(n){
        unlist(strsplit(n, split = "\\."))[2]})
    
    if(group.by == "uniprot" | group.by == "GO"){
        annofuns <- load.annofuns(db = group.by, pathways$species)
        annofuns <- annofuns[!is.na(annofuns$funs),]
        annots <- annofuns[,c(2,3)]
    }else if(group.by == "genes"){
        gens <- lapply(names(subgraphs), function(name){
            sg <- subgraphs[[name]]
            l <- unique(unlist(V(sg)$genesList))
            l <- l[!is.na(l)]
            l <- l[!l == "/"]
            cbind(name, l)
        })
        gens <- gens[sapply(gens, ncol) > 1]
        annots <- do.call("rbind", gens)
    }else{
        stop("Parameter `group.by` not recognized")
    }
    
    pseudo <- get.pseudo.pathigraphs(subgraphs, annots, verbose = verbose)
    
    pseudo.meta <- NULL
    pseudo.meta$pathigraphs <- pseudo
    pseudo.meta$species <- pathways$species
    pseudo.meta$all.labelids <- pathways$all.labelids
    
    return(pseudo.meta)
}



get.pseudo.pathigraphs <- function(subgraphs, annots, verbose = TRUE){
    
    categories <- unique(annots[,2])
    pseudo_pathigraphs <- list()
    for(i in 1:length(categories)){
        if(verbose == TRUE)
            cat("Processing ", categories[i], " (", i, " of ",
                length(categories), ")\n", sep = "")
        selnames <- annots[which(annots[,2] == categories[i]),1]
        selsub <- subgraphs[selnames]
        cat("    found ", length(selsub), " subgraphs...\n", sep = "")
        path.id <- paste0("term_", i)
        cpfs <- create.pathigraph.from.subgraphs(selsub, categories[i], path.id)
        pseudo_pathigraphs[[path.id]] <- cpfs
    }
    
    return(pseudo_pathigraphs)
    
}


create.pathigraph.from.subgraphs <- function(selsub, path.name, path.id){
    
    # Color in red node including the gene, when necessary
    if(any(sapply(selsub, function(g) path.name %in% unlist(V(g)$genesList))))
        for(i in 1:length(selsub)){
            V(selsub[[i]])$label.color[which(sapply(V(selsub[[i]])$genesList, 
                                                    function(gl){
                                                        path.name %in% gl
                                                    }))] <- "#B36100"
            V(selsub[[i]])$stroke.color[which(sapply(V(selsub[[i]])$genesList, 
                                                     function(gl){
                                                         path.name %in% gl
                                                     }))] <- "#B36100"
        }
    
    # overlapping between paths
    imat <- data.matrix(mat.or.vec(nr = length(selsub), nc = length(selsub)))
    rownames(imat) <- names(selsub)
    colnames(imat) <- names(selsub)
    for(i in 1:length(selsub)){
        for(j in 1:length(selsub)){
            imat[i,j] <- length(intersect(V(selsub[[i]])$name,
                                          V(selsub[[j]])$name))
        }
    }
    ga <- graph.adjacency(imat > 0)
    cga <- clusters(ga)
    cga$groups <- unique(cga$membership)
    
    # recalculate Y coordinates
    path_coords <- do.call("rbind", lapply(selsub, function(x){
        range(V(x)$nodeY)
    }))
    init_path_coords <- cbind(path_coords, 0, 0, 0, 0, 0)
    group_coords <- do.call("rbind", by(path_coords, cga$membership, range))
    
    ymargin <- 50
    lasty <- 0
    for(i in 1:cga$no){
        m <- cga$groups[i]
        indexes <- which(cga$membership == m)
        for(j in indexes){
            V(selsub[[j]])$nodeY <- V(selsub[[j]])$nodeY -
                group_coords[i,1] + lasty
            init_path_coords[j,3] <- group_coords[i,1]
            init_path_coords[j,4] <- lasty
            init_path_coords[j,5] <- m
            init_path_coords[j,6] <- path_coords[j,1] -
                group_coords[i,1] + lasty
            init_path_coords[j,7] <- path_coords[j,2] -
                group_coords[i,1] + lasty
        }
        lasty <- lasty + (group_coords[i,2] - group_coords[i,1]) + ymargin
    }
    
    # create general graphs
    
    node_atts <- unique(do.call("rbind", lapply(selsub, function(x) {
        data.frame(name = V(x)$name,
                   label = V(x)$label,
                   shape = V(x)$shape,
                   x = V(x)$nodeX,
                   y = V(x)$nodeY,
                   width = V(x)$width,
                   height = V(x)$height,
                   label.color = V(x)$label.color,
                   label.cex = V(x)$label.cex,
                   genesList = sapply(V(x)$genesList, paste, collapse = ","),
                   stringsAsFactors = FALSE)
    })))
    rownames(node_atts) <- node_atts$name
    
    edges <- unique(do.call("rbind", lapply(selsub, function(x) {
        el <- get.edgelist(x)
        data.frame(source = el[,1],
                   relation = E(x)$relation,
                   target = el[,2])
    })))
    
    supergraph <- graph.data.frame(edges[,c(1,3)], directed = TRUE)
    E(supergraph)$relation <- edges[,2]
    for(k in 1:ncol(node_atts)){
        field <- colnames(node_atts)[k]
        supergraph <- set.vertex.attribute(supergraph,
                                           name = field,
                                           value = node_atts[V(supergraph)$name,
                                                             field])
    }
    V(supergraph)$nodeX <- V(supergraph)$x
    V(supergraph)$nodeY <- V(supergraph)$y
    V(supergraph)$genesList <- sapply(V(supergraph)$genesList, strsplit,
                                      split = ",")
    
    pathigraph <- list()
    pathigraph$graph <- supergraph
    pathigraph$path.name <- path.name
    pathigraph$path.id <- path.id
    pathigraph$effector.subgraphs <- selsub
    pathigraph$path_coords <- path_coords
    pathigraph$init_path_coords <- init_path_coords
    pathigraph$group_coords <- group_coords
    pathigraph$ga <- ga
    pathigraph$cga <- cga
    
    return(pathigraph)
    
}
