
# SimplifyKEGGNet
# Given a KEGG pathway, it transforms it in a graph, and constructs a new graph with the simplified net.
# 17-06-2014 Marta R. Hidalgo, mhidalgo@cipf.es


KEGGpathway2Graph <- function(pathway, genesOnly=TRUE, expandGenes=TRUE) {
  stopifnot(is(pathway, "KEGGPathway"))
  cat(". ")
  pathway <- splitKEGGgroup(pathway)

  if(expandGenes) {
    pathway <- expandKEGGPathway(pathway)
  }

  knodes <- nodes(pathway)
  kedges <- unique(pathway@edges) ## to avoid duplicated edges

  node.entryIDs <- getEntryID(knodes)
  edge.entryIDs <- getEntryID(kedges)

  ## V as nodes, edL as edges
  V <- node.entryIDs
  edL <- vector("list",length=length(V))
  names(edL) <- V

  if(is.null(nrow(edge.entryIDs))) {## no edge found
    for(i in seq(along=edL)) {
      edL[[i]] <- list()
    }
  } else {
    for(i in 1:length(V)) {
      id <- node.entryIDs[i]
      hasRelation <- id == edge.entryIDs[,"Entry1ID"]
      if(!any(hasRelation)) {
        edL[[i]] <- list(edges=NULL)
      } else {
        entry2 <- unname(edge.entryIDs[hasRelation, "Entry2ID"])
        edL[[i]] <- list(edges=entry2)
      }
    }
  }
  gR <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

  ## set node and edge data - as KEGGNode and KEGGEdge
  ## attention: KEGGEdges may be more than graph edges, due to non-genes
  names(kedges) <- sapply(kedges, function(x) paste(getEntryID(x),collapse="~"))

  env.node <- new.env()
  env.edge <- new.env()
  assign("nodes", knodes, envir=env.node)
  assign("edges", kedges, envir=env.edge)

  nodeDataDefaults(gR, "KEGGNode") <- env.node
  edgeDataDefaults(gR, "KEGGEdge") <- env.edge

  if(genesOnly) {
    gR <- subGraphByNodeType(gR,"gene")
  }

  return(gR)
}


###################################################################
#
# MODIFIED BY mhidalgo@cipf.es
#
# This function was initially part of the KEGGGraph package.
# This function has been modified. Instead of splitting the components of the KEGGGroups in different nodes,
# adds a node representing each group.
#
####################################################################
splitKEGGgroup <- function(pathway) {
  pnodes <- pathway@nodes
  pedges <- pathway@edges
  cat(".")
  if(length(pedges)==0) return(pathway)

  types <- sapply(pnodes, getType)
  if(any(types == "group")) {
    isGroup <- names(pnodes)[types == "group"]
    edgeEntry <- sapply(pedges,getEntryID)
    # groupAsID <- edgeEntry[1L,] %in% isGroup | edgeEntry[2L,] %in% isGroup

    ##########################################
    # ADDED BY mhidalgo@cipf.es
    # new.edges <- list()
    for( g in isGroup ){
      comps <- getComponent(pnodes[[g]])
      comps.node <- lapply(comps, search.node, pnodes)
      new.node <- create.collapsed.KEGGnode(comps.node)
      pnodes[[g]] <- new.node
      names(pnodes)[names(pnodes) == g] <- new.node@entryID
      from.g.as.id <- edgeEntry[1,] == g
      to.g.as.id <- edgeEntry[2,] == g
      to.comps.as.id <- edgeEntry[2,] %in% comps
      neweds.from <- sapply(pedges[from.g.as.id], function(x){x@entry1ID <- new.node@entryID; return(x)})
      neweds.to <- sapply(pedges[to.g.as.id], function(x){x@entry2ID <- new.node@entryID; return(x)})
      neweds <- sapply(pedges[to.comps.as.id], function(x){x@entry2ID <- new.node@entryID; return(x)})
      pedges[from.g.as.id] <- neweds.from
      pedges[to.g.as.id] <- neweds.to
      pedges <- c(pedges, neweds)
      # Remove nodes from the group which do not go anywhere
      edgeEntry <- sapply(pedges,getEntryID)
      torem <- comps[!comps%in%edgeEntry[1,]]
      pnodes <- pnodes[!names(pnodes)%in%torem]
      edtorem <- apply(edgeEntry, 2, function(x){x[1]%in%torem || x[2]%in%torem})
      pedges <- pedges[!edtorem]
      edgeEntry <- sapply(pedges,getEntryID)
    }

    pathway@edges <- pedges
    pathway@nodes <- pnodes
    ##########################################

    ##########################################
    # REMOVED BY mhidalgo@cipf.es
    #     newly <- list()
    #     for (e in pedges[groupAsID]) {
    #       entryIDs <- getEntryID(e)
    #       node1comps <- getComponent(pnodes[[ entryIDs[1] ]])
    #       node2comps <- getComponent(pnodes[[ entryIDs[2] ]])
    #       if(length(node1comps) == 1 && is.na(node1comps)) next;
    #       if(length(node2comps) == 1 && is.na(node2comps)) next;
    #       expandmodel <- expand.grid(node1comps, node2comps)
    #       enews <- list()
    #       for (j in 1:nrow(expandmodel)) {
    #         enews[[j]] <- e
    #         entryID(enews[[j]]) <- c(as.character(expandmodel[j,1L]),as.character(expandmodel[j,2L]))
    #       }
    #       newly <- append(newly, enews)
    #     }
    #    newEdges <- append(newEdges, newly)
    #    newEdges <- pedges[!groupAsID]
    #    edges(pathway) <- newEdges
    #############################################

  }
  return(pathway)
}


parse.kegg.to.igraph <- function(name.pathway, folder, comp.list,
                                 genesOnly = FALSE){
  kegg.pathway <- parseKGML(paste0(folder, "/", name.pathway, ".xml"))

  # ori.pathway <- KEGGpathway2Graph(kegg.pathway, expandGenes=FALSE, genesOnly = genesOnly)
  new.pathway <- KEGGpathway2Graph(kegg.pathway, expandGenes=FALSE, genesOnly = genesOnly)

  new.pathway <- delete.unwanted.edges(new.pathway)
  new.pathway <- collapse.nodes( new.pathway )
  new.pathway <- delete.isolated.nodes( new.pathway )
  nodeDataDefaults( new.pathway, "Title" ) <- kegg.pathway@pathwayInfo@title

  graph <- kegg2igraph(name.pathway, new.pathway, comp.list)

  return(graph)
}




# Returns an igraph created from a KEGGGraph
kegg2igraph <- function( name.pathway, kgraph, comp.list ){

  stopifnot(is(kgraph@nodeData@defaults$KEGGNode$nodes[[1]], "KEGGNode"))
  weights <- add.weight.info( kgraph )

  sif <- do.call("rbind", lapply(names(kgraph@edgeData@data), function(x){
    nodes <- unlist(strsplit(x, "\\|"))
    actin <- ifelse(weights[[sub("\\|", "~", x)]] == 1, "activation", "inhibition")
    return(c(nodes[1], actin, nodes[2]))
  }))

  ig <- graph.data.frame(sif[,c(1,3)], directed=T)

  # Add attributes
  E(ig)$relation <- unlist(sapply(sif[,2], function(x){ifelse(x=="activation", 1, -1)}))
  E(ig)$curved <- F
  coor <- get.xy.vals( kgraph )
  V(ig)$nodeX <- coor[match(V(ig)$name,kgraph@nodes),1]
  V(ig)$nodeY <- max(coor[,2]) - coor[match(V(ig)$name,kgraph@nodes),2]
  shapes <- get.shape( kgraph )
  V(ig)$shape <- shapes[match(V(ig)$name,kgraph@nodes)]
  V(ig)$label.cex <- 0.5
  V(ig)$label.color <- "black"
  width <- sapply(kgraph@nodeData@defaults$KEGGNode$nodes, function(x){ x@graphics@width })
  height <- sapply(kgraph@nodeData@defaults$KEGGNode$nodes, function(x){ x@graphics@height })
  V(ig)$width <- width[match(V(ig)$name,kgraph@nodes)]
  V(ig)$height <- height[match(V(ig)$name,kgraph@nodes)]
  type <- sapply(kgraph@nodeData@defaults$KEGGNode$nodes, function(x){ x@type })
  V(ig)$type <- type[match(V(ig)$name,kgraph@nodes)]
  glist <- get.genes.list.type(kgraph)
  V(ig)$genesList <- glist[match(V(ig)$name,kgraph@nodes)]
  labels <- get.name.nodes( kgraph )
  V(ig)$label <- labels[match(V(ig)$name,kgraph@nodes)]
  species <- paste0(unlist(strsplit(name.pathway, split = ""))[1:3], 
                    collapse = "")
  hgnc <- hipathia:::load_entrez_hgnc(species)
  V(ig)$tooltip <- sapply(1:length(V(ig)), function(i){
      if(V(ig)$shape[i] == "rectangle"){
          paste(sapply(V(ig)$genesList[[i]], function(gen){
              if(!is.na(gen) & !gen == "/"){
                  paste0("<a target='_blank' ",
                         "href='http://www.genome.jp/dbget-bin/www_bget?",
                         species, ":", gen, "'>", gen, "</a> (", 
                         get.hgnc.gene(hgnc, gen), ")")
              }else{
                  ""
              }
          }), collapse="<br>")
          # }else if(V(ig)$shape[i] == "rectangle"){
          #     if(grepl("\n", V(ig)$label[i])){
          #         gsub("\n", "<br>", V(ig)$label[i])
          #     }else{""}
      }else if(V(ig)$shape[i] == "circle"){
          link <- kgraph@nodeData@defaults$KEGGNode$nodes[[V(ig)$name[i]]]@link
          paste0("<a target='_blank' href='", link, "'>", V(ig)$label[i], "</a>")
          # if(grepl("BR:", V(ig)$label[i]) == TRUE){
          #     n <- gsub("BR:", "", unlist(strsplit(V(ig)$label[i], split = "\\("))[1])
          #     paste0("<a target='_blank' ",
          #            "href='http://www.genome.jp/kegg-bin/get_htext?", n, "'>", 
          #            n, "</a>")
          # }else{
          #     n <- gsub("\\*", "", V(ig)$label[i])
          #     paste0("<a target='_blank' ",
          #            "href='http://www.genome.jp/dbget-bin/www_bget?", n, "'>", 
          #            n, "</a>")
          # }
      }
  })
  
  comp.is <- which(V(ig)$type == "compound")
  for(i in comp.is){
      cpd <- gsub("\\*", "", V(ig)$label[i])
      if(cpd %in% comp.list$simple.ID)
          V(ig)$label[i] <- gsub(cpd, comp.list[comp.list$simple.ID == cpd, 2],
                                 V(ig)$label[i])
  }
  V(ig)$name <- paste("N", name.pathway, V(ig)$name, sep="-")

  return(ig)
}


get.hgnc.gene <- function(hgnc, gen){
    trans <- hgnc[hgnc[,1] == gen, 2]
    locs <- grepl("LOC", trans)
    if(sum(!locs) > 0)
        trans <- trans[!locs]
    if (length(trans) > 0){
        tran <- trans[1]
    }else{
        tran <- ""
    }
    return(tran)
}

get.shape <- function(kgraph){
  shapes <- sapply(kgraph@nodeData@defaults$KEGGNode$nodes, function(x){
    x@graphics@type
  })
  return(shapes)
}



# # Returns a list where each entry is the first name of each protein included in each node
# get.name.nodes <- function( kgraph ){
#   nodes <- kgraph@nodeData@defaults$KEGGNode$nodes
#   node.names <- NULL
#   for(i in 1:length(nodes)){
#     if( grepl( "[", nodes[[i]]@graphics@name, fixed = TRUE  )){
#       names.cor <- unlist(strsplit(nodes[[i]]@graphics@name, "[", fixed = TRUE))
#       names.cor <- names.cor[which(names.cor != "")]
#     }else{
#       names.cor <- nodes[[i]]@graphics@name[[1]]
#     }
#     names <- NULL
#     for( j in 1:length(names.cor)){
#       names.cor[[j]] <- gsub("]", "", names.cor[[j]], fixed=TRUE)
#       names.cor[[j]] <- gsub("...", "", names.cor[[j]], fixed=TRUE)
#       if(nodes[[i]]@type == "other"){
#           namesj <- unlist(strsplit(names.cor[[j]], "\\("))
#           names <- paste(names, gsub("\\)", "", namesj[[2]]), sep = " ")
#       }else{
#           namesj <- unlist(strsplit(names.cor[[j]], ","))
#           names <- paste( names, namesj[[1]], sep = " ")
#       }
#     }
#     names <- sub( " ", "", names, fixed = TRUE )
#     while( names %in% node.names )
#       names <- paste( names, "*", sep="")
#     node.names <- c( node.names, names )
#   }
#   return( node.names )
# }


# Returns a list where each entry is the first name of each protein included in each node
get.name.nodes <- function( kgraph ){
    nodes <- kgraph@nodeData@defaults$KEGGNode$nodes
    node.names <- NULL
    for(node in nodes){
        # print(node@graphics@name)
        type <- node@type
        if(grepl("\\,", type)){
            names <- unlist(strsplit(node@graphics@name, "], [", fixed = TRUE))
            names <- gsub("\\[", "", gsub("\\]", "", names))
            types <- unlist(strsplit(type, split = "\\,"))
            ns <- sapply(1:length(names), function(i)
                get.name.node(names[i], types[i]))
            new <- paste(ns, collapse = " ")
        }else{
            new <- get.name.node(node@graphics@name, type)
        }
        while(new %in% node.names)
            new <- paste( new, "*", sep="")
        node.names <- c( node.names, new )
    }
    return( node.names )
}


# Returns a list where each entry is the first name of each protein included in each node
get.name.node <- function(name, type){
    name <- gsub("...", "", name, fixed=TRUE)
    if(type == "compound"){
        comps <- load_compounds(comp.file)
        if(name %in% comps$simple.ID){
            names <- comps[comps$simple.ID == name, 2]
        }else{
            names <- name
        }
    }else if(type == "other"){
        namesj <- unlist(strsplit(name, "\\("))
        names <- gsub("\\)", "", namesj[[2]])
    }else if(type == "gene"){
        namesj <- unlist(strsplit(name, ","))
        names <- namesj[[1]]
    }else{
        names <- name
    }
    return(names)
}



# Returns a list where each entry is the list of genes included in each node
get.genes.list.type <- function( kgraph ){
  nodes <- kgraph@nodeData@defaults$KEGGNode$nodes
  genes <- list()
  for( i in 1: length(nodes)){
    names <- NULL
    types <- unlist(strsplit(nodes[[i]]@type, split=","))
    if("gene" %in% types){
      for( j in 1:length(nodes[[i]]@name)){
        if( nodes[[i]]@name[[j]] == "/" ){
          names <- c(names, "/")
        }else{
          splited <- unlist(strsplit(nodes[[i]]@name[[j]], ":", fixed = TRUE))
          if(splited[[1]] == "cpd" | splited[[1]] == "BR"){
            if(is.null(names) || !is.na(names[length(names)]))
              names <- c(names, NA)
          }else{
            names <- c( names, splited[[2]])
          }
        }
      }
    }else{
      names <- NA
    }
    genes[[length(genes)+1]] <- names
  }
  return( genes )
}



get.species <- function(kgraph){
  i <- 1
  while(unlist(strsplit(kgraph@nodeData@defaults$KEGGNode$nodes[[i]]@name[[1]], ":", fixed = TRUE))[1] == "cpd"){
    i <- i+1
  }
  return(unlist(strsplit(kgraph@nodeData@defaults$KEGGNode$nodes[[i]]@name[[1]], ":", fixed = TRUE))[1])
}


# Returns a list where each entry is the pair of XY coordinates of the node
get.xy.vals <- function( kgraph ){
  nodes <- kgraph@nodeData@defaults$KEGGNode$nodes
  x <- sapply(nodes, function(node){node@graphics@x})
  y <- sapply(nodes, function(node){node@graphics@y})
  coor <- cbind(x, y)
  names(coor) <- get.name.nodes( kgraph )
  return(coor)
}



# Changes the weight of inhibition edges to "-1"
add.weight.info <- function( kgraph ){
  types <- sapply(kgraph@edgeData@defaults$KEGGEdge$edges, function(e){ifelse(!is.null(e@subtype$subtype), e@subtype$subtype@name, "activation")})
  types[which(types=="compound")] <- sapply(which(types=="compound"), function(x){ifelse(length(kgraph@edgeData@defaults$KEGGEdge$edges[[x]]@subtype)>1, kgraph@edgeData@defaults$KEGGEdge$edges[[x]]@subtype[[2]]@name, "activation")})
  weights <- sapply(types, function(x){ifelse(x == "inhibition" || x == "ubiquitination" || x == "methylation" || x == "repression", -1, 1)})
  return(weights)
}





pathigraph2IG <- function(graph){

  sif <- do.call("rbind", lapply(names(graph@edgeData), function(x){
    nodes <- unlist(strsplit(x, "\\|"))
    actin <- ifelse(graph@edgeData@data[[x]]$weight == 1, "activation", "inhibition")
    return(c(nodes[1], actin, nodes[2]))
  }))

  ig <- graph.data.frame(sif[,c(1,3)], directed=T)

  # Add attributes
  E(ig)$relation <- unlist(sapply(sif[,2], function(x){ifelse(x=="activation", 1, -1)}))
  E(ig)$curved <- F
  V(ig)$nodeX <- graph@nodeData@defaults$Coordinates[match(V(ig)$name,graph@nodes),1]
  V(ig)$nodeY <- max(graph@nodeData@defaults$Coordinates[,2]) - graph@nodeData@defaults$Coordinates[match(V(ig)$name,graph@nodes),2]
  V(ig)$shape <- "rectangle"
  V(ig)$label.cex <- 0.5
  V(ig)$label.color <- "black"
  V(ig)$size <- 15
  V(ig)$size2 <- 5
  glist <- graph@nodeData@defaults$GenesList
  V(ig)$genesList <- glist[V(ig)$name]
  V(ig)$label <- graph@nodeData@defaults$NodeNames[match(V(ig)$name,graph@nodes)]

  return(ig)

}



pathigraph2sif <- function(graph, pathway.name){
  sif <- do.call("rbind", lapply(names(graph@edgeData), function(x){
    nodes <- unlist(strsplit(x, "\\|"))
    actin <- ifelse(graph@edgeData@data[[x]]$weight == 1, "activation", "inhibition")
    return(c(nodes[1], actin, nodes[2]))
  }))
  att <- data.frame(cbind(graph@nodes, graph@nodeData@defaults$Coordinates[,1], graph@nodeData@defaults$Coordinates[,2], "rectangle", 0.5, "black", 15, 5, sapply(graph@nodeData@defaults$GenesList, paste, collapse=",")))
  write.table(sif, file=paste("files/sifs/", pathway.name, ".sif", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
  write.table(att, file=paste("files/sifs/", pathway.name, ".att", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
}


delete.isolated.nodes <- function( graph ){
  i <- 1
  while( i <= length( graph@nodeData@defaults$KEGGNode$nodes )){
    nodei <- graph@nodeData@defaults$KEGGNode$nodes[[i]]
    flag <- TRUE
    edges <- graph@edgeData@defaults$KEGGEdge$edges
    for( j in 1:length( edges )){
      if( edges[[j]]@entry1ID == nodei@entryID || edges[[j]]@entry2ID == nodei@entryID )
        flag <- FALSE
    }
    if( flag ){
      graph <- removeNode( nodei@entryID, graph )
      graph@nodeData@defaults$KEGGNode$nodes <- graph@nodeData@defaults$KEGGNode$nodes[-i]
      i <- i-1
    }
    i <- i+1
  }
  return( graph )
}



collapse.nodes <- function( graph ){
    graph <- removeInnecessaryBindings(graph)
  new.graph <- graph
  i <- 1
  while( i <= length( new.graph@nodeData@defaults$KEGGNode$nodes )){
    node1 <- new.graph@nodeData@defaults$KEGGNode$nodes[[i]]
    bind.edges1 <- binding.edges.from.to( node1, new.graph@edgeData@defaults$KEGGEdge$edges )
    while( length(bind.edges1) >= 1 ){
      edges <- new.graph@edgeData@defaults$KEGGEdge$edges
      node2 <- find.opposite( node1@entryID, bind.edges1[[1]], new.graph@nodeData@defaults$KEGGNode$nodes )
      bind.edges2 <- binding.edges.from.to.except( node2, node1, edges )
      flag <- TRUE
      if( length(bind.edges2) == 1 && length( not.binding.edges.from.to( node2, edges )) == 0 ){
        n3 <- find.opposite(node2@entryID, bind.edges2[[1]], new.graph@nodeData@defaults$KEGGNode$nodes)
        if( length( binding.edges.from.to.except(n3, node2, edges)) == 0 && length( not.binding.edges.from.to( n3, edges)) == 0 ){
          # CASE 3 NODES connected as a string with two of them not having any not-binding edge
          node3 <- n3
          flag <- FALSE
        }
      }
      if( flag )
        node3 <- nodo.comun( node1, node2, bind.edges1, bind.edges2, new.graph@nodeData@defaults$KEGGNode$nodes, edges )
      if( !is.null(node3) ){

        # CASE 3 NODES
        bedges1 <- binding.edges.from.to.except2( node1, node2, node3, edges )
        bedges2 <- binding.edges.from.to.except2( node2, node1, node3, edges )
        bedges3 <- binding.edges.from.to.except2( node3, node2, node1, edges )
        if( length(bedges1) > 0 ){
          new.graph <- create.new.node( node1, new.graph )
          new.node1 <- search.node( paste( node1@entryID, "_2", sep="" ), getKEGGnodeData( new.graph ))
        }else{
          new.node1 <- node1
          i <- i-1
        }
        if( length(bedges2) > 0 ){
          new.graph <- create.new.node( node2, new.graph )
          new.node2 <- search.node( paste( node2@entryID, "_2", sep="" ), getKEGGnodeData( new.graph ))
        }else{
          new.node2 <- node2
        }
        if( length(bedges3) > 0 ){
          new.graph <- create.new.node( node3, new.graph )
          new.node3 <- search.node( paste( node3@entryID, "_2", sep="" ), getKEGGnodeData( new.graph ))
        }else{
          new.node3 <- node3
        }
        new.graph <- delete.redundant.edges( node1, node2, new.graph )
        new.graph <- delete.redundant.edges( node1, node3, new.graph )
        new.graph <- delete.redundant.edges( node2, node3, new.graph )
        bind.edges1 <- delete.processed.edges( node1, node2, bind.edges1 )
        bind.edges1 <- delete.processed.edges( node1, node3, bind.edges1 )
        bind.edges1 <- delete.processed.edges( node2, node3, bind.edges1 )
        # Combine the nodes
        new.graph <- combineNodes( c(new.node1@entryID, new.node2@entryID, new.node3@entryID), new.graph, paste( node1@entryID, node2@entryID, node3@entryID, sep=" "))

        # Add KEGGNode information related to new node
        new.kegg.node <- create.collapsed.KEGGnode(c(new.node1,new.node2,new.node3))
        new.graph@nodeData@defaults$KEGGNode$nodes <- c( new.graph@nodeData@defaults$KEGGNode$nodes, new.kegg.node )
        names(new.graph@nodeData@defaults$KEGGNode$nodes)[[length(new.graph@nodeData@defaults$KEGGNode$nodes)]] <- paste( node1@entryID, node2@entryID, node3@entryID, sep=" ")

        # Remove duplicated nodes
        nombres <- names(new.graph@nodeData@defaults$KEGGNode$nodes)
        to.delete <- !nombres[]==new.node1@entryID & !nombres[]==new.node2@entryID & !nombres[]==new.node3@entryID
        new.graph@nodeData@defaults$KEGGNode$nodes <- new.graph@nodeData@defaults$KEGGNode$nodes[to.delete]

        # Update KEGGEdge information related to new node
        n <- unlist(sapply(names(new.graph@edgeL), function(x){if(length(new.graph@nodes[new.graph@edgeL[[x]]$edges])>0){paste0(x, "|", new.graph@nodes[new.graph@edgeL[[x]]$edges])}}))
        new.graph@edgeData@data <- new.graph@edgeData@data[n]
        pedges <- new.graph@edgeData@defaults$KEGGEdge$edges
        edgeEntry <- sapply(pedges, getEntryID )
        asID1 <- edgeEntry[1L,]==new.node1@entryID | edgeEntry[1L,]==new.node2@entryID | edgeEntry[1L,]==new.node3@entryID
        asID2 <- edgeEntry[2L,]==new.node1@entryID | edgeEntry[2L,]==new.node2@entryID | edgeEntry[2L,]==new.node3@entryID
        names.edges <- names(pedges)
        for( nam in names.edges[asID1] ){
          e <- pedges[[nam]]
          e@entry1ID <- new.kegg.node@entryID
          pedges[[nam]] <- e
        }
        for( nam in names.edges[asID2] ){
          e <- pedges[[nam]]
          e@entry2ID <- new.kegg.node@entryID
          pedges[[nam]] <- e
        }
        names.edges[asID1] <- sapply(names.edges[asID1], function(x, name) paste(name, unlist(strsplit(x, "~", fixed = TRUE))[2], sep="~") , new.kegg.node@entryID)
        names.edges[asID2] <- sapply(names.edges[asID2], function(x, name) paste(unlist(strsplit(x, "~", fixed = TRUE))[1], name, sep="~") , new.kegg.node@entryID)
        names(pedges) <- names.edges
        pedges <- pedges[gsub("\\|", "~", n)]
        new.graph@edgeData@defaults$KEGGEdge$edges <- pedges

        # CASE 2 NODES
      }else{
        bedges1 <- binding.edges.from.to.except( node1, node2, edges )
        bedges2 <- binding.edges.from.to.except( node2, node1, edges )
        if( length(bedges1) > 0 ){
          new.graph <- create.new.node( node1, new.graph )
          new.node1 <- search.node( paste( node1@entryID, "_2", sep="" ), getKEGGnodeData( new.graph ))
        }else{
          new.node1 <- node1
          i <- i-1
        }
        if( length(bedges2) > 0 ){
          new.graph <- create.new.node( node2, new.graph )
          new.node2 <- search.node( paste( node2@entryID, "_2", sep="" ), getKEGGnodeData( new.graph ))
        }else{
          new.node2 <- node2
        }
        new.graph <- delete.redundant.edges( node1, node2, new.graph )
        bind.edges1 <- delete.processed.edges( node1, node2, bind.edges1 )
        # Combine the two nodes in just one
        new.graph <- combineNodes( c(new.node1@entryID, new.node2@entryID ), new.graph, paste( node1@entryID, node2@entryID, sep=" " ))
        # Add KEGGNode information related to new node
        new.kegg.node <- create.collapsed.KEGGnode(c(new.node1,new.node2))
        new.graph@nodeData@defaults$KEGGNode$nodes <- c( new.graph@nodeData@defaults$KEGGNode$nodes, new.kegg.node )
        names(new.graph@nodeData@defaults$KEGGNode$nodes)[[length(new.graph@nodeData@defaults$KEGGNode$nodes)]] <- paste( node1@entryID, node2@entryID, sep=" ")
        # Remove duplicated nodes
        nombres <- names(new.graph@nodeData@defaults$KEGGNode$nodes)
        to.delete <- !nombres[]==new.node1@entryID & !nombres[]==new.node2@entryID
        new.graph@nodeData@defaults$KEGGNode$nodes <- new.graph@nodeData@defaults$KEGGNode$nodes[to.delete]

        # Update KEGGEdge information related to new node
        n <- unlist(sapply(names(new.graph@edgeL), function(x){if(length(new.graph@nodes[new.graph@edgeL[[x]]$edges])>0){paste0(x, "|", new.graph@nodes[new.graph@edgeL[[x]]$edges])}}))
        new.graph@edgeData@data <- new.graph@edgeData@data[n]
        pedges <- new.graph@edgeData@defaults$KEGGEdge$edges
        edgeEntry <- sapply(pedges, getEntryID )
        asID1 <- edgeEntry[1L,]==new.node1@entryID | edgeEntry[1L,]==new.node2@entryID
        asID2 <- edgeEntry[2L,]==new.node1@entryID | edgeEntry[2L,]==new.node2@entryID
        names.edges <- names(pedges)
        for( nam in names.edges[asID1] ){
          e <- pedges[[nam]]
          e@entry1ID <- new.kegg.node@entryID
          pedges[[nam]] <- e
        }
        for( nam in names.edges[asID2] ){
          e <- pedges[[nam]]
          e@entry2ID <- new.kegg.node@entryID
          pedges[[nam]] <- e
        }
        names.edges[asID1] <- sapply(names.edges[asID1], function(x, name) paste(name, unlist(strsplit(x, "~", fixed = TRUE))[2], sep="~") , new.kegg.node@entryID)
        names.edges[asID2] <- sapply(names.edges[asID2], function(x, name) paste(unlist(strsplit(x, "~", fixed = TRUE))[1], name, sep="~") , new.kegg.node@entryID)
        names(pedges) <- names.edges
        pedges <- pedges[gsub("\\|", "~", n)]
        new.graph@edgeData@defaults$KEGGEdge$edges <- pedges
      }
    }
    i <- i+1
  }
  return(new.graph)
}

create.collapsed.KEGGnode <- function( node.list ){
  new.ID <- NULL
  new.name <- NULL
  new.type <- NULL
  new.link <- NULL
  new.reaction <- NULL
  new.map <- NULL
  new.graphics <- NULL

  for( node in node.list ){
    if( grepl("_", node@entryID, fixed=TRUE))
      node@entryID <- unlist(strsplit(node@entryID, "_", fixed = TRUE))[1]
    if( is.null(new.ID)){
      new.ID <- paste( new.ID, node@entryID, sep="")
    }
    else{
      new.ID <- paste( new.ID, node@entryID, sep=" ")
    }
    if( is.null(new.name)){
      new.name <- c(new.name, node@name)
    }
    else{
      new.name <- c( new.name, "/", node@name )
    }
    # TODO: Probar qué pasa si ponemos type "group"
    if( is.null( new.type)){
      new.type <- node@type
    }
    else{
      new.type <- paste(new.type, node@type, sep=",")
    }
    if( is.null(new.link)){
      new.link <- node@link
    }
    else{
      for( k in 1:length(node@name))
        new.link <- paste( new.link, "+", node@name[[k]], sep="")
    }
    if( is.null(new.reaction)){
      new.reaction <- node@reaction
    }
    else{
      new.reaction <- paste(new.reaction, node@reaction, sep=",")
    }
    if( is.null(new.map)){
      new.map <- node@map
    }
    else{
      new.map <- paste(new.map, node@map, sep=",")
    }
    if(is.null(new.graphics)){
      new.graphics <- node@graphics
      new.graphics@name <- paste("[", new.graphics@name, "]", sep="")
      new.graphics@x <- as.integer(new.graphics@x )
      new.graphics@y <- as.integer(new.graphics@y )
    }
    else{
      display.name2 <- sub("...", "", node@graphics@name, fixed=TRUE)
      display.name1 <- new.graphics@name
      new.graphics@name <- paste( display.name1, ", [", display.name2, "]", sep="")
    }
  }

  new.kegg.node <- new("KEGGNode", entryID = new.ID, name = new.name, type = new.type, link = new.link, reaction = new.reaction, map = new.map, graphics = new.graphics )
  return(new.kegg.node)
}


# Returns the list of binding edges, without the binding edge between "node1" and "node2"
delete.processed.edges <- function( node1, node2, bind.edges ){
  i <- 1
  while( i <= length(bind.edges) ){
    edge <- bind.edges[[i]]
    if(( edge@entry1ID == node1@entryID && edge@entry2ID == node2@entryID ) || ( edge@entry1ID == node2@entryID && edge@entry2ID == node1@entryID )){
      bind.edges <- bind.edges[-i]
      i <- i-1
    }
    i <- i+1
  }
  return( bind.edges )
}


# Returns TRUE if "edge" is included in "list", FALSE otherwise
contains.kegg.edge <- function( list, edge ){
  for( e in list )
    if( e@entry1ID == edge@entry1ID && e@entry2ID == edge@entry2ID )
      return(TRUE)
  return(FALSE)
}

# Returns a graph from where KEGGEdges from "node1" to "node2" or viceversa have been removed
delete.redundant.edges <- function( node1, node2, graph ){
  new.graph <- graph
  i <- 1
  while( i <= length(new.graph@edgeData@defaults$KEGGEdge$edges)){
    edgei <- new.graph@edgeData@defaults$KEGGEdge$edges[[i]]
    if(( edgei@entry1ID == node1@entryID && edgei@entry2ID == node2@entryID ) || (edgei@entry2ID == node1@entryID && edgei@entry1ID == node2@entryID )){
        new.graph <- removeEdge( edgei@entry1ID, edgei@entry2ID, new.graph )
        new.graph@edgeData@defaults$KEGGEdge$edges <- new.graph@edgeData@defaults$KEGGEdge$edges[-i]
        i <- i-1
    }
    i <- i+1
  }
  return(new.graph)
}


nodo.comun <- function( node1, node2, edges1, edges2, vertices, edges ){
  for( edge1 in edges1 ){
    node3 <- find.opposite( node1@entryID, edge1, vertices )
    if( node3@entryID != node2@entryID ){
      for( edge2 in edges2 ){
        node4 <- find.opposite( node2@entryID, edge2, vertices )
        if( node4@entryID != node1@entryID && node4@entryID == node3@entryID )
          return( node4 )
      }
    }
  }
  return(NULL)
}

#node <- node2

create.new.node <- function( node, graph ){
  # Create node data
  new.kegg.node <- new("KEGGNode", entryID = paste( node@entryID, "_2", sep="" ), name = node@name, type = node@type, link = node@link, reaction = node@reaction, map = node@map, graphics = node@graphics )
  new.graph <- graph
  new.graph <- graph::addNode( paste( node@entryID, "_2", sep="" ), graph )
  new.graph@nodeData@defaults$KEGGNode$nodes <- c( new.graph@nodeData@defaults$KEGGNode$nodes, new.kegg.node )
  names(new.graph@nodeData@defaults$KEGGNode$nodes)[[length(new.graph@nodeData@defaults$KEGGNode$nodes)]] <- paste( node@entryID, "_2", sep="" )
  # Create Edge data
  edges <- getKEGGedgeData( new.graph )
  for( i in 1:length(edges)){
    if( edges[[i]]@entry1ID == node@entryID && !is.null(edges[[i]]@subtype$subtype) && edges[[i]]@subtype$subtype@name != "binding/association"){
      new.kegg.edge <- new( "KEGGEdge", entry1ID = paste( node@entryID, "_2", sep="" ), entry2ID = edges[[i]]@entry2ID, type = edges[[i]]@type, subtype = edges[[i]]@subtype )
      new.graph <- addEdge( paste( node@entryID, "_2", sep="" ), edges[[i]]@entry2ID, new.graph )
      new.graph@edgeData@defaults$KEGGEdge$edges <- c( new.graph@edgeData@defaults$KEGGEdge$edges, new.kegg.edge )
      names(new.graph@edgeData@defaults$KEGGEdge$edges)[[length(new.graph@edgeData@defaults$KEGGEdge$edges)]] <- paste( paste( node@entryID, "_2", sep="" ), "~", edges[[i]]@entry2ID, sep="" )
    }else if( edges[[i]]@entry2ID == node@entryID && !is.null(edges[[i]]@subtype$subtype) && edges[[i]]@subtype$subtype@name != "binding/association"){
      new.kegg.edge <- new( "KEGGEdge", entry1ID = edges[[i]]@entry1ID, entry2ID = paste( node@entryID, "_2", sep="" ), type = edges[[i]]@type, subtype = edges[[i]]@subtype )
      new.graph <- addEdge( edges[[i]]@entry1ID, paste( node@entryID, "_2", sep="" ), new.graph )
      new.graph@edgeData@defaults$KEGGEdge$edges <- c( new.graph@edgeData@defaults$KEGGEdge$edges, new.kegg.edge )
      names(new.graph@edgeData@defaults$KEGGEdge$edges)[[length(new.graph@edgeData@defaults$KEGGEdge$edges)]] <- paste( edges[[i]]@entry1ID, "~", paste( node@entryID, "_2", sep="" ), sep="" )
    }
  }
  return( new.graph )
}


find.opposite <- function( name.node, edge, vertices ){
  if( edge@entry1ID == name.node ){
    node = search.node( edge@entry2ID, vertices )
  }else{
    node = search.node( edge@entry1ID, vertices)
  }
  return(node)
}

delete.unwanted.edges <- function( graph ){
  vertices <- graph@nodeData@defaults$KEGGNode$nodes
  i <- 1
  while( i <= length(graph@edgeData@defaults$KEGGEdge$edges)){
    edgei <- graph@edgeData@defaults$KEGGEdge$edges[[i]]
    if( !is.null(edgei@subtype$subtype) && edgei@subtype$subtype@name == "binding/association" ){
      iflag <- FALSE
      node1 <- search.node( edgei@entry1ID, vertices )
      node2 <- search.node( edgei@entry2ID, vertices )
      edges1 <- binding.edges.from.to.except(node1, node2, graph@edgeData@defaults$KEGGEdge$edges)
      if( length(edges1) > 0 ){
        for( j in 1:length(edges1)){
          # define node3  as the node different to node1 linked by edges1[[j]]
          if( edges1[[j]]@entry1ID == node1@entryID ){
            node3 = search.node( edges1[[j]]@entry2ID, vertices )
          }else{
            node3 = search.node( edges1[[j]]@entry1ID, vertices )
          }
          edgeIDX <- exist.edge(graph@edgeData@defaults$KEGGEdge$edges, node3, node2)
          if( !is.null(edgeIDX) && graph@edgeData@defaults$KEGGEdge$edges[[edgeIDX]]@subtype$subtype@name != "binding/association" ){
            iflag <- TRUE
          }
        }
      }

      edges2 <- binding.edges.from.to.except(node2, node1, graph@edgeData@defaults$KEGGEdge$edges)
      if( length(edges2) > 0 ){
        for( j in 1:length(edges2)){
          # define node3  as the node different to node2 linked by edges2[[j]]
          if( edges2[[j]]@entry1ID == node2@entryID ){
            node3 = search.node( edges2[[j]]@entry2ID, vertices )
          }else{
            node3 = search.node( edges2[[j]]@entry1ID, vertices)
          }
          edgeIDX <- exist.edge(graph@edgeData@defaults$KEGGEdge$edges, node3, node1)
          if( !is.null(edgeIDX) && graph@edgeData@defaults$KEGGEdge$edges[[edgeIDX]]@subtype$subtype@name != "binding/association" ){
            iflag <- TRUE
          }
        }
      }

      if( iflag ){
        graph <- removeEdge(edgei@entry1ID, edgei@entry2ID, graph)
        graph@edgeData@defaults$KEGGEdge$edges <- graph@edgeData@defaults$KEGGEdge$edges[-i]
        i <- i-1
      }
    }
    i <- i+1
  }
  return(graph)
}


exist.edge <- function( edges, node1, node2 ){
  for( i in 1:length(edges))
    if ( edges[[i]]@entry1ID == node1@entryID && edges[[i]]@entry2ID == node2@entryID)
      return(i)
  return(NULL)

}

binding.edges.from.to.except <- function( node1, node2, edges ){
  edges1 <- NULL
  for( i in 1:length(edges)){
    if( !is.null(edges[[i]]@subtype$subtype) && edges[[i]]@subtype$subtype@name == "binding/association" && ((edges[[i]]@entry1ID == node1@entryID && edges[[i]]@entry2ID != node2@entryID) || (edges[[i]]@entry2ID == node1@entryID && edges[[i]]@entry1ID != node2@entryID)))
      edges1 <- c( edges1, edges[[i]])
  }
  return(edges1)
}


binding.edges.from.to.except2 <- function( node1, node2, node3, edges ){
  edges1 <- NULL
  for( i in 1:length(edges)){
    if( !is.null(edges[[i]]@subtype$subtype) && edges[[i]]@subtype$subtype@name == "binding/association" && ((edges[[i]]@entry1ID == node1@entryID && edges[[i]]@entry2ID != node2@entryID && edges[[i]]@entry2ID != node3@entryID) || (edges[[i]]@entry2ID == node1@entryID && edges[[i]]@entry1ID != node2@entryID && edges[[i]]@entry1ID != node3@entryID)))
      edges1 <- c( edges1, edges[[i]])
  }
  return(edges1)
}

# Returns the list of edges of type "binding/association" arriving to or leaving from "node"
binding.edges.from.to <- function( node, edges ){
  edges1 <- NULL
  for( i in 1:length(edges)){
    if( !is.null(edges[[i]]@subtype$subtype) && edges[[i]]@subtype$subtype@name == "binding/association" && ( edges[[i]]@entry1ID == node@entryID || edges[[i]]@entry2ID == node@entryID ))
      edges1 <- c( edges1, edges[[i]])
  }
  return(edges1)
}

# Returns the list of edges of type different to "binding/association" arriving to or leaving from "node"
not.binding.edges.from.to <- function( node, edges ){
  edges1 <- NULL
  for( i in 1:length(edges)){
    if( !is.null(edges[[i]]@subtype$subtype) && edges[[i]]@subtype$subtype@name != "binding/association" && ( edges[[i]]@entry1ID == node@entryID || edges[[i]]@entry2ID == node@entryID ))
      edges1 <- c( edges1, edges[[i]])
  }
  return(edges1)
}


search.node <- function(name.node, vertices){
  if( length( vertices) > 0 ){
    for( i in 1:length(vertices)){
      if( contains.name.node( vertices[[i]], name.node ))
        return(vertices[[i]])
    }
  }
  return(NULL)
}


contains.name.node <- function( vertex, name.node ){
  name.vertex <- vertex@entryID
  if( name.vertex == name.node ){
    return(TRUE)}
  else {
    names.list <- unlist(strsplit(name.vertex, " "))
    for( i in 1:length(names.list))
      if( names.list[i] == name.node )
        return(TRUE)
  }
  return(FALSE)
}


removeInnecessaryBindings <- function(graph, verbose = TRUE){
    es <- graph@edgeData@defaults$KEGGEdge$edges
    types <- sapply(es, function(x) if(!is.null(x@subtype$subtype)){x@subtype$subtype@name}else{x@type})
    nodes1 <- as.character(sapply(es, function(x) x@entry1ID))
    nodes2 <- as.character(sapply(es, function(x) x@entry2ID))
    all.nodes <- unique(c(nodes1, nodes2))
    all.nodes.unlist <- lapply(all.nodes, function(x) unlist(strsplit(x, split = " ")))
    es.mat <- data.frame(nodes1, nodes2, types, stringsAsFactors = F)
    es.mat$exist <- apply(es.mat[,c(1,2)], 1, function(v){
        v1 <- unlist(strsplit(v[1], split = " "))
        v2 <- unlist(strsplit(v[2], split = " "))
        both <- unique(c(v1, v2))
        any(sapply(all.nodes.unlist, function(x) all(both %in% x)))
    })
    if(any(es.mat$exist)){
        toremove <- which(es.mat$exist == TRUE & es.mat$type == "binding/association")
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
