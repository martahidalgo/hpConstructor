
# 
# pathways_per_cathegory <- function(species, cats = NULL, subcats = NULL, 
#                                    paths = NULL){
#     df <- load_cats(species) ##########################
#     ids <- NULL
#     if(!is.null(cats))
#         ids <- union(ids, df[df$Cathegory %in% cats, "Pathway.id"])
#     if(!is.null(subcats))
#         ids <- union(ids, df[df$Subcathegor %in% subcats, "Pathway.id"])
#     if(!is.null(paths))
#         ids <- union(ids, df[df$Pathway.name %in% paths, "Pathway.id"])
#     
#     return(paste0(species, ids))
# }
# 
# 
# available_pathways<- function(species, paths = TRUE, cats = FALSE, 
#                               subcats = FALSE){
#     nfeat <- sum(c(paths, cats, subcats))
#     if(nfeat > 1)
#         stop("Please, select only one feature: paths, cats or subcats")
#     if(nfeat < 1)
#         stop("Please, select one feature: paths, cats or subcats")
#     
#     df <- load_cats(species) ##########################
#     metaginfo <- load_mgi(species)
#     all_paths <- unique(metaginfo$all.labelids[,"path.name"])
#     df <- df[df$Pathway.name %in% all_paths,]
#     if(paths == TRUE){
#         return(all_paths)
#     }else if(subcats == TRUE){
#         return(unique(df$Subcathegory))
#     }else if(cats == TRUE){
#         return(unique(df$Cathegory))
#     }
# }


load_cats <- function(){
    file <- "~/appl/hpConstructor/code_develop/pathway_categories/pathway_cathegories.RData"
    cats_df <- load(file, envir = environment())
    cats <- get(cats_df)
    return(cats)
}

load_index <- function(species){
    file <- "~/appl/hpConstructor/code_develop/pathway_categories/pathway_index.RData"
    index_df <- load(file, envir = environment())
    index <- get(index_df)
    return(index[[species]])
}

available_cathegories <- function(species, subs = FALSE){
    df <- get_pathway_cathegories(species)
    if(subs == TRUE){
        return(unique(df[,3]))
    }else{
        return(unique(df[,4]))
    }
}


get_pathway_cathegories <- function(species){
    cats <- load_cats()
    index <- load_index(species)
    index$Pathway.id <- gsub(species, "", index$path.id)
    return(cats[cats$Pathway.id %in% index$Pathway.id,])
}


get_pathways_from_cathegories <- function(species, cat_list = NULL, sub_list = NULL){
    if(is.null(cat_list) & is.null(sub_list))
        stop("At least one of the parameters (cat_list, sub_list) must be provided")
    cats <- get_pathway_cathegories(species)
    df <- cats[cats$Subcathegory %in% sub_list | cats$Cathegory %in% cat_list,]
    return(df)
}


load_pathways_per_cathegory <- function(species, cat_list = NULL, 
                                        sub_list = NULL, path_list = NULL){
    if(is.null(cat_list) & is.null(sub_list))
        stop("At least one of the parameters (cat_list, sub_list) must be provided")
    df <- get_pathways_from_cathegories(species, cat_list, sub_list)
    paths <- c(path_list, paste0(species, df$Pathway.id))
    mgi <- load_pathways(species, paths)
    return(mgi)
}
