#' Identify sister taxa from node in a phylogeny
#' @importFrom ape extract.clade
#' @export
sister_taxa <- function(tree, node){
  mat.edge <- as.matrix(tree$edge)
  sister.clades <- mat.edge[mat.edge[,1] == node, ]
  if(is.na(tree$tip.label[sister.clades[1,2]]) & is.na(tree$tip.label[sister.clades[2,2]])){
    sister.names.1 <- extract.clade(tree, sister.clades[1,2])$tip.label
    sister.names.2 <- extract.clade(tree, sister.clades[2,2])$tip.label
    sister.df <- data.frame(ancestral.node = rep(node, (length(sister.names.1)+length(sister.names.2))),
                            sister.clade = c(rep(1, length(sister.names.1)), rep(2, length(sister.names.2))),
                            taxa = c(sister.names.1, sister.names.2))}
  if(is.na(tree$tip.label[sister.clades[1,2]]) & is.character(ifelse(is.na(tree$tip.label[sister.clades[2,2]]), 0, tree$tip.label[sister.clades[2,2]]))){
    sister.names.1 <- extract.clade(tree, sister.clades[1,2])$tip.label
    sister.taxon.2 <- tree$tip.label[sister.clades[2,2]]
    sister.df <- data.frame(ancestral.node = rep(node, length(sister.names.1)+1),
                            sister.clade = c(rep(1, length(sister.names.1)), 2),
                            taxa = c(sister.names.1, sister.taxon.2))}
  if(is.character(ifelse(is.na(tree$tip.label[sister.clades[1,2]]), 0, tree$tip.label[sister.clades[1,2]])) & is.na(tree$tip.label[sister.clades[2,2]])){
    sister.taxon.1 <- tree$tip.label[sister.clades[1,2]]
    sister.names.2 <- extract.clade(tree, sister.clades[2,2])$tip.label
    sister.df <- data.frame(ancestral.node = rep(node, length(sister.names.2)+1),
                            sister.clade = c(1, rep(2, length(sister.names.2))),
                            taxa = c(sister.taxon.1, sister.names.2))}
  if(is.character(ifelse(is.na(tree$tip.label[sister.clades[1,2]]), 0, tree$tip.label[sister.clades[1,2]])) & is.character(ifelse(is.na(tree$tip.label[sister.clades[2,2]]), 0, tree$tip.label[sister.clades[2,2]]))){
    sister.taxon.1 <- tree$tip.label[sister.clades[1,2]]
    sister.taxon.2 <- tree$tip.label[sister.clades[2,2]]
    sister.df <- data.frame(ancestral.node = rep(node, 2), sister.clade = c(1,2),
                            taxa = c(sister.taxon.1, sister.taxon.2))}
  return(sister.df)
}
