#' Identify the sister taxa at a node in a phylogeny
#'
#' For an internal node, returns the tips descending from each of its two
#' immediate descendant clades, labelled as sister clade 1 or 2. This is useful,
#' for example, for comparing the hosts or traits of sister lineages in a
#' parasite phylogeny (Ellis and Bensch 2018). One or several nodes may be
#' supplied.
#'
#' @param tree A phylogeny of class \code{phylo} (see \pkg{ape}).
#' @param node An internal node number, or a vector of node numbers. For a vector,
#'   results for each node are stacked into one data frame. Each value must be an
#'   internal (non-tip) node of \code{tree}; a tip number or an out-of-range value
#'   is an error.
#' @return A \code{data.frame} with columns \code{ancestral.node} (the node
#'   supplied), \code{sister.clade} (1 or 2, labelling the node's two immediate
#'   descendant clades), and \code{taxa} (tip label). One row per descending tip.
#' @references
#' Ellis VA, Bensch S (2018). Host specificity of avian haemosporidian parasites
#' is unrelated among sister lineages but shows phylogenetic signal across larger
#' clades. International Journal for Parasitology 48: 897-902.
#' \doi{10.1016/j.ijpara.2018.05.005}
#' @examples
#' tree <- ape::read.tree(text = "((A,B),(C,(D,E)));")
#' sister_taxa(tree, node = 8)
#' @importFrom ape extract.clade
#' @export
sister_taxa <- function(tree, node){
  ## --- validate inputs so bad calls fail with a clear message rather than a
  ## cryptic subscript error deep inside the internal helper ---
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylogeny of class 'phylo' (see ape::read.tree).",
         call. = FALSE)
  }
  if (missing(node) || !is.numeric(node) || length(node) < 1L || any(is.na(node))) {
    stop("`node` must be one or more internal node numbers.", call. = FALSE)
  }
  ## internal nodes are exactly the parent nodes in the edge matrix (tips never
  ## appear in column 1); anything else has no descendant clades to return
  internal_nodes <- unique(tree$edge[, 1])
  bad <- setdiff(node, internal_nodes)
  if (length(bad) > 0) {
    stop("Not an internal node: ", paste(bad, collapse = ", "),
         ". Internal nodes of this tree are ", min(internal_nodes), "-",
         max(internal_nodes), ".", call. = FALSE)
  }

  sister_taxa_internal <- function(tree, node){
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
  if(length(node)==1){
    out <- sister_taxa_internal(tree = tree, node = node)
    return(out)
  } else {
    out <- do.call("rbind", lapply(node, function(x)sister_taxa_internal(tree = tree, node = x)))
    return(out)
  }
}
