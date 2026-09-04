#' Identify the sister taxa at a node in a phylogeny
#'
#' For an internal node, returns the tips descending from each of its immediate
#' descendant clades, labeled as sister clade 1, 2, and so on. This is useful,
#' for example, for comparing the hosts or traits of sister lineages in a
#' parasite phylogeny (Ellis and Bensch 2018). One or several nodes may be
#' supplied.
#'
#' A bifurcating node -- the usual case -- gives two sister clades. At a polytomy
#' every descendant is sister to every other, so a node with three or more
#' children gives that many groups rather than two. (Before version 1.1.2 only
#' the first two were returned, with nothing to say the rest had been dropped.)
#' A node with a single descendant has no sister groups and is an error.
#'
#' @param tree A phylogeny of class \code{phylo} (see \pkg{ape}).
#' @param node An internal node number, or a vector of node numbers. For a vector,
#'   results for each node are stacked into one data frame. Each value must be an
#'   internal (non-tip) node of \code{tree}; a tip number or an out-of-range value
#'   is an error.
#' @return A \code{data.frame} with columns \code{ancestral.node} (the node
#'   supplied), \code{sister.clade} (1, 2, ... labeling the node's immediate
#'   descendant clades in edge order), and \code{taxa} (tip label). One row per
#'   descending tip.
#' @references
#' Ellis VA, Bensch S (2018). Host specificity of avian haemosporidian parasites
#' is unrelated among sister lineages but shows phylogenetic signal across larger
#' clades. International Journal for Parasitology 48: 897-902.
#' \doi{10.1016/j.ijpara.2018.05.005}
#' @examples
#' tree <- ape::read.tree(text = "((A,B),(C,(D,E)));")
#' sister_taxa(tree, node = 8)
#'
#' ## at a polytomy, every descendant is returned as its own sister group
#' poly <- ape::read.tree(text = "((A,B),(C,D,E));")
#' sister_taxa(poly, node = 8)
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

  ## Descendants of one node, as a data frame of one row per tip.
  ##
  ## Each child edge of the node defines one sister group: a tip contributes
  ## itself, an internal node contributes every tip in its clade. The groups are
  ## numbered in edge order, so on a bifurcating node this is the familiar
  ## sister.clade 1 and 2.
  ##
  ## Written as a loop over the children rather than as separate tip/clade cases.
  ## Before version 1.1.2 the two children were read out of the edge matrix by
  ## position -- sister.clades[1, 2] and [2, 2] -- in four near-identical blocks
  ## covering clade+clade, clade+tip, tip+clade and tip+tip. That silently
  ## discarded the third and later children of a polytomy: on ((A,B),(C,D,E)) the
  ## node above C, D and E returned C and D only, with no indication that E had
  ## been dropped. A node with a single child made the subscript itself invalid
  ## and produced "incorrect number of dimensions".
  sister_taxa_internal <- function(tree, node){
    children <- tree$edge[tree$edge[, 1] == node, 2]
    if (length(children) < 2L) {
      stop("Node ", node, " has ", length(children),
           " descendant edge(s), so it has no sister groups.", call. = FALSE)
    }
    ## a child number <= Ntip is a tip; anything larger is an internal node
    n_tips <- length(tree$tip.label)
    groups <- lapply(children, function(child) {
      if (child <= n_tips) tree$tip.label[child] else extract.clade(tree, child)$tip.label
    })
    data.frame(
      ancestral.node = rep(node, sum(lengths(groups))),
      sister.clade   = rep(seq_along(groups), lengths(groups)),
      taxa           = unlist(groups, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  }

  if(length(node)==1){
    out <- sister_taxa_internal(tree = tree, node = node)
    return(out)
  } else {
    out <- do.call("rbind", lapply(node, function(x)sister_taxa_internal(tree = tree, node = x)))
    return(out)
  }
}
