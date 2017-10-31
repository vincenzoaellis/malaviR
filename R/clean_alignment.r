#' Identify and remove repeated haplotypes from (MalAvi) sequence alignment
#' @importFrom pegas haplotype
#' @importFrom stringi stri_list2matrix
#' @importFrom dplyr group_by
#' @importFrom dplyr sample_n
#' @importFrom dplyr filter
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export

clean_alignment <- function(alignment, separate_by_genus = FALSE, haplotype_format_wide = TRUE){
  if(!class(alignment) == "DNAbin") stop("The alignment should be of class 'DNAbin'.")

  ## identify haplotypes
  h <- haplotype(alignment) # identify unique haplotypes
  h.rep <- subset(h, minfreq = 2) # select haplotypes with more than one lineage
  if(length(h.rep) == 0) stop("The alignment has no repeated haplotypes")

  ## get lineage names for haplotypes representing more than one lineage
  seq.rep <- list() # get the lineage names for each repeated haplotype
  for(i in 1:dim(h.rep)[1]){
    seq.rep[[i]] <- rownames(alignment)[attr(h.rep, "index")[[i]]]
  }

  ## data frame of haplotypes and associated lineage names in wide format
  seq.rep.df <- stri_list2matrix(seq.rep, byrow=TRUE) %>%
    as.data.frame %>%
    cbind(haplotype = 1:dim(h.rep)[1])

  ## data frame of haplotypes and associated lineage names in long format
  seq.rep.df.g <- seq.rep.df %>% gather(lin_number, Lineage_Name, -haplotype) %>%
    filter(!is.na(Lineage_Name))

  ## select one lineage per haplotype at random
  seq.selected <- seq.rep.df.g %>%
    group_by(haplotype) %>%
    sample_n(1)

  ## drop lineages that were not selected from the alignment
  lins.to.drop <- setdiff(seq.rep.df.g$Lineage_Name, seq.selected$Lineage_Name)
  alignment.n <- alignment[!rownames(alignment)%in%lins.to.drop, ]

  if(separate_by_genus == TRUE){
    alignment.plas <- alignment.n[grep("P_", rownames(alignment.n)), ]
    alignment.haem <- alignment.n[grep("H_", rownames(alignment.n)), ]
    alignment.leuc <- alignment.n[grep("L_", rownames(alignment.n)), ]

    if(haplotype_format_wide == TRUE){
      list.out <- list(repeated_haplotypes = seq.rep.df, selected_lineages = seq.selected$Lineage_Name,
                       alignment_clean_Plasmodium = alignment.plas,
                       alignment_clean_Haemoproteus = alignment.haem,
                       alignment_clean_Leucocytozoon = alignment.leuc)
      return(list.out)
    } else{
      list.out <- list(repeated_haplotypes = seq.rep.df.g, selected_lineages = seq.selected$Lineage_Name,
                       alignment_clean_Plasmodium = alignment.plas,
                       alignment_clean_Haemoproteus = alignment.haem,
                       alignment_clean_Leucocytozoon = alignment.leuc)
      return(list.out)
    }
  } else{
    if(haplotype_format_wide == TRUE){
      list.out <- list(repeated_haplotypes = seq.rep.df, selected_lineages = seq.selected$Lineage_Name,
                       alignment_clean = alignment.n)
      return(list.out)
    } else{
      list.out <- list(repeated_haplotypes = seq.rep.df.g, selected_lineages = seq.selected$Lineage_Name,
                       alignment_clean = alignment.n)
      return(list.out)
    }
  }
}
