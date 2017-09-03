#' Modify lineage names from sequence alignment so they match names in database
#' @importFrom dplyr mutate
#' @importFrom dplyr recode
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export

clean_names <- function(lin.names, keep.genus = FALSE){
  out <- data.frame(do.call("rbind", strsplit(as.character(lin.names), "_")))
  if(keep.genus == FALSE){
    return(as.character(out[,2]))
  } else{
    out.n <- out %>% select(X1, X2) %>% rename(parasiteGenus = X1, Lineage_Name = X2) %>%
      mutate(parasiteGenus = recode(parasiteGenus, P = "Plasmodium", H = "Haemoproteus",
                                    L = "Leucocytozoon"))
    return(out.n)
  }
}
