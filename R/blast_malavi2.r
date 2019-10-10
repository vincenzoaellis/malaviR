#' BLAST to MalAvi
#' @importFrom magrittr %>%
#' @importFrom rvest set_values
#' @importFrom rvest html_session
#' @importFrom rvest html_form
#' @importFrom rvest submit_form
#' @importFrom httr content
#' @importFrom stringr str_extract
#' @importFrom stringr str_replace
#' @export

blast_malavi <- function(sequence, evalue = 1e-80, hits = 5, print.alignments = FALSE){

  ## set up html form
  base.url <- "http://130.235.244.92/Malavi/blast.html"
  form <- html_session(base.url) %>%
    html_form()

  ## set form values
  form <- set_values(form[[1]],"sequence" = sequence,
                     "evalue" = evalue,
                     "hits" = hits)

  ## submit form
  result <- submit_form(html_session(base.url), form)

  ## result
  x <- content(result$response, as = "text")

  if(!is.na(str_extract(x, "No hits found"))){
    warning("No hits found: check your input sequence")
    out.empty <- data.frame(Lineage = NA, Score = NA, Identities = NA, Gaps = NA, Strand = NA,
                            Coverage = NA, Perfect.Match = NA)
  } else{

    ## construct output table
    out.df <- data.frame(Lineage = str_extract(strsplit(x, "&gt")[[1]][2:(hits+1)], "[A-Z]+[0-9]+"),
                         Score = str_extract(strsplit(x, "&gt")[[1]][2:(hits+1)], "Score = [0-9]+") %>%
                           str_replace("Score = ", ""),
                         Identities = str_extract(strsplit(x, "&gt")[[1]][2:(hits+1)], "Identities = [0-9]+/[0-9]+") %>%
                           str_replace("Identities = ", ""),
                         Gaps = str_extract(strsplit(x, "&gt")[[1]][2:(hits+1)], "Gaps = [0-9]+/[0-9]+") %>%
                           str_replace("Gaps = ", ""),
                         Strand = str_extract(strsplit(x, "&gt")[[1]][2:(hits+1)], "Strand=[A-Z][a-z]+/[A-Z][a-z]+" %>%
                                                str_replace("Strand=", "")))

    out.df$Coverage <- str_extract(out.df$Identities, "/[0-9]+") %>%
      str_replace("/", "") %>%
      as.numeric

    out.df$Perfect.Match <- ifelse((str_extract(out.df$Identities, "[0-9]+/") %>%
                                      str_replace("/", "") %>%
                                      as.numeric)/out.df$Coverage == 1, "Yes", "No")
  }

  if(print.alignments == TRUE){
    cat(x)
  }

  if(exists("out.df")){
    return(out.df)
  } else{
    return(out.empty)
  }
}
