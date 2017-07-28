#' Download MalAvi alignment
#' @importFrom ape read.dna
#' @export

extract_alignment <- function(alignment = "all seqs"){
  alignment.names <- c("all seqs", "long seqs", "morpho seqs")
  alignment.urls <- c("http://mbio-serv2.mbioekol.lu.se/Malavi/PHP/get_all_sequences.php",
                      "http://mbio-serv2.mbioekol.lu.se/Malavi/PHP/get_long_sequences.php",
                      "http://mbio-serv2.mbioekol.lu.se/Malavi/PHP/get_MS_name_sequences_.php")
  names(alignment.urls) <- alignment.names
  if(!(alignment %in% alignment.names)){
    stop('Please choose one of the following alignment names: "all seqs", "long seqs", "morpho seqs"')
    return(c(cat("Please choose one of the following alignment names: "),
             cat(alignment.names, sep = ", ")))
  } else{
    return(read.dna(alignment.urls[alignment], format = "fasta"))
  }
}
