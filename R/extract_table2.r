#' Download MalAvi data tables
#' @importFrom magrittr %>%
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr
#' @importFrom data.table fread
#' @export
extract_table <- function(table = "Hosts and Sites Table"){
  base.url <- "http://130.235.244.92/"
  table.names <- c("Hosts and Sites Table", "Table of References", "Grand Lineage Summary",
                   "Parasite Summary Per Host", "Table of Lineage Names",
                   "Morpho Species Summary", "Vector Data Table", "Other Genes Table",
                   "Database Summary Report")
  table.urls <- c("http://130.235.244.92/bcgi/malaviReport.cgi?report4=Hosts+And+Sites+Table",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report6=Table+of+References",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report1=Grand+Lineage+Summary",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report5=Parasite+Summary+Per+Host",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report7=Table+of+Lineage+Names",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report2=Morpho+Species+Summary",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report3=Vector+Data+Table",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report9=Other+Genes+Table",
                  "http://130.235.244.92/bcgi/malaviReport.cgi?report8=Database+Summary+Report")
  names(table.urls) <- table.names
  if(!(table %in% c(table.names, "all"))){
    stop('Please choose one of the following table names:
         "Hosts and Sites Table", "Table of References", "Grand Lineage Summary",
         "Parasite Summary Per Host", "Table of Lineage Names", "Morpho Species Summary",
         "Vector Data Table", "Other Genes Table", "Database Summary Report", or "all"')
  }
  if(table == "all"){
    data.list <- list()
    for(i in 1:length(table.urls)){
      data.url <- read_html(table.urls[i]) %>%
        html_nodes("a") %>%
        html_attr("href")
      data.url.n <- sub("../", "", data.url)
      data.list[[i]] <- fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE)
    }
    names(data.list) <- table.names
    return(data.list)
  } else {
    data.url <- read_html(table.urls[table]) %>%
      html_nodes("a") %>%
      html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
}
