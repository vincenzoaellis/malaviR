#' Download MalAvi data tables
#' @importFrom magrittr %>%
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr
#' @importFrom data.table fread
#' @export
extract_table <- function(table = "Hosts and Sites Table"){
  base.url <- "http://mbio-serv2.mbioekol.lu.se/"
  table.names <- c("Hosts and Sites Table", "Table of References", "Grand Lineage Summary",
                   "Parasite Summary Per Host", "Table of Lineage Names",
                   "Morpho Species Summary", "Vector Data Table", "Other Genes Table",
                   "Database Summary Report")
  table.urls <- c("http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report4=Hosts+And+Sites+Table",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report6=Table+of+References",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report1=Grand+Lineage+Summary",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report5=Parasite+Summary+Per+Host",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report7=Table+of+Lineage+Names",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report2=Morpho+Species+Summary",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report3=Vector+Data+Table",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report9=Other+Genes+Table",
                  "http://mbio-serv2.mbioekol.lu.se/bcgi/malaviReport.cgi?report8=Database+Summary+Report")
  names(table.urls) <- table.names
  if(!(table %in% c(table.names, "all"))){
    stop('Please choose one of the following table names:
         "Hosts and Sites Table", "Table of References", "Grand Lineage Summary",
         "Parasite Summary Per Host", "Table of Lineage Names", "Morpho Species Summary",
         "Vector Data Table", "Other Genes Table", "Database Summary Report", or "all"')
  }
  if(table == table.names[1]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[2]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[3]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[4]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[5]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[6]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[7]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[8]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == table.names[9]){
    data.url <- read_html(table.urls[table]) %>% html_nodes("a") %>% html_attr("href")
    data.url.n <- sub("../", "", data.url)
    return(fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE))
  }
  if(table == "all"){
    data.list <- list()
    for(i in 1:length(table.urls)){
      data.url <- read_html(table.urls[i]) %>% html_nodes("a") %>% html_attr("href")
      data.url.n <- sub("../", "", data.url)
      data.list[[i]] <- fread(paste(base.url, data.url.n, sep = ""), data.table = FALSE)
    }
    names(data.list) <- table.names
    return(data.list)
  }
}

