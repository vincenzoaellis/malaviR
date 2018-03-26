#' Get version number and upload date from MalAvi website
#' @importFrom magrittr %>%
#' @importFrom xml2 read_html
#' @importFrom rvest html_text
#' @importFrom stringr str_extract
#' @export

malavi_version <- function(){
  read_html("http://mbio-serv2.mbioekol.lu.se/Malavi/JavaScript/footer.js") %>%
    html_text() %>%
    str_extract(pattern = "Version [0-9].[0-9].[0-9], [A-Z][a-z]+ [0-9], [0-9]+")
}
