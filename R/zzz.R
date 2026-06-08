.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  packageStartupMessage(
    "malaviR ", v, " has undergone a major update. MalAvi is no longer downloaded\n",
    "from the web; the database is now bundled with the package.\n",
    "Please review the README and the function documentation before using.\n",
    "Please report any problems at\n",
    "https://github.com/vincenzoaellis/malaviR/issues"
  )
}
