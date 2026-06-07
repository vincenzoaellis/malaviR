.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  packageStartupMessage(
    "malaviR ", v, " has undergone a major update. MalAvi is no longer downloaded\n",
    "from the web; the database is now bundled with the package as a dated snapshot.\n",
    "Please review the documentation (the \"Using malaviR\" vignette and the README)\n",
    "before using. Please report any problems at\n",
    "https://github.com/vincenzoaellis/malaviR/issues"
  )
}
