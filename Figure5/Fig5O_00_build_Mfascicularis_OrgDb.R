suppressMessages(library(AnnotationForge))
suppressMessages(library(biomaRt))
setwd("monkey/orgdb_build")
message("== start (rebuildCache=FALSE) == ", Sys.time())
makeOrgPackageFromNCBI(
  version="0.1", author="admin", maintainer="admin <admin@example.org>",
  outputDir=".",
  tax_id="9541", genus="Macaca", species="fascicularis",
  NCBIFilesDir="ncbi",
  rebuildCache=FALSE)
message("== build done == ", Sys.time())
