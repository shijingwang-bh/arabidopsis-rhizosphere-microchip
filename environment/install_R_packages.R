# Explicit optional installation; the reproduction runner never installs packages.
cran<-c('car','dplyr','jsonlite','vegan','statmod','ggplot2','patchwork','magick','ggrepel','circlize','BiocManager')
missing<-cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
bioc<-c('edgeR','limma','fgsea','BiocParallel')
missing<-bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))BiocManager::install(missing,ask=FALSE,update=FALSE)
cat('Installation uses the repositories available on your system. Consult the\nrecorded package-version tables for the environment actually verified.\n')
