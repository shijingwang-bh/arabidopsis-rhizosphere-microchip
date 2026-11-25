args<-commandArgs(TRUE);stopifnot(length(args)==1L)
pkg<-normalizePath(args[1],winslash='/')
packages<-c('car','dplyr','jsonlite','vegan','edgeR','limma','statmod','fgsea','BiocParallel',
            'ggplot2','patchwork','magick','ggrepel','circlize')
stopifnot(all(vapply(packages,requireNamespace,logical(1),quietly=TRUE)))
result<-data.frame(package=c('R',packages),version=c(as.character(getRversion()),
                 vapply(packages,function(p)as.character(packageVersion(p)),character(1))))
write.csv(result,file.path(pkg,'environment/verified_package_versions.csv'),row.names=FALSE)
cat('Recorded',nrow(result),'verified R/package versions\n')
