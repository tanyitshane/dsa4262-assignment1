#!/usr/bin/env Rscript
suppressMessages(library(bambu))
args <- commandArgs(trailingOnly = TRUE)
useAnno <- as.logical(args[1]); fa <- args[2]; gtf <- args[3]; bams <- args[-(1:3)]
if (useAnno) {
  se <- bambu(reads = bams, annotations = prepareAnnotations(gtf), genome = fa, ncore = 1)
} else {
  se <- bambu(reads = bams, annotations = NULL, genome = fa, NDR = 1,
              opt.discovery = list(min.readFractionByEqClass = 0.2), ncore = 1)
}
writeBambuOutput(se, path = ".")
saveRDS(se, "se.rds")
