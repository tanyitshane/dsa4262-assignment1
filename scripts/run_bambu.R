library(bambu)
fa   <- "reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
gtf  <- "reference/Homo_sapiens.GRCh38.91.gtf"
bams <- list.files("bam", pattern = "\\.bam$", full.names = TRUE)

if (file.exists("bambu/annotations.rds")) {
  annotations <- readRDS("bambu/annotations.rds")
} else {
  annotations <- prepareAnnotations(gtf)
  saveRDS(annotations, "bambu/annotations.rds")
}

se <- bambu(reads = bams, annotations = annotations, genome = fa,
            ncore = 8, rcOutDir = "bambu/rc")
saveRDS(se, "bambu/se.rds")
writeBambuOutput(se, path = "bambu/")

sink("bambu/show_se.txt"); show(se); sink()
print(se)
