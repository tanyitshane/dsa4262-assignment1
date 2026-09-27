suppressMessages(library(bambu))
se <- readRDS("bambu/se.rds")
genes <- c(RPS28P7 = "ENSG00000227097", OAZ1 = "ENSG00000104904", CHCHD6 = "ENSG00000159685")
for (n in names(genes)) {
  pdf(paste0("bambu/plot_", n, ".pdf"), width = 10, height = 6)
  p <- plotBambu(se, type = "annotation", gene_id = genes[[n]])
  if (!is.null(p)) print(p)
  dev.off()
}
