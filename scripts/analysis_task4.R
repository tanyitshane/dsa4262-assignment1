suppressMessages(library(bambu))
se  <- readRDS("bambu/se.rds")
rd  <- rowData(se)
cnt <- assays(se)$counts
tot <- rowSums(cnt)
gn  <- read.delim("bambu/gene_names.tsv", header = FALSE, col.names = c("GENEID", "name"))
nm  <- function(ids) gn$name[match(ids, gn$GENEID)]
good <- !grepl("Hct116_directRNA_replicate6", colnames(se))   # failed sample (1% mapped)
cat("bambu version:", as.character(packageVersion("bambu")), "\n")

cat("\n== 4.2a novel transcripts\n"); print(table(novel = rd$novelTranscript))
cat("novel genes:", length(unique(rd$GENEID[rd$novelGene])), "\n")

cat("\n== 4.2b top 5 genes\n")
seg <- transcriptToGeneExpression(se)
g <- data.frame(GENEID = rownames(seg), name = nm(rownames(seg)),
                totalCounts = round(rowSums(assays(seg)$counts)),
                meanCPM_3good = round(rowMeans(assays(seg)$CPM[, good])))
cat("by total counts:\n"); print(head(g[order(-g$totalCounts), ], 5), row.names = FALSE)
cat("by mean CPM (excl. failed sample):\n"); print(head(g[order(-g$meanCPM_3good), ], 5), row.names = FALSE)

stats <- function(t) cat("genes:", length(t), " min:", min(t), " max:", max(t),
                         " mean:", round(mean(t), 2), "\n")
cat("\n== 4.2c transcripts per gene (all, extended annotations)\n"); stats(table(rd$GENEID))
cat("\n== 4.2d transcripts per gene (expressed: total counts >= 10)\n")
cat("expressed transcripts:", sum(tot >= 10), "\n"); stats(table(rd$GENEID[tot >= 10]))

nex <- lengths(rowRanges(se))
cat("\n== 4.2e single-exon novel transcripts:", sum(rd$novelTranscript & nex == 1), "\n")

cat("\n== 4.2f top transcripts with only unique reads\n")
onlyU <- rowSums(abs(cnt - assays(se)$uniqueCounts)) < 1e-6 & tot > 0
f <- data.frame(TXNAME = rd$TXNAME, GENEID = rd$GENEID, name = nm(rd$GENEID), total = round(tot, 1))[onlyU, ]
print(head(f[order(-f$total), ], 5), row.names = FALSE)

cat("\n== 4.4 candidates: most expressed novel transcripts\n")
nov <- which(rd$novelTranscript); r <- unlist(range(rowRanges(se)[nov]))
cls <- if ("txClassDescription" %in% colnames(rd)) rd$txClassDescription[nov] else NA
d <- data.frame(TXNAME = rd$TXNAME[nov], GENEID = rd$GENEID[nov], name = nm(rd$GENEID[nov]),
                class = cls, exons = nex[nov], total = round(tot[nov], 1),
                locus = paste0(seqnames(r), ":", start(r), "-", end(r), "(", strand(r), ")"))
print(head(d[order(-d$total), ], 10), row.names = FALSE)

cat("\n== warnings (for 5.3)\n"); print(metadata(se)$warnings)
