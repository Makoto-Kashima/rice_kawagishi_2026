# =============================================================================
# OsMYB80 / TDR anther RNA-seq: DEG, HOT genes and Fig. 5
#
#   Rscript OsMYB80_TDR_RNAseq_analysis.R
#
# Input (GEO supplementary files, placed in input/ or downloaded automatically
# when GEO_ACCESSION is set):
#   OsMYB80_TDR_anther_Salmon_NumReads.txt.gz     Salmon NumReads, transcripts x samples
#   OsMYB80_TDR_anther_sample_metadata.txt.gz     one row per sample
# Bundled in data/:
#   HOT_genes_Fig5.csv                            36 HOT genes shown in Fig. 5
#   high_temperature_repressed_71genes.csv        71 genes (Endo et al. 2009), optional
#   GO/                                           GO annotation for enrichment analysis
# Output: results/
#
# Upstream processing (fastp, BWA-MEM, Salmon 1.4.0; RAP-DB IRGSP-1.0
# representative transcripts 2022-09-01) followed Kashima et al. (2021)
# Front. Bioinform. 1:777299; see README.md.
#
# Package versions used for the published results: R 4.1, Seurat 4.2.0,
# SeuratObject 4.1.2, TCC 1.34.0, edgeR 3.36.0, clusterProfiler 4.2.2
# (Bioconductor 3.14). Newer edgeR versions give different TCC results.
# =============================================================================

suppressMessages({
  library(Seurat)
  library(TCC)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(AnnotationHub)
})

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
INPUT_DIR     <- Sys.getenv("INPUT_DIR", unset = "input")
GEO_ACCESSION <- Sys.getenv("GEO_ACCESSION", unset = "")   # e.g. "GSE000000"
DATA_DIR      <- "data"
RESULTS_DIR   <- "results"

GENOTYPE_LEVELS <- c("WT", "myb80", "tdr")
MUTANTS         <- c("myb80", "tdr")
SCALE_FACTOR    <- 1e6    # NormalizeData(scale.factor = 10^6)
TCC_ITERATION   <- 3      # calcNormFactors(iteration = 3)
TCC_FDR         <- 0.05   # estimateDE(FDR = 0.05)
LOG2FC_CUTOFF   <- 1      # |log2FC| > 1: used for the Fig. 5 overlap (see README)
RUN_GO          <- TRUE

COUNT_FILE <- "OsMYB80_TDR_anther_Salmon_NumReads.txt.gz"
META_FILE  <- "OsMYB80_TDR_anther_sample_metadata.txt.gz"

dir.create(INPUT_DIR, showWarnings = FALSE)
for (d in c("", "intermediate", "gene_sets", "GO")) {
  dir.create(file.path(RESULTS_DIR, d), showWarnings = FALSE, recursive = TRUE)
}
out <- function(...) file.path(RESULTS_DIR, ...)
strip_suffix <- function(x) unique(sub("\\.[0-9]+$", "", x))   # feature -> RAP locus

# -----------------------------------------------------------------------------
# 1. Load GEO count matrix and sample metadata
# -----------------------------------------------------------------------------
for (f in c(COUNT_FILE, META_FILE)) {
  if (!file.exists(file.path(INPUT_DIR, f))) {
    if (!nzchar(GEO_ACCESSION)) stop("Missing ", file.path(INPUT_DIR, f), " (set GEO_ACCESSION to download)")
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/%s",
                   substr(GEO_ACCESSION, 1, nchar(GEO_ACCESSION) - 3), GEO_ACCESSION, f)
    download.file(url, file.path(INPUT_DIR, f), mode = "wb")
  }
}

meta <- read.delim(file.path(INPUT_DIR, META_FILE), stringsAsFactors = FALSE)
meta <- meta[order(meta$analysis_order), ]      # column order used in the original analysis
rownames(meta) <- meta$sample_title

tab <- read.delim(file.path(INPUT_DIR, COUNT_FILE), check.names = FALSE, stringsAsFactors = FALSE)
# Features: one per RAP representative transcript (transcripts are not summed per gene).
# Os04t0470600-01 -> Os04g0470600; further transcripts of the same locus get ".1", ".2", ...
# (make.unique, as done by Seurat's CreateSeuratObject in the original analysis).
features <- data.frame(feature_id    = tab$feature_id_in_analysis,
                       transcript_id = tab$transcript_id,
                       gene_id       = tab$RAP_locus, stringsAsFactors = FALSE)
stopifnot(identical(features$feature_id,
                    make.unique(gsub("t", "g", sub("-.*", "", features$transcript_id)))))
counts <- as.matrix(tab[, meta$sample_title])
rownames(counts) <- features$feature_id
message(sprintf("Counts: %d features (%d loci) x %d samples", nrow(counts),
                length(unique(features$gene_id)), ncol(counts)))
print(table(factor(meta$genotype, levels = GENOTYPE_LEVELS)))

# -----------------------------------------------------------------------------
# 2. Seurat: normalisation, variable features, scaling, PCA
# -----------------------------------------------------------------------------
seurat.myb80.tdr <- CreateSeuratObject(counts = counts, meta.data = meta)
seurat.myb80.tdr <- NormalizeData(seurat.myb80.tdr, normalization.method = "LogNormalize",
                                  scale.factor = SCALE_FACTOR)
seurat.myb80.tdr <- FindVariableFeatures(seurat.myb80.tdr)
seurat.myb80.tdr <- ScaleData(seurat.myb80.tdr)
seurat.myb80.tdr <- RunPCA(seurat.myb80.tdr, verbose = FALSE)
seurat.myb80.tdr$genotype <- factor(seurat.myb80.tdr$genotype, levels = GENOTYPE_LEVELS)
saveRDS(seurat.myb80.tdr, out("intermediate", "seurat_myb80_tdr.rds"))

pdf(out("QC_seurat_myb80_tdr.pdf"), width = 9, height = 7)
print(VlnPlot(seurat.myb80.tdr, features = c("nCount_RNA", "nFeature_RNA"), group.by = "genotype", pt.size = 1))
print(ElbowPlot(seurat.myb80.tdr, ndims = 50))
print(DimPlot(seurat.myb80.tdr, reduction = "pca", group.by = "genotype", pt.size = 2) |
      DimPlot(seurat.myb80.tdr, reduction = "pca", group.by = "genotype", pt.size = 2, dims = c(1, 3)))
dev.off()

# -----------------------------------------------------------------------------
# 3. Differential expression with TCC (myb80 vs WT, tdr vs WT, separately)
#
# TCC's m.value is the log-ratio relative to the group of the FIRST sample
# column, not the first factor level. Here the first myb80-vs-WT column is a
# myb80 sample (m = log2 WT/myb80); the first tdr-vs-WT column is WT
# (m = log2 tdr/WT). Direction is therefore determined explicitly from the
# mean TCC-normalised counts, and log2FC_mut_vs_WT is m.value oriented
# accordingly (checked for every feature).
# -----------------------------------------------------------------------------
counts_all <- as.matrix(GetAssayData(seurat.myb80.tdr, slot = "counts"))
genotype <- as.character(seurat.myb80.tdr$genotype)

run_tcc <- function(mutant) {
  keep  <- genotype %in% c("WT", mutant)
  group <- factor(genotype[keep], levels = c("WT", mutant))
  tcc <- new("TCC", counts_all[, keep], group)
  tcc <- calcNormFactors(tcc, iteration = TCC_ITERATION)
  tcc <- estimateDE(tcc, FDR = TCC_FDR)
  res <- getResult(tcc, sort = FALSE)
  stopifnot(identical(res$gene_id, rownames(counts_all)))

  norm <- getNormalizedData(tcc)
  mean_WT  <- rowMeans(norm[, group == "WT", drop = FALSE])
  mean_mut <- rowMeans(norm[, group == mutant, drop = FALSE])
  lfc <- unname(log2(mean_mut) - log2(mean_WT))
  ok  <- is.finite(lfc) & is.finite(res$m.value) & lfc != 0
  orientation <- sign(sum(sign(lfc[ok]) * sign(res$m.value[ok])))
  stopifnot(isTRUE(all.equal(orientation * res$m.value[ok], lfc[ok], tolerance = 1e-6)))
  message(sprintf("%s vs WT: first column = %s; log2FC(mutant/WT) = %+d * m.value",
                  mutant, as.character(group[1]), orientation))

  r <- data.frame(feature_id       = res$gene_id,
                  gene_id          = features$gene_id,
                  transcript_id    = features$transcript_id,
                  mean_norm_WT     = mean_WT,
                  mean_norm_mutant = mean_mut,
                  a.value          = res$a.value,
                  m.value_TCC_raw  = res$m.value,
                  log2FC_mut_vs_WT = orientation * res$m.value,
                  p.value          = res$p.value,
                  q.value          = res$q.value,
                  rank             = res$rank,
                  estimatedDEG     = res$estimatedDEG,
                  stringsAsFactors = FALSE)
  r$direction <- ifelse(r$estimatedDEG == 1, ifelse(r$mean_norm_mutant > r$mean_norm_WT, "up", "down"), "NS")
  r$abs_log2FC_gt1 <- abs(r$log2FC_mut_vs_WT) > LOG2FC_CUTOFF
  list(result = r, norm.factors = tcc$norm.factors)
}

deg <- lapply(setNames(MUTANTS, MUTANTS), run_tcc)
saveRDS(deg, out("intermediate", "TCC_results.rds"))
for (m in MUTANTS) write.csv(deg[[m]]$result, out(sprintf("DEG_%s.csv", m)), row.names = FALSE)

deg_summary <- do.call(rbind, lapply(MUTANTS, function(m) {
  r <- deg[[m]]$result
  do.call(rbind, lapply(c("FDR=0.05", "FDR=0.05 & |log2FC|>1"), function(crit) {
    sel <- r$estimatedDEG == 1 & (crit == "FDR=0.05" | r$abs_log2FC_gt1)
    do.call(rbind, lapply(c("up", "down"), function(dirn) {
      f <- r$feature_id[sel & r$direction == dirn]
      data.frame(comparison = paste0(m, "_vs_WT"), criterion = crit, direction = dirn,
                 n_features = length(f), n_genes = length(strip_suffix(f)))
    }))
  }))
}))
write.csv(deg_summary, out("DEG_summary.csv"), row.names = FALSE)
print(deg_summary)

# -----------------------------------------------------------------------------
# 4. HOT genes and Fig. 5
#
# Downregulated = TCC estimatedDEG == 1 (FDR = 0.05, i.e. q < 0.05) and log2FC < 0.
# Reported for two definitions, because the Fig. 5 / main-text overlap
# (1,992 genes) is obtained only with the additional |log2FC| > 1 cutoff:
#   FDR      : FDR = 0.05
#   FDR_FC1  : FDR = 0.05 and log2FC < -1
# and at two levels: features (transcripts; as in the original lists) and
# unique RAP loci.
# -----------------------------------------------------------------------------
hot <- read.csv(file.path(DATA_DIR, "HOT_genes_Fig5.csv"), stringsAsFactors = FALSE)
ht71_file <- file.path(DATA_DIR, "high_temperature_repressed_71genes.csv")
ht71 <- if (file.exists(ht71_file)) read.csv(ht71_file, stringsAsFactors = FALSE)$RAP.locus else NULL
if (is.null(ht71)) message("NOTE: ", ht71_file, " not found; Venn regions involving the 71 genes = NA")

down_set <- function(m, fc) with(deg[[m]]$result, feature_id[estimatedDEG == 1 & log2FC_mut_vs_WT < -fc])

venn <- NULL
sets <- list()
for (def in c("FDR", "FDR_FC1")) {
  fc <- if (def == "FDR") 0 else LOG2FC_CUTOFF
  for (lvl in c("feature", "gene")) {
    A <- down_set("myb80", fc); B <- down_set("tdr", fc)
    if (lvl == "gene") { A <- strip_suffix(A); B <- strip_suffix(B) }
    AB <- intersect(A, B)
    sets[[paste(def, lvl, sep = "_")]] <- list(myb80_down = A, tdr_down = B, both_down = AB)
    C <- if (lvl == "gene") ht71 else NULL     # the 71 genes are RAP locus IDs
    na_if_noC <- function(x) if (is.null(C)) NA_integer_ else x
    venn <- rbind(venn, data.frame(
      definition = def, level = lvl,
      region = c("myb80_down (total)", "tdr_down (total)", "both_down (total)",
                 "myb80 only", "tdr only", "myb80 & tdr only", "all three (HOT)",
                 "myb80 & HT only", "tdr & HT only", "HT only", "HT (total)"),
      count = c(length(A), length(B), length(AB),
                length(setdiff(setdiff(A, B), C)), length(setdiff(setdiff(B, A), C)), length(setdiff(AB, C)),
                na_if_noC(length(intersect(AB, C))),
                na_if_noC(length(setdiff(intersect(A, C), B))),
                na_if_noC(length(setdiff(intersect(B, C), A))),
                na_if_noC(length(setdiff(C, union(A, B)))),
                na_if_noC(length(C))),
      note = ifelse(is.null(C) & seq_len(11) %in% 4:6, "HT set not applied", "")))
  }
}
write.csv(venn, out("Fig5_venn_counts.csv"), row.names = FALSE)
print(venn[venn$level == "feature", 1:4])
for (nm in names(sets)) for (s in names(sets[[nm]])) {
  writeLines(sort(sets[[nm]][[s]]), out("gene_sets", sprintf("%s_%s.txt", s, nm)))
}

# status of the 36 HOT genes in Fig. 5
hot_status <- hot[, c("order", "RAP.locus", "Symbol", "cluster_22K_KMC")]
for (m in MUTANTS) {
  r <- deg[[m]]$result
  i <- match(hot_status$RAP.locus, r$feature_id)            # primary transcript feature
  hot_status[[paste0(m, "_log2FC")]]        <- r$log2FC_mut_vs_WT[i]
  hot_status[[paste0(m, "_q.value")]]       <- r$q.value[i]
  hot_status[[paste0(m, "_down_FDR")]]      <- r$estimatedDEG[i] == 1 & r$log2FC_mut_vs_WT[i] < 0
  hot_status[[paste0(m, "_down_FDR_FC1")]]  <- r$estimatedDEG[i] == 1 & r$log2FC_mut_vs_WT[i] < -LOG2FC_CUTOFF
  hot_status[[paste0(m, "_down_FDR_FC1_anyTranscript")]] <-
    hot_status$RAP.locus %in% strip_suffix(down_set(m, LOG2FC_CUTOFF))
}
hot_status$both_down_FDR_FC1 <- hot_status$myb80_down_FDR_FC1 & hot_status$tdr_down_FDR_FC1
hot_status$both_down_FDR_FC1_anyTranscript <-
  hot_status$myb80_down_FDR_FC1_anyTranscript & hot_status$tdr_down_FDR_FC1_anyTranscript
hot_status$in_HT71 <- if (is.null(ht71)) NA else hot_status$RAP.locus %in% ht71
write.csv(hot_status, out("HOT_genes_status.csv"), row.names = FALSE)
message(sprintf("Fig. 5 HOT genes meeting the definition: %d/%d (primary transcript), %d/%d (any transcript)",
                sum(hot_status$both_down_FDR_FC1), nrow(hot_status),
                sum(hot_status$both_down_FDR_FC1_anyTranscript), nrow(hot_status)))

if (!is.null(ht71)) {
  s <- sets$FDR_FC1_gene
  hot_recomputed <- sort(intersect(s$both_down, ht71))
  write.csv(data.frame(RAP.locus = hot_recomputed, in_Fig5_list = hot_recomputed %in% hot$RAP.locus),
            out("HOT_genes_recomputed.csv"), row.names = FALSE)
  vd <- VennDiagram::venn.diagram(list(osmyb80 = s$myb80_down, tdr = s$tdr_down, `HT-repressed` = ht71),
                                  filename = NULL, disable.logging = TRUE)
  pdf(out("Fig5_venn.pdf")); grid::grid.draw(vd); dev.off()
}

# heat map (Fig. 5): all samples, HOT genes scaled with ScaleData
g1 <- DoHeatmap(ScaleData(seurat.myb80.tdr, features = hot$RAP.locus, verbose = FALSE),
                group.by = "genotype", features = hot$RAP.locus)
pdf(out("Fig5_heatmap.pdf")); print(g1); dev.off()

# alternative present in the original code (not used for Fig. 5):
# samples with auricle distance in (-2.01, 2.01); ScaleData on default (variable) features
sub <- subset(seurat.myb80.tdr, subset = auricle_distance > -2.01 & auricle_distance < 2.01)
sub <- ScaleData(sub, verbose = FALSE)
n_dropped <- length(setdiff(hot$RAP.locus, rownames(GetAssayData(sub, slot = "scale.data"))))
g2 <- suppressWarnings(DoHeatmap(sub, group.by = "genotype", features = hot$RAP.locus))
pdf(out("Fig5_heatmap_auricle_subset_ALTERNATIVE.pdf"))
print(g2 + labs(caption = sprintf("auricle distance in (-2.01, 2.01): n = %d; HOT genes not scaled: %d",
                                  ncol(sub), n_dropped)))
dev.off()

# -----------------------------------------------------------------------------
# 5. GO enrichment
#   A  (original analysis): enrichGO with AnnotationHub OrgDb AH96211
#      (org.Oryza_sativa_(japonica_cultivar-group).eg.sqlite; Bioconductor 3.14),
#      RAP -> Entrez via data/GO/riceIDtable.csv, BP, default cutoffs.
#        common       = myb80 down & tdr down (FDR = 0.05)
#        tdr.specific = tdr down, not myb80 down
#   B1 (HOT genes): same as A
#   B2 (HOT genes): enricher with RAP-DB (2021-11-11) + Oryzabase GO-BP annotation
#      (data/GO/goid2gene_BP.txt, godb_BP.txt); universe = all quantified loci
# -----------------------------------------------------------------------------
if (RUN_GO) {
  IDtable <- read.csv(file.path(DATA_DIR, "GO", "riceIDtable.csv"))
  to_entrez <- function(ids) {
    eid <- IDtable[match(ids, IDtable$rapdb), "entrezgene"]
    as.character(eid[!is.na(eid)])
  }
  save_go <- function(res, name) {
    df <- as.data.frame(res)
    write.csv(df, out("GO", paste0(name, ".csv")), row.names = FALSE)
    if (nrow(df) > 0) {
      pdf(out("GO", paste0(name, "_dotplot.pdf")), width = 8, height = 6)
      print(dotplot(res, showCategory = 30, font.size = 8, title = name))
      dev.off()
    }
    data.frame(analysis = name, n_input = length(res@gene), n_significant_terms = nrow(df))
  }

  myb80_down <- down_set("myb80", 0); tdr_down <- down_set("tdr", 0)
  common       <- myb80_down[is.element(myb80_down, tdr_down)]
  tdr.specific <- tdr_down[!is.element(tdr_down, myb80_down)]

  rice <- AnnotationHub()[["AH96211"]]
  run_enrichGO <- function(ids) enrichGO(gene = to_entrez(ids), OrgDb = rice, ont = "BP",
                                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                                         minGSSize = 10, maxGSSize = 500)
  term2gene <- read.delim(file.path(DATA_DIR, "GO", "goid2gene_BP.txt"), quote = "", stringsAsFactors = FALSE)
  term2name <- read.delim(file.path(DATA_DIR, "GO", "godb_BP.txt"), quote = "", stringsAsFactors = FALSE)

  go_summary <- rbind(
    save_go(run_enrichGO(common),        "A_enrichGO_AH96211_common_down"),
    save_go(run_enrichGO(tdr.specific),  "A_enrichGO_AH96211_tdr_specific_down"),
    save_go(run_enrichGO(hot$RAP.locus), "B1_enrichGO_AH96211_HOT_genes"),
    save_go(enricher(gene = hot$RAP.locus, universe = unique(features$gene_id),
                     TERM2GENE = term2gene, TERM2NAME = term2name,
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                     minGSSize = 10, maxGSSize = 500),
            "B2_enricher_RAPDB_Oryzabase_HOT_genes"))
  write.csv(go_summary, out("GO", "GO_summary.csv"), row.names = FALSE)
  print(go_summary)
}

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Results in ", RESULTS_DIR, "/")
