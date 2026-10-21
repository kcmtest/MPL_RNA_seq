# Differential Expression Workshop
# Synthetic RNA-seq counts: DESeq2 + edgeR
# No external data download required.

set.seed(2026)

required <- c(
  "DESeq2",
  "edgeR",
  "ggplot2",
  "dplyr",
  "tibble",
  "pheatmap"
)

missing_pkgs <- required[
  !vapply(
    required,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall Bioconductor packages DESeq2/edgeR and CRAN packages ",
    "ggplot2/dplyr/tibble/pheatmap before running this script."
  )
}

library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)

# -------------------------------------------------------------------
# 1. Metadata
# -------------------------------------------------------------------

metadata <- data.frame(
  sample = paste0("S", 1:8),
  condition = factor(
    rep(c("Control", "Treated"), each = 4),
    levels = c("Control", "Treated")
  ),
  batch = factor(
    rep(c("Batch1", "Batch2"), times = 4)
  )
)

rownames(metadata) <- metadata$sample

# -------------------------------------------------------------------
# 2. Synthetic counts
# -------------------------------------------------------------------

n_genes <- 1000
n_samples <- nrow(metadata)
gene_ids <- paste0("Gene", sprintf("%04d", 1:n_genes))

base_mean <- exp(
  rnorm(
    n_genes,
    mean = log(80),
    sd = 1
  )
)

up_genes <- 1:80
down_genes <- 81:160

true_log2fc <- rep(0, n_genes)
true_log2fc[up_genes] <- 1.5
true_log2fc[down_genes] <- -1.5

library_factor <- c(
  0.85, 1.05, 0.95, 1.15,
  0.90, 1.10, 1.00, 1.20
)

names(library_factor) <- metadata$sample

batch_log2fc <- 0.25

batch_factor <- ifelse(
  metadata$batch == "Batch2",
  2^batch_log2fc,
  1
)

dispersion <- runif(
  n_genes,
  min = 0.08,
  max = 0.35
)

count_matrix <- matrix(
  0L,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    gene_ids,
    metadata$sample
  )
)

for (j in seq_len(n_samples)) {

  condition_effect <- if (
    metadata$condition[j] == "Treated"
  ) {
    2^true_log2fc
  } else {
    rep(1, n_genes)
  }

  mu <- base_mean *
    condition_effect *
    library_factor[j] *
    batch_factor[j]

  count_matrix[, j] <- rnbinom(
    n_genes,
    mu = mu,
    size = 1 / dispersion
  )
}

stopifnot(
  identical(
    colnames(count_matrix),
    rownames(metadata)
  )
)

write.csv(
  count_matrix,
  "synthetic_gene_counts.csv"
)

write.csv(
  metadata,
  "synthetic_metadata.csv"
)

# -------------------------------------------------------------------
# 3. Filter
# -------------------------------------------------------------------

keep <- rowSums(
  count_matrix >= 10
) >= 3

counts_filtered <- count_matrix[
  keep,
  ,
  drop = FALSE
]

# -------------------------------------------------------------------
# 4. DESeq2
# -------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ batch + condition
)

dds <- DESeq(dds)

vsd <- vst(
  dds,
  blind = FALSE
)

plotPCA(
  vsd,
  intgroup = c(
    "condition",
    "batch"
  )
)

res_deseq <- results(
  dds,
  contrast = c(
    "condition",
    "Treated",
    "Control"
  ),
  alpha = 0.05
)

deseq_df <- as.data.frame(
  res_deseq
) |>
  rownames_to_column("gene_id") |>
  mutate(
    status = case_when(
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange >= 1 ~ "Up",

      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange <= -1 ~ "Down",

      TRUE ~ "Not significant"
    )
  )

print(table(deseq_df$status))

plotMA(
  res_deseq,
  alpha = 0.05,
  ylim = c(-4, 4)
)

volcano_df <- deseq_df |>
  mutate(
    minus_log10_padj =
      -log10(padj)
  )

print(
  ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = minus_log10_padj,
      shape = status
    )
  ) +
    geom_point(alpha = 0.65) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = 2
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2
    ) +
    theme_classic()
)

# -------------------------------------------------------------------
# 5. edgeR
# -------------------------------------------------------------------

dge <- DGEList(
  counts = counts_filtered,
  samples = metadata
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

design <- model.matrix(
  ~ batch + condition,
  data = metadata
)

dge <- estimateDisp(
  dge,
  design
)

fit <- glmQLFit(
  dge,
  design
)

qlf <- glmQLFTest(
  fit,
  coef = "conditionTreated"
)

edge_df <- topTags(
  qlf,
  n = Inf
)$table |>
  rownames_to_column("gene_id") |>
  mutate(
    status = case_when(
      FDR < 0.05 &
        logFC >= 1 ~ "Up",

      FDR < 0.05 &
        logFC <= -1 ~ "Down",

      TRUE ~ "Not significant"
    )
  )

print(table(edge_df$status))

# -------------------------------------------------------------------
# 6. Compare frameworks
# -------------------------------------------------------------------

comparison <- deseq_df |>
  select(
    gene_id,
    DESeq2_log2FC = log2FoldChange,
    DESeq2_padj = padj,
    DESeq2_status = status
  ) |>
  inner_join(
    edge_df |>
      select(
        gene_id,
        edgeR_log2FC = logFC,
        edgeR_FDR = FDR,
        edgeR_status = status
      ),
    by = "gene_id"
  )

print(
  cor(
    comparison$DESeq2_log2FC,
    comparison$edgeR_log2FC,
    use = "complete.obs"
  )
)

deseq_sig <- comparison |>
  filter(
    DESeq2_status != "Not significant"
  ) |>
  pull(gene_id)

edger_sig <- comparison |>
  filter(
    edgeR_status != "Not significant"
  ) |>
  pull(gene_id)

cat(
  "\nDESeq2 significant:", length(deseq_sig),
  "\nedgeR significant:", length(edger_sig),
  "\nOverlap:", length(intersect(deseq_sig, edger_sig)),
  "\n"
)

write.csv(
  deseq_df,
  "synthetic_DESeq2_results.csv",
  row.names = FALSE
)

write.csv(
  edge_df,
  "synthetic_edgeR_results.csv",
  row.names = FALSE
)

write.csv(
  comparison,
  "synthetic_DESeq2_edgeR_comparison.csv",
  row.names = FALSE
)
