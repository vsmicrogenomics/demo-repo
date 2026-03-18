#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

## ==============================
## User-configurable parameters
## ==============================

seurat_rds_path <- "path/to/seurat_object.rds"  # <- set me
output_dir <- "edger_pseudobulk_out"            # <- set me

min_cells_per_pseudobulk <- 20
min_replicates_per_condition <- 2
assay_name <- "RNA"

## Fixed for this analysis (requirements)
sample_col <- "orig.ident"
condition_col <- "condition"
cell_type_col <- "cell_type_l1"
condition_levels <- c("Healthy", "IPF") # explicit levels so contrast is IPF - Healthy

## ==============================
## Helpers
## ==============================

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

assert_cols_present <- function(df, cols, df_name = "data.frame") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns in ", df_name, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

assert_nonempty <- function(x, msg) {
  if (length(x) == 0) stop(msg, call. = FALSE)
}

safe_filename <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

## ==============================
## Load Seurat object
## ==============================

cat("[1/6] Loading Seurat object.\n")
if (!file.exists(seurat_rds_path)) stop("Seurat object not found at: ", seurat_rds_path, call. = FALSE)
obj <- readRDS(seurat_rds_path)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.", call. = FALSE)

assays <- names(obj@assays)
if (!(assay_name %in% assays)) {
  stop(
    "Missing assay '", assay_name, "'. Available assays: ",
    paste(assays, collapse = ", "),
    call. = FALSE
  )
}

## ==============================
## Extract counts + align metadata by cell names
## ==============================

cat("[2/6] Extracting raw counts and validating metadata.\n")
counts <- GetAssayData(obj, assay = assay_name, slot = "counts")
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

assert_nonempty(rownames(counts), "Counts matrix has no gene names (rownames).")
assert_nonempty(colnames(counts), "Counts matrix has no cell names (colnames).")

meta <- obj[[]]
assert_cols_present(meta, c(sample_col, condition_col, cell_type_col), df_name = "Seurat metadata")

if (is.null(rownames(meta)) || anyNA(rownames(meta))) {
  stop("Seurat metadata is missing rownames (cell names); cannot align to counts.", call. = FALSE)
}

missing_cells <- setdiff(colnames(counts), rownames(meta))
if (length(missing_cells) > 0) {
  stop(
    "Metadata missing for ", length(missing_cells), " cells found in counts. Example cell: ",
    missing_cells[[1]],
    call. = FALSE
  )
}
meta <- meta[colnames(counts), , drop = FALSE]

meta[[sample_col]] <- as.character(meta[[sample_col]])
meta[[condition_col]] <- as.character(meta[[condition_col]])
meta[[cell_type_col]] <- as.character(meta[[cell_type_col]])

if (anyNA(meta[[sample_col]]) || anyNA(meta[[condition_col]]) || anyNA(meta[[cell_type_col]])) {
  stop("NA values found in required metadata columns (orig.ident/condition/cell_type_l1).", call. = FALSE)
}

present_conditions <- sort(unique(meta[[condition_col]]))
if (!all(condition_levels %in% present_conditions)) {
  stop(
    "Condition levels missing. Need: ", paste(condition_levels, collapse = ", "),
    ". Present: ", paste(present_conditions, collapse = ", "),
    call. = FALSE
  )
}

## ==============================
## Pseudo-bulk aggregation (sparse-aware)
## ==============================

cat("[3/6] Aggregating to pseudo-bulk by orig.ident + cell_type_l1.\n")
pb_id <- paste(meta[[sample_col]], meta[[cell_type_col]], sep = "__")

pb_tab <- table(pb_id)
keep_pb <- pb_tab >= min_cells_per_pseudobulk
if (!any(keep_pb)) {
  stop(
    "After filtering pseudo-bulk groups with < ", min_cells_per_pseudobulk,
    " cells, no groups remain.",
    call. = FALSE
  )
}

keep_cells <- keep_pb[pb_id]
keep_cells[is.na(keep_cells)] <- FALSE

counts_f <- counts[, keep_cells, drop = FALSE]
meta_f <- meta[keep_cells, , drop = FALSE]
pb_id_f <- pb_id[keep_cells]

cat("  - Cells before: ", ncol(counts), "\n", sep = "")
cat("  - Cells after : ", ncol(counts_f), "\n", sep = "")
cat("  - Pseudo-bulk samples after filter: ", sum(keep_pb), "\n", sep = "")

# Build pseudo-bulk sample metadata from filtered cell metadata (no pb_id string parsing assumptions)
pb_meta_raw <- data.frame(
  pb_id = pb_id_f,
  orig.ident = meta_f[[sample_col]],
  cell_type_l1 = meta_f[[cell_type_col]],
  condition = meta_f[[condition_col]],
  stringsAsFactors = FALSE
)

# Verify each pb_id maps to exactly one orig.ident, one cell_type, and one condition
uniq_len <- function(x) length(unique(x))
bad_samples <- names(tapply(pb_meta_raw$orig.ident, pb_meta_raw$pb_id, uniq_len))[
  tapply(pb_meta_raw$orig.ident, pb_meta_raw$pb_id, uniq_len) != 1
]
bad_types <- names(tapply(pb_meta_raw$cell_type_l1, pb_meta_raw$pb_id, uniq_len))[
  tapply(pb_meta_raw$cell_type_l1, pb_meta_raw$pb_id, uniq_len) != 1
]
bad_conds <- names(tapply(pb_meta_raw$condition, pb_meta_raw$pb_id, uniq_len))[
  tapply(pb_meta_raw$condition, pb_meta_raw$pb_id, uniq_len) != 1
]
if (length(bad_samples) + length(bad_types) + length(bad_conds) > 0) {
  ex <- c(bad_samples, bad_types, bad_conds)[[1]]
  stop(
    "Inconsistent metadata within pseudo-bulk group (pb_id). Example pb_id: ", ex,
    call. = FALSE
  )
}

pb_levels <- sort(unique(pb_meta_raw$pb_id))
pb_meta <- unique(pb_meta_raw[, c("pb_id", "orig.ident", "cell_type_l1", "condition")])
pb_meta <- pb_meta[match(pb_levels, pb_meta$pb_id), , drop = FALSE]
pb_meta$condition <- factor(pb_meta$condition, levels = condition_levels)

pb_ncells <- as.integer(table(pb_id_f)[pb_meta$pb_id])
pb_meta$n_cells <- pb_ncells

# Sparse membership matrix (cells x pseudo-bulk) and aggregation: (genes x cells) %*% (cells x groups)
pb_index <- match(pb_id_f, pb_levels)
if (anyNA(pb_index)) stop("Internal error: pb_id mapping produced NA indices.", call. = FALSE)

G <- sparseMatrix(
  i = seq_len(length(pb_id_f)),
  j = pb_index,
  x = 1,
  dims = c(length(pb_id_f), length(pb_levels)),
  dimnames = list(colnames(counts_f), pb_levels)
)
pb_counts <- counts_f %*% G # dgCMatrix genes x pseudo-bulk

## ==============================
## edgeR QL DE per cell type
## ==============================

cat("[4/6] Running edgeR quasi-likelihood DE per cell type.\n")
make_dir(output_dir)
out_per_ct <- file.path(output_dir, "per_cell_type")
make_dir(out_per_ct)

summary_rows <- list()
de_rows <- list()

cell_types <- sort(unique(pb_meta$cell_type_l1))
if (length(cell_types) == 0) stop("No cell types found after pseudo-bulk filtering.", call. = FALSE)

for (ct in cell_types) {
  ct_mask <- pb_meta$cell_type_l1 == ct
  ct_meta <- pb_meta[ct_mask, , drop = FALSE]
  if (nrow(ct_meta) == 0) next

  # Replicate counts are pseudo-bulk samples (orig.ident) per condition
  repl_tab <- table(ct_meta$condition)
  healthy_n <- if ("Healthy" %in% names(repl_tab)) as.integer(repl_tab[["Healthy"]]) else 0L
  ipf_n <- if ("IPF" %in% names(repl_tab)) as.integer(repl_tab[["IPF"]]) else 0L

  summary_rows[[ct]] <- data.frame(
    cell_type_l1 = ct,
    n_cells = sum(ct_meta$n_cells, na.rm = TRUE),
    n_pseudobulk = nrow(ct_meta),
    n_repl_Healthy = healthy_n,
    n_repl_IPF = ipf_n,
    skipped = FALSE,
    skip_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  if (healthy_n < min_replicates_per_condition || ipf_n < min_replicates_per_condition) {
    summary_rows[[ct]]$skipped <- TRUE
    summary_rows[[ct]]$skip_reason <- paste0(
      "Insufficient replicates (Healthy=", healthy_n,
      ", IPF=", ipf_n, "; need >= ", min_replicates_per_condition, " each)"
    )
    cat("  - Skipping '", ct, "': insufficient replicates (Healthy=", healthy_n, ", IPF=", ipf_n, ").\n", sep = "")
    next
  }

  y_counts_sparse <- pb_counts[, ct_mask, drop = FALSE]
  if (ncol(y_counts_sparse) == 0) {
    summary_rows[[ct]]$skipped <- TRUE
    summary_rows[[ct]]$skip_reason <- "No pseudo-bulk samples"
    cat("  - Skipping '", ct, "': no pseudo-bulk samples.\n", sep = "")
    next
  }

  # Convert only this small pseudo-bulk slice to dense; avoid densifying full scRNA matrix
  y_counts <- as.matrix(y_counts_sparse)

  y <- DGEList(counts = y_counts)
  ct_meta$condition <- factor(ct_meta$condition, levels = condition_levels)

  design <- model.matrix(~ 0 + condition, data = ct_meta)
  if (!all(c("conditionHealthy", "conditionIPF") %in% colnames(design))) {
    stop(
      "Design matrix columns unexpected for cell type '", ct, "'. Got: ",
      paste(colnames(design), collapse = ", "),
      call. = FALSE
    )
  }

  keep_genes <- filterByExpr(y, design = design)
  if (!any(keep_genes)) {
    summary_rows[[ct]]$skipped <- TRUE
    summary_rows[[ct]]$skip_reason <- "No genes passed filterByExpr"
    cat("  - Skipping '", ct, "': no genes passed filterByExpr.\n", sep = "")
    next
  }

  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  con <- makeContrasts(IPF_vs_Healthy = conditionIPF - conditionHealthy, levels = design)
  qlf <- glmQLFTest(fit, contrast = con[, "IPF_vs_Healthy"])
  tab <- topTags(qlf, n = Inf)$table

  tab$gene <- rownames(tab)
  tab$cell_type_l1 <- ct
  tab$n_pseudobulk <- ncol(y$counts)
  tab$n_repl_Healthy <- healthy_n
  tab$n_repl_IPF <- ipf_n

  out_ct_file <- file.path(out_per_ct, paste0(safe_filename(ct), "_edgeR_QLF_IPF_vs_Healthy.csv"))
  write.csv(tab, out_ct_file, row.names = FALSE, quote = TRUE)
  cat("  - Wrote DE for '", ct, "' (genes=", nrow(tab), ") -> ", out_ct_file, "\n", sep = "")

  de_rows[[ct]] <- tab
}

## ==============================
## Write outputs
## ==============================

cat("[5/6] Writing outputs.\n")
summary_df <- do.call(rbind, summary_rows)
if (is.null(summary_df) || nrow(summary_df) == 0) stop("No summary rows produced; check inputs/filters.", call. = FALSE)

summary_file <- file.path(output_dir, "summary_pseudobulk_and_replicates.csv")
write.csv(summary_df, summary_file, row.names = FALSE, quote = TRUE)
cat("  - Summary table -> ", summary_file, "\n", sep = "")

if (length(de_rows) == 0) {
  cat("[6/6] Done. All cell types were skipped; see summary for reasons.\n")
  quit(save = "no", status = 0)
}

de_all <- do.call(rbind, de_rows)
combined_file <- file.path(output_dir, "edgeR_pseudobulk_DE_all_cell_types_IPF_vs_Healthy.csv")
write.csv(de_all, combined_file, row.names = FALSE, quote = TRUE)
cat("  - Combined DE results -> ", combined_file, "\n", sep = "")
cat("[6/6] Done.\n")