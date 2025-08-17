library(dplyr)
library(Seurat)
library(BiocParallel); register(SnowParam(4))
library(scRepertoire)
library(stringr)  # for str_extract

root <- "C:/Users/sunny/OneDrive/Documents/GSE200996_RAW"

# 1) Find files
h5_files <- list.files(root, pattern = "^GSM\\d+.*\\.h5$", full.names = TRUE)
vdj_gz   <- list.files(root, recursive = TRUE,
                       pattern = "^GSM\\d+_filtered_contig_annotations.*\\.csv\\.gz$",
                       full.names = TRUE)

# 2) Helpers (all vectorized)
get_patient   <- function(x) str_extract(x, "P\\d+")
get_timepoint <- function(x) str_extract(x, "(?i)(pre|post)-Tx")  # case-insensitive
get_compart   <- function(x) {
  xl <- tolower(x)
  ifelse(grepl("tumor", xl), "Tumor",
         ifelse(grepl("pbmc", xl), "PBMC", "Unknown"))
}

# 3) Build maps using the shared key patient_timepoint_compartment
bh5  <- basename(h5_files)
bvdj <- basename(vdj_gz)

h5_map <- tibble(
  path = h5_files,
  base = bh5,
  patient = get_patient(bh5),
  timepoint = get_timepoint(bh5),
  compartment = get_compart(bh5)
) |> filter(!is.na(patient), !is.na(timepoint), !is.na(compartment)) |>
  mutate(key = paste(patient, timepoint, compartment, sep = "_"))

# VDJ: drop pooled PBMC files like "P01-P02_..." and sorted subsets
vdj_map <- tibble(
  path = vdj_gz,
  base = bvdj,
  patient = get_patient(bvdj),
  timepoint = get_timepoint(bvdj),
  compartment = get_compart(bvdj),
  pooled = grepl("P\\d+-P\\d+", bvdj) | grepl("sorted_", bvdj, ignore.case = TRUE)
) |>
  filter(!pooled, !is.na(patient), !is.na(timepoint), !is.na(compartment)) |>
  mutate(key = paste(patient, timepoint, compartment, sep = "_"))

pairs <- inner_join(h5_map, vdj_map, by = "key", suffix = c(".h5", ".vdj"))
nrow(pairs); head(pairs |> select(key, base.h5, base.vdj), 10)

# 4) Load matched samples
infer_type <- get_compart  # already vectorized

load_one <- function(h5_path, sample_id){
  expr <- Read10X_h5(h5_path)
  so <- CreateSeuratObject(expr, project = sample_id, min.cells = 3, min.features = 200)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  subset(so, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt < 20)
}

seurat_list <- list(); tcr_list <- list()
for(i in seq_len(nrow(pairs))){
  key <- pairs$key[i]
  so <- load_one(pairs$path.h5[i], sample_id = key)
  so$Key  <- key
  so$Type <- infer_type(pairs$base.h5[i])  # single string here; also fine if vector
  seurat_list[[key]] <- so
  
  tcr <- read.csv(gzfile(pairs$path.vdj[i]), stringsAsFactors = FALSE)
  tcr$Key <- key
  tcr_list[[key]] <- tcr
}

table(sapply(seurat_list, function(x) unique(x$Type)))

# 5) Prefer Tumor samples; fallback to all if none labeled
tumor_ids <- names(seurat_list)[sapply(seurat_list, function(x) any(x$Type == "Tumor"))]
if(length(tumor_ids) == 0) tumor_ids <- names(seurat_list)

obj <- Reduce(function(x,y) merge(x,y), seurat_list[tumor_ids])

obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData()
obj <- RunPCA(obj) |> FindNeighbors(dims = 1:30) |> FindClusters(resolution = 0.6) |> RunUMAP(dims = 1:30)
DimPlot(obj, group.by = "seurat_clusters", label = TRUE)

# Peek at one tumor VDJ table you already loaded into tcr_list
ex <- tcr_list[[ tumor_ids[1] ]]
table(ex$chain, useNA = "ifany")        # should show TRA/TRB
table(ex$productivity, useNA = "ifany") # or `productive` depending on file
head(colnames(ex))

# Helper to make anything logical TRUE/FALSE
to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) {
    lx <- tolower(trimws(x))
    return(lx %in% c("true","t","1","yes","y"))
  }
  rep(TRUE, length(x))
}

filter_good_contigs <- function(df) {
  has <- function(nm) nm %in% names(df)
  get <- function(nm, default) if (has(nm)) df[[nm]] else default
  
  is_cell         <- to_logical(get("is_cell",         TRUE))
  high_confidence <- to_logical(get("high_confidence", TRUE))
  chain           <- if (has("chain")) df$chain else NA_character_
  cdr3            <- get("cdr3",    NA_character_)
  cdr3_nt         <- get("cdr3_nt", NA_character_)
  
  keep <- is_cell &
    high_confidence &
    !is.na(chain) & chain %in% c("TRA","TRB") &
    ( (!is.na(cdr3) & nzchar(cdr3)) | (!is.na(cdr3_nt) & nzchar(cdr3_nt)) )
  
  df[keep, , drop = FALSE]
}

tcr_good <- lapply(tcr_list[tumor_ids], filter_good_contigs)

keep_paired <- function(df) {
  has_ab <- tapply(df$chain, df$barcode, function(v) any(v=="TRA") & any(v=="TRB"))
  paired <- names(has_ab)[has_ab]
  df[df$barcode %in% paired, , drop = FALSE]
}
tcr_paired <- lapply(tcr_good, keep_paired)

comb_tcr <- combineTCR(
  input.data = tcr_paired,
  samples    = tumor_ids,
  ID         = rep("HNSCC", length(tumor_ids)),
  removeNA   = TRUE,
  filterNonproductive = FALSE,
  filterMulti = TRUE
)

# 1) Inspect columns in your combined TCR (one element is enough)
cols <- names(comb_tcr[[1]])
cols

# Use gene-level clonotypes and group by sample (both exist in your data)
obj <- combineExpression(
  input.data = comb_tcr,
  sc.data    = obj,
  cloneCall  = "gene",     # maps to CTgene
  chain      = "both",
  group.by   = "sample",   # column in comb_tcr[[i]]
  proportion = TRUE,
  filterNA   = TRUE
)

# Re-prefix each sample BEFORE merge so cell names are unique and match scRepertoire
seurat_list_pref <- seurat_list
for (k in names(seurat_list_pref)) {
  seurat_list_pref[[k]] <- RenameCells(seurat_list_pref[[k]], add.cell.id = k)
}

# Merge tumor samples (or all if you prefer)
tumor_ids <- names(seurat_list_pref)[sapply(seurat_list_pref, function(x) any(x$Type=="Tumor"))]
if (length(tumor_ids) == 0) tumor_ids <- names(seurat_list_pref)

obj <- Reduce(function(x,y) merge(x,y), seurat_list_pref[tumor_ids])

# Re-run the standard processing quickly
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData()
obj <- RunPCA(obj) |> FindNeighbors(dims = 1:30) |> FindClusters(resolution = 0.6) |> RunUMAP(dims = 1:30)



# 1) Remove the extra "_HNSCC_" from every TCR table's barcode
comb_tcr_fix <- lapply(comb_tcr, function(df) {
  df$barcode <- gsub("_HNSCC_", "_", df$barcode, fixed = TRUE)
  df
})

# 2) Try attaching again (old scRepertoire API)
obj <- combineExpression(
  input.data = comb_tcr_fix,
  sc.data    = obj,
  cloneCall  = "gene",      # you have CTgene/CTnt/CTaa; "gene" is fine
  chain      = "both",
  group.by   = "sample",
  proportion = TRUE,
  filterNA   = TRUE
)
head(colnames(obj), 5)
head(comb_tcr_fix[[1]]$barcode, 5)
unique(comb_tcr_fix[[1]]$sample)[1:5]


cells <- colnames(obj)
bcs   <- unique(unlist(lapply(comb_tcr_fix, `[[`, "barcode")))
match_rate <- length(intersect(cells, bcs)) / length(bcs)
match_rate

obj <- combineExpression(
  input.data = comb_tcr_fix,
  sc.data    = obj,
  cloneCall  = "gene",   # you have CTgene/CTnt/CTaa
  chain      = "both",
  group.by   = NULL,     # <- avoid using 'sample' since some are NA
  proportion = TRUE,
  filterNA   = TRUE
)

length(intersect(cells, bcs))
length(bcs)
length(cells)

p_homeo <- clonalHomeostasis(comb_tcr_fix)
p_div   <- clonalDiversity(comb_tcr_fix, cloneCall = "gene")

library(dplyr)
tcr_all <- do.call(rbind, comb_tcr_fix)

top_tbl <- tcr_all %>%
  count(CTgene, sort = TRUE, name = "n_cells") %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  slice_head(n = 10)
top_tbl

library(dplyr)

# Get all clonotype data in one table
tcr_all <- do.call(rbind, comb_tcr_fix)

# Identify top 10 by cell count
top10_clones <- tcr_all %>%
  count(CTgene, sort = TRUE, name = "n_cells") %>%
  slice_head(n = 10) %>%
  pull(CTgene)

# Add a column to Seurat object metadata: "Top10" vs "Other"
obj$TopClone <- ifelse(obj$CTgene %in% top10_clones, obj$CTgene, "Other")

# Plot on UMAP
DimPlot(obj, group.by = "TopClone", label = FALSE) +
  ggtitle("Top 10 Clonotypes Highlighted")


cells <- trimws(colnames(obj))
bcs   <- trimws(unique(unlist(lapply(comb_tcr_fix, `[[`, "barcode"))))

length(cells); length(bcs)
sum(bcs %in% cells)                           # exact matches
sum(sub("-\\d+$","", bcs) %in% sub("-\\d+$","", cells))   # ignore trailing -1/-2
head(cells, 5); head(bcs, 5)




library(dplyr)

## 1) Rebuild a clean TCR table (one row per cell)
tcr_all <- bind_rows(comb_tcr_fix) %>%
  mutate(hasTRB = grepl("TRB", TCR1) | grepl("TRB", TCR2)) %>%
  arrange(desc(hasTRB)) %>%
  mutate(barcode_core = sub("-\\d+$","", barcode)) %>%   # strip -1/-2
  distinct(barcode_core, .keep_all = TRUE) %>%
  mutate(CTgene = as.character(CTgene))

nrow(tcr_all)  # ~ number of TCR-tagged cells

## 2) Build a mapping for ALL Seurat cells
cells_tbl <- tibble(
  cell       = colnames(obj),
  barcode_core = sub("-\\d+$","", colnames(obj))
)

## 3) Left-join to bring CTgene onto Seurat cells
map_tbl <- cells_tbl %>%
  left_join(tcr_all %>% select(barcode_core, CTgene), by = "barcode_core")

sum(!is.na(map_tbl$CTgene))         # <-- this MUST be >> 0
head(map_tbl, 3)

## 4) Write CTgene into Seurat metadata (overwrite if exists)
obj$CTgene <- NA_character_
obj$CTgene[match(map_tbl$cell, colnames(obj))] <- map_tbl$CTgene

## 5) Sanity checks
sum(!is.na(obj$CTgene))             # number of cells with a clonotype
table(is.na(obj$CTgene))[c("FALSE","TRUE")]
head(obj$CTgene, 10)


library(dplyr)
# Top 10 clonotypes (by cell count)
clone_sizes <- tcr_all %>% count(CTgene, name = "n_cells")
top10 <- clone_sizes %>% arrange(desc(n_cells)) %>% slice_head(n = 10) %>% pull(CTgene)

# Tag cells and plot
obj$TopClone <- ifelse(obj$CTgene %in% top10, obj$CTgene, "Other")
DimPlot(obj, group.by = "TopClone", label = FALSE) + ggtitle("Top 10 Clonotypes Highlighted")


table_high <- table(obj$CTgene)
length(top10)                                   # should be 10
sum(obj$CTgene %in% top10)                      # cells in top10 clones
sort(table_high[top10], decreasing = TRUE)      # sizes of those clones

topN <- 50
topN_clones <- names(sort(table_high, decreasing = TRUE))[1:topN]
obj$TopClone <- ifelse(obj$CTgene %in% topN_clones, obj$CTgene, "Other")
DimPlot(obj, group.by = "TopClone", label = FALSE) + ggtitle(paste0("Top ", topN, " Clonotypes"))

# we already computed clone_sizes and obj$cloneSize_n earlier
min_cells <- 10   # adjust (e.g., 10, 20, 50)
big_clones <- clone_sizes$CTgene[clone_sizes$n_cells >= min_cells]
obj$ExpandedClone <- ifelse(obj$CTgene %in% big_clones, obj$CTgene, "Other")
DimPlot(obj, group.by = "ExpandedClone", label = FALSE) +
  ggtitle(paste0("Clonotypes with \u2265", min_cells, " cells"))



FeaturePlot(obj, features = "cloneSize_n") + ggtitle("Clone size per cell")

obj$expanded_flag <- ifelse(obj$cloneSize_n >= 10, "Expanded (≥10)", "Not expanded")
DimPlot(obj, group.by = "expanded_flag", label = FALSE) + ggtitle("Expanded vs Not")


# If you stored a per-cell sample key (e.g., in obj$Key or obj$GSM), use that:
sample_col <- "Key"   # change if your metadata uses a different name

# top 10 clonotypes within each sample
library(dplyr)
top_by_sample <- tcr_all %>%
  count(sample, CTgene, name = "n_cells") %>%
  group_by(sample) %>% slice_max(n_cells, n = 10, with_ties = FALSE) %>% ungroup()

obj$TopCloneBySample <- "Other"
for (s in unique(top_by_sample$sample)) {
  top_s <- top_by_sample$CTgene[top_by_sample$sample == s]
  idx_s <- which(obj[[sample_col]] == s & obj$CTgene %in% top_s)
  obj$TopCloneBySample[idx_s] <- obj$CTgene[idx_s]
}
DimPlot(obj, group.by = "TopCloneBySample", label = FALSE, split.by = sample_col)












library(Seurat)
library(dplyr)
library(ggplot2)

dir.create("outputs", showWarnings = FALSE)

# Rebuild a clean one-row-per-barcode TCR table
tcr_all <- bind_rows(comb_tcr_fix) %>%
  mutate(hasTRB = grepl("TRB", TCR1) | grepl("TRB", TCR2),
         barcode_core = sub("-\\d+$","", barcode)) %>%
  arrange(desc(hasTRB)) %>%
  distinct(barcode_core, .keep_all = TRUE) %>%
  mutate(CTgene = as.character(CTgene))

# Clone sizes (global)
clone_sizes <- tcr_all %>% count(CTgene, name = "n_cells")

# Add per-cell clone size & bucket (if not already present)
if (is.null(obj$cloneSize_n)) {
  obj$cloneSize_n <- clone_sizes$n_cells[ match(obj$CTgene, clone_sizes$CTgene) ]
  obj$cloneSize_n[is.na(obj$cloneSize_n)] <- 0
}
if (is.null(obj$cloneType)) {
  obj$cloneType <- cut(
    obj$cloneSize_n,
    breaks = c(-Inf, 1, 5, 20, 100, Inf),
    labels = c("Singleton","Small (2-5)","Medium (6-20)","Large (21-100)","Hyperexpanded (>100)")
  )
}

# Derive a per-cell sample label from the cell name prefix
# (example cell: "P18_pre-Tx_Tumor_AAACCTG...-1")
if (is.null(obj$Sample)) {
  obj$Sample <- sub("_[ACGT]+-\\d+$","", colnames(obj))
}

# Extract simple factors from Sample for easy faceting
obj$Time       <- ifelse(grepl("pre-Tx",  obj$Sample, ignore.case = TRUE),"pre","post")
obj$Compartment<- ifelse(grepl("PBMC",    obj$Sample, ignore.case = TRUE),"PBMC",
                         ifelse(grepl("Tumor", obj$Sample, ignore.case = TRUE),"Tumor","NA"))

# Homeostasis (proportion in each bucket)
homeo_tbl <- as.data.frame(prop.table(table(obj$cloneType)))
colnames(homeo_tbl) <- c("Bucket","Prop")
gg_homeo <- ggplot(homeo_tbl, aes(Bucket, Prop)) +
  geom_col() + coord_flip() + theme_bw() +
  ggtitle("Clonal homeostasis (overall)")
ggsave("outputs/homeostasis_overall.png", gg_homeo, width = 6, height = 4, dpi = 300)

# Shannon diversity of clonotypes (global)
shannon <- function(x){ p <- x/sum(x); -sum(p*log(p)) }
div_shannon <- shannon(clone_sizes$n_cells)
div_shannon


# Map CTgene to per-cell, then summarize by Sample
per_sample_sizes <- tcr_all %>%
  select(CTgene, barcode_core) %>%
  right_join(
    data.frame(cell = colnames(obj),
               barcode_core = sub("-\\d+$","", colnames(obj)),
               Sample = obj$Sample,
               stringsAsFactors = FALSE),
    by = "barcode_core") %>%
  filter(!is.na(CTgene)) %>%
  count(Sample, CTgene, name = "n_cells")

# Diversity per sample
div_by_sample <- per_sample_sizes %>%
  group_by(Sample) %>%
  summarize(shannon = shannon(n_cells), .groups = "drop")

write.csv(div_by_sample, "outputs/diversity_by_sample.csv", row.names = FALSE)

# If Sample isn't present, derive it from the prefix of each cell name
# Example cell: "P18_pre-Tx_Tumor_AAACCTGAGAGCCCAA-1"  ->  Sample = "P18_pre-Tx_Tumor"
if (!"Sample" %in% colnames(obj@meta.data)) {
  obj$Sample <- sub("_[ACGT]+-\\d+$", "", colnames(obj))
}

# (Optional) sanity check
head(obj$Sample); length(unique(obj$Sample))


library(dplyr)

per_sample_sizes <- tcr_all %>%
  select(CTgene, barcode_core) %>%
  right_join(
    data.frame(cell = colnames(obj),
               barcode_core = sub("-\\d+$","", colnames(obj)),
               Sample = obj$Sample,
               stringsAsFactors = FALSE),
    by = "barcode_core"
  ) %>%
  filter(!is.na(CTgene)) %>%
  count(Sample, CTgene, name = "n_cells")

# Preview & save
head(per_sample_sizes)
write.csv(per_sample_sizes, "outputs/per_sample_clonotype_sizes.csv", row.names = FALSE)

shannon <- function(x){ p <- x/sum(x); -sum(p*log(p)) }

div_by_sample <- per_sample_sizes |>
  dplyr::group_by(Sample) |>
  dplyr::summarise(shannon = shannon(n_cells), .groups = "drop")

print(div_by_sample)
write.csv(div_by_sample, "outputs/diversity_by_sample.csv", row.names = FALSE)


homeo_by_sample <- data.frame(Sample = obj$Sample, Bucket = obj$cloneType) |>
  dplyr::filter(!is.na(Bucket)) |>
  dplyr::count(Sample, Bucket) |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(Prop = n/sum(n)) |>
  dplyr::ungroup()

gg_homeo_s <- ggplot2::ggplot(homeo_by_sample, ggplot2::aes(Prop, Bucket, fill = Bucket)) +
  ggplot2::geom_col() + ggplot2::facet_wrap(~Sample, scales = "free_x") +
  ggplot2::theme_bw() + ggplot2::ggtitle("Clonal homeostasis by sample")

ggplot2::ggsave("outputs/homeostasis_by_sample.png", gg_homeo_s, width = 12, height = 8, dpi = 300)


top_by_sample <- per_sample_sizes |>
  dplyr::group_by(Sample) |>
  dplyr::slice_max(n_cells, n = 10, with_ties = FALSE) |>
  dplyr::ungroup()

print(head(top_by_sample, 20))
write.csv(top_by_sample, "outputs/top_clonotypes_by_sample.csv", row.names = FALSE)


thr <- 10  # expansion threshold (cells per clonotype)
obj$ExpandedFlag <- ifelse(obj$cloneSize_n >= thr, paste0("Expanded (≥",thr,")"), "Not expanded")

p_split <- DimPlot(obj, group.by = "ExpandedFlag", split.by = "Sample", ncol = 3) +
  ggtitle(paste0("Expanded vs Not by Sample (≥",thr,")"))

ggplot2::ggsave("outputs/umap_expanded_by_sample.png", p_split, width = 14, height = 10, dpi = 300)

saveRDS(obj, "outputs/hnscc_with_tcr_annotated.rds")


library(Seurat); library(dplyr); library(ggplot2)
dir.create("outputs_bio", showWarnings = FALSE)

# ---- marker panels (compact & periodontal-relevant) ----
marker_list <- list(
  Treg       = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  Th17       = c("RORC","CCR6","IL17A","IL17F","IL23R"),
  Th1        = c("TBX21","IFNG","CXCR3","STAT4"),
  Th2        = c("GATA3","IL4","IL13","CCR4"),
  Tfh        = c("CXCR5","BCL6","IL21","PDCD1"),
  CD8_cyt    = c("CD8A","CD8B","PRF1","GZMB","NKG7","KLRD1"),
  Exhaustion = c("PDCD1","LAG3","TIGIT","CTLA4","HAVCR2","TOX"),
  MAIT       = c("KLRB1","TRAV1-2","SLC4A10","ZBTB16","DPP4","IL7R"),
  GammaDelta = c("TRDC","TRGC1","TRGC2"),
  Bcell      = c("MS4A1","CD79A","CD79B","CD74","BANK1"),
  Plasma     = c("MZB1","XBP1","JCHAIN","SDC1") # SDC1 = CD138
)

# Add one module score per panel
scores <- AddModuleScore(obj, features = marker_list, name = "sig_", search = TRUE)
# Seurat returns columns sig_1 ... sig_n in obj@meta.data; map them back:
score_names <- paste0("sig_", seq_along(marker_list))
for (i in seq_along(marker_list)) {
  obj[[names(marker_list)[i]]] <- obj@meta.data[[score_names[i]]]
}

# Simple rules to call a single "Subset" per cell
subset_matrix <- FetchData(obj, vars = names(marker_list))
best <- colnames(subset_matrix)[max.col(subset_matrix, ties.method = "first")]
obj$Subset <- best

# helpful refinements
obj$Subset[which(obj$Bcell > obj$Plasma & obj$Bcell > 0.1)]  <- "Bcell"
obj$Subset[which(obj$Plasma > obj$Bcell & obj$Plasma > 0.1)] <- "Plasma"
# CD8 cytotoxic vs Exhausted: if both high, prefer Exhaustion
obj$Subset[which(obj$Exhaustion > obj$CD8_cyt & obj$Exhaustion > 0.1)] <- "Exhaustion"

# QC: visualize
p_subset <- DimPlot(obj, group.by = "Subset", label = TRUE) + ggtitle("Immune subsets (called)")
ggsave("outputs_bio/umap_subsets.png", p_subset, width = 8, height = 6, dpi = 300)






thr <- 10                                  # same threshold you used before
obj$ExpandedFlag <- ifelse(obj$cloneSize_n >= thr, "Expanded","NotExpanded")

tab_sub <- obj@meta.data %>%
  filter(!is.na(Subset)) %>%
  count(Subset, ExpandedFlag) %>%
  tidyr::pivot_wider(names_from = ExpandedFlag, values_from = n, values_fill = 0)

# Fisher’s exact per subset vs rest of cells
enrichment <- lapply(split(obj@meta.data, obj$Subset), function(df){
  a <- sum(df$cloneSize_n >= thr)                    # Expanded in this subset
  b <- nrow(df) - a                                  # Not expanded in this subset
  rest <- obj@meta.data[obj$Subset != unique(df$Subset), ]
  c <- sum(rest$cloneSize_n >= thr)                  # Expanded elsewhere
  d <- nrow(rest) - c
  ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
  data.frame(Subset = unique(df$Subset), OR = ft$estimate, p = ft$p.value,
             Expanded_in_subset = a, Cells_in_subset = nrow(df))
}) %>% bind_rows() %>%
  arrange(p)

write.csv(enrichment, "outputs_bio/enrichment_expanded_by_subset.csv", row.names = FALSE)

# bar showing % expanded by subset
prop_exp <- obj@meta.data %>%
  group_by(Subset) %>%
  summarize(pct_expanded = mean(cloneSize_n >= thr)*100, n = n()) %>%
  arrange(desc(pct_expanded))
g_prop <- ggplot(prop_exp, aes(reorder(Subset, pct_expanded), pct_expanded)) +
  geom_col() + coord_flip() + theme_bw() +
  ylab("% Expanded (≥10)") + xlab("") + ggtitle("Expanded clonotypes by subset")
ggsave("outputs_bio/expanded_pct_by_subset.png", g_prop, width = 6, height = 5, dpi = 300)







# Ensure these exist as before:
# obj$Sample like "P18_pre-Tx_Tumor"
if (!"Patient" %in% colnames(obj@meta.data)) obj$Patient <- sub("_.*","", obj$Sample)
if (!"Time"    %in% colnames(obj@meta.data)) obj$Time    <- ifelse(grepl("pre", obj$Sample, TRUE),"pre","post")

# overall per patient
patient_prop <- obj@meta.data %>%
  group_by(Patient, Time) %>%
  summarize(pct_expanded = mean(cloneSize_n >= thr)*100, .groups="drop")

# paired test across patients with both timepoints
paired_pts <- intersect(patient_prop$Patient[patient_prop$Time=="pre"],
                        patient_prop$Patient[patient_prop$Time=="post"])
wide <- tidyr::pivot_wider(filter(patient_prop, Patient %in% paired_pts),
                           names_from = Time, values_from = pct_expanded)
wilcox_res <- wilcox.test(wide$pre, wide$post, paired = TRUE)
wilcox_res

write.csv(patient_prop, "outputs_bio/patient_pct_expanded_pre_vs_post.csv", row.names = FALSE)

# per subset per patient (which biology changes?)
by_subset <- obj@meta.data %>%
  group_by(Patient, Time, Subset) %>%
  summarize(pct_expanded = mean(cloneSize_n >= thr)*100, .groups="drop")
write.csv(by_subset, "outputs_bio/patient_subset_pct_expanded.csv", row.names = FALSE)

# plot a few key subsets across pre/post (e.g., CD8_cyt, Th17, Treg)
focus <- c("CD8_cyt","Exhaustion","Th17","Treg")
g_sub <- by_subset %>%
  filter(Subset %in% focus, Patient %in% paired_pts) %>%
  ggplot(aes(Time, pct_expanded, group = Patient)) +
  geom_line(alpha=.4) + geom_point() +
  facet_wrap(~Subset, scales="free_y") + theme_bw() +
  ylab("% Expanded (≥10)") + ggtitle("Pre vs Post by subset (paired by patient)")
ggsave("outputs_bio/pre_post_expanded_by_subset.png", g_sub, width = 9, height = 6, dpi = 300)





# RANKL & OPG
genes_check <- c("TNFSF11","TNFRSF11B")  # RANKL, OPG
present <- genes_check[genes_check %in% rownames(obj)]
g_rankl_umap <- FeaturePlot(obj, features = present) + ggtitle("RANKL/OPG")
ggsave("outputs_bio/feature_RANKL_OPG.png", g_rankl_umap, width = 8, height = 5, dpi = 300)

# Which subsets express RANKL? Are expanded clones higher?
df_rankl <- FetchData(obj, vars = c("TNFSF11","Subset","ExpandedFlag"))
df_rankl$Subset_simple <- ifelse(df_rankl$Subset %in% c("Bcell","Plasma"), df_rankl$Subset,
                                 ifelse(grepl("^T", df_rankl$Subset) | df_rankl$Subset %in% c("CD8_cyt","Exhaustion","MAIT","GammaDelta","Th1","Th2","Th17","Tfh","Treg"),
                                        "Tcell","Other"))
g_vln <- ggplot(df_rankl %>% filter(!is.na(TNFSF11), Subset_simple %in% c("Tcell","Bcell","Plasma")),
                aes(ExpandedFlag, TNFSF11)) +
  geom_violin(fill="grey90") + geom_boxplot(width=.15, outlier.size=0.5) +
  facet_wrap(~Subset_simple, scales="free_y") + theme_bw() +
  ylab("TNFSF11 (RANKL) expr") + xlab("") +
  ggtitle("RANKL in expanded vs not (T/B lineages)")
ggsave("outputs_bio/rankl_expanded_vs_not.png", g_vln, width = 8, height = 5, dpi = 300)

# quick test within T cells
t_T <- with(subset(df_rankl, Subset_simple=="Tcell"),
            wilcox.test(TNFSF11 ~ ExpandedFlag))
t_B <- with(subset(df_rankl, Subset_simple=="Bcell"),
            wilcox.test(TNFSF11 ~ ExpandedFlag))
t_T; t_B



# color by subset + expansion
obj$Subset_Exp <- paste(obj$Subset, ifelse(obj$cloneSize_n >= thr, "Expanded","Other"), sep=":")
p_combo <- DimPlot(obj, group.by = "Subset_Exp", label = FALSE) +
  ggtitle("Expanded clones within functional subsets")
ggsave("outputs_bio/umap_subset_expanded.png", p_combo, width = 9, height = 6, dpi = 300)

# differential expression inside CD8 T cells: expanded vs not
cd8_cells <- WhichCells(obj, expression = CD8A > 0 | CD8B > 0)
obj_cd8 <- subset(obj, cells = cd8_cells)
Idents(obj_cd8) <- ifelse(obj_cd8$cloneSize_n >= thr, "Expanded","NotExpanded")
de_cd8 <- FindMarkers(obj_cd8, ident.1 = "Expanded", ident.2 = "NotExpanded",
                      test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1)
write.csv(de_cd8, "outputs_bio/DE_CD8_Expanded_vs_Not.csv")



# MAIT often TRAV1-2; γδ cells lack TRA/TRB & express TRDC
obj$MAIT_CT <- grepl("TRAV1-2", obj$CTgene %||% "")
obj$GammaDelta_CT <- grepl("TRD", obj$CTgene %||% "")
table(obj$MAIT_CT, obj$Subset)       # cross-check with MAIT score
table(obj$GammaDelta_CT, obj$Subset)





FeaturePlot(obj, features = c("CD4", "CD8A", "CCR6", "RORC", "IL17A"))



library(Seurat)
library(dplyr)
library(ggplot2)

# --- 0) set expansion flag if you haven't already ---
thr <- 10
if (!"ExpandedFlag" %in% colnames(obj@meta.data)) {
  obj$ExpandedFlag <- ifelse(obj$cloneSize_n >= thr, "Expanded","NotExpanded")
}

# --- 1) fetch expression for marker genes from the RNA assay ---
DefaultAssay(obj) <- "RNA"
genes <- c("CD4","CD8A","CD8B","CCR6","RORC","IL23R","IL17A")
present <- intersect(genes, rownames(obj))
stopifnot(length(present) > 0)

expr <- FetchData(obj, vars = present, slot = "data")  # log-normalized values

# --- 2) make simple scores ---
CD8_cols   <- intersect(c("CD8A","CD8B"), colnames(expr))
Th17_cols  <- intersect(c("RORC","CCR6","IL23R","IL17A"), colnames(expr))

obj$CD8_score  <- if (length(CD8_cols))  rowMeans(expr[, CD8_cols, drop=FALSE],  na.rm=TRUE) else 0
obj$CD4_score  <- if ("CD4" %in% colnames(expr)) expr[, "CD4"] else 0
obj$Th17_score <- if (length(Th17_cols)) rowMeans(expr[, Th17_cols, drop=FALSE], na.rm=TRUE) else 0

# data-driven thresholds (70th percentile works well across datasets)
qCD4  <- quantile(obj$CD4_score,  0.70, na.rm=TRUE)
qCD8  <- quantile(obj$CD8_score,  0.70, na.rm=TRUE)
qTh17 <- quantile(obj$Th17_score, 0.75, na.rm=TRUE)  # slightly stricter

# --- 3) assign subset label (CD8 beats CD4 if both high; Th17 is a CD4 sublineage) ---
obj$Subset_simple <- "Other"
obj$Subset_simple[obj$CD8_score >= qCD8 & obj$CD4_score < qCD4] <- "CD8_T"
obj$Subset_simple[obj$CD4_score >= qCD4]                        <- "CD4_T"
obj$Subset_simple[obj$CD4_score >= qCD4 & obj$Th17_score >= qTh17] <- "Th17"

# quick viz
p1 <- DimPlot(obj, group.by = "Subset_simple", label = TRUE) + ggtitle("CD4/CD8/Th17 calls")
ggsave("subset_calls_umap.png", p1, width=7, height=5, dpi=300)

# --- 4) quantify expansion by subset ---
sub_exp <- obj@meta.data %>%
  filter(Subset_simple %in% c("CD8_T","CD4_T","Th17")) %>%
  group_by(Subset_simple) %>%
  summarize(
    n_cells = n(),
    n_exp   = sum(ExpandedFlag == "Expanded"),
    pct_expanded = 100 * n_exp / n_cells
  ) %>% arrange(desc(pct_expanded))
sub_exp
write.csv(sub_exp, "subset_expansion_summary.csv", row.names = FALSE)

# stacked bar: proportion expanded per subset
sub_exp_long <- obj@meta.data %>%
  filter(Subset_simple %in% c("CD8_T","CD4_T","Th17")) %>%
  count(Subset_simple, ExpandedFlag) %>%
  group_by(Subset_simple) %>%
  mutate(prop = n/sum(n))

g_bar <- ggplot(sub_exp_long, aes(Subset_simple, prop, fill = ExpandedFlag)) +
  geom_col() + coord_flip() + theme_bw() +
  ylab("Proportion of cells") + xlab("") +
  ggtitle("Expanded vs Not by subset (≥10 cells/clone)")
ggsave("subset_expanded_bar.png", g_bar, width=6, height=4, dpi=300)

# --- 5) statistical enrichment: does a subset over-represent expanded clones? ---
enrich <- lapply(split(obj@meta.data, obj$Subset_simple), function(df){
  a <- sum(df$ExpandedFlag == "Expanded")
  b <- nrow(df) - a
  rest <- obj@meta.data[obj$Subset_simple != unique(df$Subset_simple), ]
  c <- sum(rest$ExpandedFlag == "Expanded")
  d <- nrow(rest) - c
  ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
  data.frame(Subset = unique(df$Subset_simple), OR = unname(ft$estimate), p = ft$p.value,
             Expanded_in_subset = a, Cells_in_subset = nrow(df))
}) %>% bind_rows() %>% arrange(p)
enrich
write.csv(enrich, "subset_expansion_enrichment.csv", row.names = FALSE)

# --- 6) optional: pre vs post within each subset (paired by patient) ---
if (!"Patient" %in% colnames(obj@meta.data)) obj$Patient <- sub("_.*","", obj$Sample)
if (!"Time"    %in% colnames(obj@meta.data)) obj$Time    <- ifelse(grepl("pre", obj$Sample, TRUE),"pre","post")

by_pt_sub <- obj@meta.data %>%
  filter(Subset_simple %in% c("CD8_T","CD4_T","Th17")) %>%
  group_by(Patient, Time, Subset_simple) %>%
  summarize(pct_expanded = 100 * mean(ExpandedFlag=="Expanded"), .groups="drop")

write.csv(by_pt_sub, "subset_pct_expanded_pre_post_by_patient.csv", row.names=FALSE)

# paired lines for a few patients having pre & post
paired <- by_pt_sub %>%
  group_by(Subset_simple) %>%
  filter(duplicated(Patient) | duplicated(Patient, fromLast = TRUE))
g_lines <- ggplot(paired, aes(Time, pct_expanded, group = Patient)) +
  geom_line(alpha=.5) + geom_point() +
  facet_wrap(~Subset_simple, scales = "free_y") + theme_bw() +
  ylab("% Expanded (≥10)") + ggtitle("Pre vs Post by subset (paired where available)")
ggsave("subset_pre_post_lines.png", g_lines, width=8, height=5, dpi=300)
