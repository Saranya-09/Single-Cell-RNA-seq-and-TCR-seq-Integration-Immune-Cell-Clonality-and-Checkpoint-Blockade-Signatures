# Single-Cell-RNA-seq-and-TCR-seq-Integration-Immune-Cell-Clonality-and-Checkpoint-Blockade-Signatures
This project integrates single-cell RNA sequencing (scRNA-seq) and T-cell receptor sequencing (TCR-seq) data to explore how immune cells behave within the tumor microenvironment. By analyzing gene expression profiles together with TCR clonotypes, we uncover insights into immune cell diversity, clonal expansion, and checkpoint blockade mechanisms.

## 🔎 Overview (what this repo does)
This project integrates single-cell RNA-seq (scRNA-seq) with paired T-cell receptor sequencing (TCR-seq) from **GSE200996** (oral squamous cell carcinoma, neoadjuvant checkpoint blockade: anti-PD-1 ± anti-CTLA-4).  
I built an end-to-end **R** workflow to:
- load and QC scRNA-seq,
- ingest VDJ contig tables,
- **harmonize barcodes** and attach clonotypes,
- quantify **clonal expansion** (singleton→hyperexpanded),
- compare **pre- vs post-treatment** per patient,
- perform **CD8 vs Th17** subset expansion analysis, and
- map **RANKL (TNFSF11) / OPG** biology relevant to bone resorption.

---

## 📦 Dataset (what it is)
- **Accession**: GSE200996  
- **Design**: ~30 patients, **pre- and post-treatment** tumor + matched blood; two arms (anti-PD-1; anti-PD-1 + anti-CTLA-4).  
- **Modalities**:
  - **scRNA-seq**: 10x h5 (raw_feature_bc_matrix*.h5)
  - **TCR-seq (VDJ)**: filtered contig annotations (*.csv.gz)

> Note: Download from GEO and keep the raw folder structure. I used Windows paths like:
> `C:/Users/sunny/OneDrive/Documents/GSE200996_RAW/`

---

## 🧰 Tools & versions (R-based)
- **R 4.4.1** (Windows 11)
- **Seurat 5** (SeuratObject 5)
- **scRepertoire** (GitHub: ncborcherding/scRepertoire)
- **BiocParallel** (SnowParam on Windows)
- **dplyr, ggplot2**, and friends

**Local setup tips (Windows)**  
To avoid DLL permission issues:
```r


Run scripts in order

01_load_qc_seurat.R

Loads all *.h5, QC filters (nFeature_RNA > 300 & < 7000, percent.mt < 20), merges samples.

02_ingest_tcr_and_attach_clonotypes.R

Reads GSM*_filtered_contig_annotations*.csv.gz, fixes barcode mismatches
(removes _HNSCC_, strips -1/-2), attaches CTgene to cells.

✅ Results: ~42,952 cells with clonotype labels.

03_clonal_expansion_and_homeostasis.R

Computes cloneSize_n, defines expansion buckets (Singleton, Small, … Hyperexpanded),
plots UMAP by clone size and homeostasis (overall + per-sample).

04_pre_post_tx_comparison.R

Derives Patient and Time (pre/post) from sample keys,
computes % expanded (≥10) per patient × timepoint, plots paired lines.

05_subset_CD8_vs_Th17_and_RANKL.R

Fetches log-normalized expression and builds scores:
CD8 (CD8A/CD8B), CD4 (CD4), Th17 (RORC, CCR6, IL23R ± IL17A).

Assigns Subset_simple ∈ {CD8_T, CD4_T, Th17, Other} using data-driven quantile thresholds.

Crosses subsets with ExpandedFlag → % expanded by subset and Fisher enrichment.

Plots RANKL (TNFSF11) / OPG and compares RANKL in expanded vs not (T/B lineages).
.libPaths("C:/Rlibs")  # user-writable
library(BiocParallel); register(SnowParam(4))
