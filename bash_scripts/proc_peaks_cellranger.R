suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))

args <- commandArgs(TRUE)
cellranger_path <- args[1]
data_name <- args[2]
method_name <- args[3]
ref_genome <- args[4]
out_path <- args[5]

if (!file.exists(out_path)) {
  stop(paste("invalid path", out_path))
}
r_data_dir <- file.path(out_path, "rdata")
r_img_dir <- file.path(out_path, "images")
if(!file.exists(r_data_dir)) {
  dir.create(r_data_dir)
}
if(!file.exists(r_img_dir)) {
  dir.create(r_img_dir)
}

count_file <- file.path(cellranger_path, "outs", "filtered_peak_bc_matrix.h5")
if (!file.exists(count_file)) {
    stop(paste(count_file, " file does not exist"))
}
counts_mat <- Read10X_h5(count_file)
m_file <- file.path(cellranger_path, "outs", "singlecell.csv")
if (!file.exists(m_file)) {
    stop(paste(m_file, " file does not exist"))
}
meta_data <- read.csv(
  file = m_file,
  header = TRUE,
  row.names = 1
)

frag_file <- file.path(cellranger_path, "outs", "fragments.tsv.gz")
if (!file.exists(frag_file)) {
    stop(paste(frag_file, " file does not exist"))
}
mapping_assay <- CreateChromatinAssay(
  counts = counts_mat,
  sep = c(":", "-"),
  genome = ref_genome,
  fragments = frag_file,
  min.features = 500
)

seurat_ob <- CreateSeuratObject(mapping_assay, 
    assay = "peaks",
    meta.data = meta_data)
seurat_ob$Sample <- paste(data_name, method_name, sep = "_")

annotations <- if(ref_genome=="hg38") GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) else GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- ref_genome

Annotation(seurat_ob) <- annotations

fragmentInfo <- CountFragments(frag_file)
rownames(fragmentInfo) <- paste0(seurat_ob$Sample, "_", fragmentInfo$CB)

# Attach cell metadata to seurat object
seurat_ob$fragments <- fragmentInfo[colnames(seurat_ob), "frequency_count"]
seurat_ob$mononucleosomal <- fragmentInfo[colnames(seurat_ob), "mononucleosomal"]
seurat_ob$nucleosome_free <- fragmentInfo[colnames(seurat_ob), "nucleosome_free"]
seurat_ob$reads_count <- fragmentInfo[colnames(seurat_ob), "reads_count"]

# Calculate FRiP
seurat_ob <- FRiP(
  object = seurat_ob,
  assay = 'peaks',
  total.fragments = "fragments"
)

seurat_ob$blacklist_fraction <- if(ref_genome=="hg38") FractionCountsInRegion(
  object = seurat_ob, 
  assay = 'peaks',
  regions = blacklist_hg38
) else FractionCountsInRegion(
  object = seurat_ob, 
  assay = 'peaks',
  regions = blacklist_mm10
) 

# Compute nucleosome signal score per cell
seurat_ob <- NucleosomeSignal(seurat_ob)

# Compute TSS enrichment
Annotation(seurat_ob) <- annotations
seurat_ob <- TSSEnrichment(seurat_ob, fast=TRUE)

seurat_ob_sub <- RunTFIDF(seurat_ob)
seurat_ob_sub <- FindTopFeatures(seurat_ob_sub, min.cutoff = 'q0')
seurat_ob_sub <- RunSVD(seurat_ob_sub)

seurat_ob_sub <- RunUMAP(object = seurat_ob_sub, reduction = 'lsi', dims = 2:30)
seurat_ob_sub <- FindNeighbors(object = seurat_ob_sub, reduction = 'lsi', dims = 2:30)
seurat_ob_sub <- FindClusters(object = seurat_ob_sub, verbose = FALSE, algorithm = 3)


save(seurat_ob_sub, file = file.path(r_data_dir, "seurat_ob_sub.rdata"))

### Images
pdf(file.path(r_img_dir, "peaks_enrichment.pdf"))
VlnPlot(
  object = seurat_ob,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf(file.path(r_img_dir, "umap.pdf"))
DimPlot(object = seurat_ob_sub, label = TRUE, label.size = 7)
dev.off()
