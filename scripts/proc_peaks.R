suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))

args <- commandArgs(TRUE)
mapping_path <- args[1]
if (!file.exists(mapping_path)) {
    stop(paste("invalid mapping path", mapping_path))
}
peaks_path <- args[2]
if (!file.exists(peaks_path)) {
    stop(paste("invalid peaks path", peaks_path))
}
p <- as.data.frame(read.table(peaks_path, header=F, sep="\t"))
p <- p[,c(1:3)]
colnames(p) <- c("chr","start","stop")
peaks <- suppressPackageStartupMessages(makeGRangesFromDataFrame(p))

data_name <- args[3]
method_name <- args[4]
ref_genome <- args[5]
out_path <- args[6]

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

annotations <- if(ref_genome=="hg38") GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) else GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- ref_genome

skip=0
if(method_name=="cell_ranger") {
  if (ref_genome == "hg38") {
    skip= 51
  } else {
    skip = 52
  }
  
}
print(paste("skip ", skip))
mapping_cells <- read_tsv(mapping_path, skip=skip,
    col_names=c("chr","start","stop","cell", "support"), 
    col_types=c("-","-","-","-","-"),
    col_select="cell") %>% 
    pull(cell) %>% 
    unique()

names(x = mapping_cells) <- paste(data_name, method_name, mapping_cells, sep="_")
mapping_frags <- CreateFragmentObject(path = mapping_path, cells = mapping_cells, max.lines=NULL)
mat <- FeatureMatrix(
  fragments = mapping_frags,
  features = peaks,
  process_n = 20000,
  sep = c("-", "-"),
  verbose = TRUE
) 
mapping_assay <- CreateChromatinAssay(mat, fragments = mapping_frags, 
    genome = ref_genome, min.features = 500)
seurat_ob <- CreateSeuratObject(mapping_assay, assay = "peaks")
seurat_ob$Sample <- paste(data_name, method_name, sep = "_")

fragmentInfo <- CountFragments(mapping_path)
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

seurat_ob_sub <- subset(x = seurat_ob,
        subset = nCount_peaks > 1000 &
        nCount_peaks < 100000 &
        FRiP > 0.15 &
        blacklist_fraction < 0.05 &
        nucleosome_signal < 4 &
        TSS.enrichment > 2
        )

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