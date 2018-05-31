# process adipose data
library(scran)
library(scater)
library(SingleCellExperiment)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

adipose <- read.table("~/Dropbox/Adipose/Adipose-counts.tsv",
                      h=T, sep="\t", stringsAsFactors=F)

adipose.meta <- read.table("~/Dropbox/Adipose/Adipose-metadata.tsv",
                           h=T, sep="\t", stringsAsFactors=F)

# get rid of the multi cell at this point
multi.cell <- adipose.meta$Sample[adipose.meta$Cell == "Multi"]
adipose <- adipose[, !colnames(adipose) %in% multi.cell]
adipose.meta <- adipose.meta[!adipose.meta$Sample%in% multi.cell, ]
# do some filtering on sparsity first
# no cells need remove, median cell sparsity is just 82%, very few cells 
# with > 90% sparsity
cell_sparsity <- apply(adipose == 0, 2, sum)/dim(adipose)[1]

gene_sparsity <- apply(adipose == 0, 1, sum)/dim(adipose)[2]

# remove genes with ?> 5% of missing values, equivalent to expression in < 50 cells
# is that reasonable?
hist(gene_sparsity, 100)
abline(vline=0.95, lty=3, col='purple', xlab="Gene Sparsity")
keep_genes <- gene_sparsity <= 0.95
adipose.nz <- adipose[keep_genes, ]

# bung it all into an SCE object and do some QC
# need better QC, try scater
spikes <- grepl(rownames(adipose.nz), pattern="ERCC")
genes <- rownames(adipose.nz)
adipose.nz <- apply(adipose.nz, 2, as.integer)
rownames(adipose.nz) <- genes
sce <- SingleCellExperiment(list("counts"=adipose.nz[, 1:(dim(adipose.nz)[2]-1)]),
                            colData=adipose.meta)
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))
isSpike(sce, "ERCC") <- spikes
