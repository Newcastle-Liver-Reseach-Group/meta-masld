library(here)
library(edgeR)
library(limma)
library(tidyverse)



# function for simple PCA plots
mypca <- function(pcs, md, grp) {
  df <- as.data.frame(pcs$x)
  pcv <- round((pcs$sdev)^2 / sum(pcs$sdev^2) * 100, 2)
  df$grp <- pull(md, grp)
  g <- ggplot(df, aes(x=PC1, y=PC2)) +
    geom_point(aes(colour=grp)) + 
    theme_bw() + 
    xlab(label = paste0("PC1 (", pcv[1], "%)")) +
    ylab(label = paste0("PC2 (", pcv[2], "%)"))
  if (grp == 'fibrosis') {
    g <- g + scale_colour_manual(values=c(F0="#3B9AB2", 
                                          F1="#78B7C5",
                                          F2="#EBCC2A",
                                          F3="#E1AF00",
                                          F4="#F21A00"))
  }
  return(g)
}

# read the cohort metadata
# where is the data?
metafile <- here('metadata', 'meta_masld.tsv')
meta_masld <- read_tsv(metafile)

# filter out the IMID samples
meta_masld <- filter(meta_masld, masld_type == "classic-NAFLD")

# read the count data
# where are the counts?
countfile <- here('data', 'counts.csv.gz')
counts <- read.csv(countfile, sep=",", row.names=1) |>
  as.matrix

# the basic design
design  <- model.matrix(~ fibrosis + BioProject, data=meta_masld)

# Import the data
dge <- DGEList(counts, samples = meta_masld)
# filter by expression
keep <- filterByExpr(dge, design = design)
dge <- dge[keep,]
# normalise (TMM)
dge <- normLibSizes(dge)

# transform (voom) 
v <- voom(dge, design = design)

# calculate and plot PCA
pca <- prcomp(t(v$E))
g1 <- mypca(pca, meta_masld, 'BioProject')

# correct experiment effect for further viz
treatment.design <- design[,1:5]
batch.design <- design[,-(1:5)]
corrected_v <- removeBatchEffect(v, design=treatment.design, covariates=batch.design)

# calculate and plot experiment-corrected PCA
cpca <- prcomp(t(corrected_v))
g2 <- mypca(cpca, meta_masld, 'BioProject')
g3 <- mypca(cpca, meta_masld, 'fibrosis')

# analysis - "mild" (F0-2) vs "advanced" (F3-4), correcting for experiment

# where is the gene annotation?
genefile <- here('metadata', 'gencode50_annotation.csv')
gene_metadata <- read.csv(genefile, sep=",", row.names=1)
design <- model.matrix(~ 0 + fibrosis + BioProject, data = meta_masld)
v <- voom(dge, design)

fit <- lmFit(v, design)
contrast <- makeContrasts(test="((fibrosisF3 + fibrosisF4)/2) - ((fibrosisF0+fibrosisF1+fibrosisF2)/3)", 
                          levels = colnames(design))
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
fit2$genes <- gene_metadata[rownames(fit2),]
diffgenes <- topTable(fit2, number = Inf, p.value = 0.05)

# volcano plot of results
diffgenes$change <- ifelse(diffgenes$logFC > log2(1.5) & diffgenes$adj.P.Val < 0.05, 'up', 
                     ifelse(diffgenes$logFC < -log2(1.5) & diffgenes$adj.P.Val < 0.05, 'down', NA))
dge_labels <- c(diffgenes$external_gene_name[1:50], rep(NA, nrow(diffgenes)-50))

g4 <- ggplot(data=diffgenes, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=change)) + 
  scale_color_manual(values=c('#3A9AB2', '#F11B00'), na.value = 'lightgrey') +
  geom_hline(yintercept=-log10(0.05), linetype=2, colour='dodgerblue') +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), linetype=2, colour='dodgerblue') + 
  theme_minimal() +
  ggrepel::geom_label_repel(aes(label=dge_labels, color=change)) +
  guides(color='none')



