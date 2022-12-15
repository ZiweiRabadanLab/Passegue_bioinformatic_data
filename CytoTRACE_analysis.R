library(CytoTRACE)
# EA.combined <- readRDS('~/Passegue_bioinformatic_data/EA.combined.rds')
datasets <- as.data.frame(as.matrix(EA.combined@assays$RNA@counts))
batches_EA <- EA.combined$orig.ident
batches <- c(1:length(batches_EA))
for(i in 1:length(unique(batches_EA))){
  sampleid <- unique(batches_EA)[i]
  batches[which(batches_EA==sampleid)] <- i
}
results <- CytoTRACE(datasets,ncores = 1, batch = batches)
#saveRDS(results,'~/Passegue_bioinformatic_data/CytoTRACE_results.rds')

umap_EA <- readRDS('~/Passegue_bioinformatic_data/EA_monocle3_UMAP_Embedding.rds') # using the UMP embedding results from monocle3
phenotype_cluster <- as.character(Idents(EA.combined)) # using Seurat clustering results
names(phenotype_cluster) <- colnames(EA.combined)

plotCytoGenes(results, numOfGenes = 10, outputDir = '~/Passegue_bioinformatic_data/')
plotCytoTRACE(
  cyto_obj = results,
  phenotype = phenotype_cluster,
  emb = umap_EA,
  outputDir = '~/Passegue_bioinformatic_data/'
)
