REACTOMEanalysis <- function(matS){
  background_genes <- keys(org.Hs.eg.db)
  signif_GeneID <- matS$GeneID
  react <- enrichPathway(gene = signif_GeneID,
                         universe = background_genes,
                         pvalueCutoff = 0.05, 
                         readable =TRUE)
  return(react)
}


