KEGGanalysis <- function(matS){
  signif_GeneID <- matS$GeneID
  KEGGparams <- enrichKEGG(gene = signif_GeneID,
                           organism = "hsa",
                           pvalueCutoff = 0.05) 
  return(KEGGparams)
}


