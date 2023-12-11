GOanalysis <- function(matS){
  background_genes <- keys(org.Hs.eg.db)
  signif_GeneID <- matS$GeneID
  params <- enrichGO(gene = signif_GeneID,
                     universe = background_genes,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     readable = TRUE)
  return(params)
}


