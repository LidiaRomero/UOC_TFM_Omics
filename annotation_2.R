annotation_2 <- function(input_case,matS, matA, variant_location_info){
  vcf <- readVcf(input_case, "hg19")
  matS$GeneID <- matA$GeneID[match(matS$Variant, matA$Variant)]
  matS$Protein_posi <- matA$Protein_posi[match(matS$Variant, matA$Variant)]
  matS$ref_AA <- matA$ref_AA[match(matS$Variant, matA$Variant)]
  matS$alt_AA <- matA$alt_AA[match(matS$Variant, matA$Variant)]
  matS$aaType <- matA$aaType[match(matS$Variant, matA$Variant)]
  matS$aaChange <- matA$aaChange[match(matS$Variant, matA$Variant)]
  matS$SampleID <- rownames(colData(vcf))[1]
  matS$VariantLoc <- variant_location_info$Location_Info[match(matS$Variant, variant_location_info$Variant_Ranges.names)]
  entrez_ids <- matS$GeneID
  gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,keys=entrez_ids, columns="SYMBOL", keytype="ENTREZID")
  matS$Gene <- gene_symbols$SYMBOL[match(matS$GeneID, gene_symbols$ENTREZID)]
  
  omim_symbol <- AnnotationDbi::select(org.Hs.eg.db,keys=entrez_ids, columns="OMIM", keytype="ENTREZID")
  matS$OMIM <- omim_symbol$OMIM[match(matS$GeneID, omim_symbol$ENTREZID)]
  
  GO_symbols <- AnnotationDbi::select(org.Hs.eg.db,keys=entrez_ids, columns="GO", keytype="ENTREZID")
  matS$GO <- GO_symbols$GO[match(matS$GeneID, omim_symbol$ENTREZID)]
  
  pp <- AnnotationDbi::select(PolyPhen.Hsapiens.dbSNP131, keys=matS$rsID,
                              cols=c("TRAININGSET", "PREDICTION", "PPH2PROB"))
  matS$polyphenpredict <- pp$PREDICTION[match(matS$rsID, pp$RSID)]
  
  ensembl_gene <- AnnotationDbi::select(org.Hs.eg.db,keys=entrez_ids, columns="ENSEMBL", keytype="ENTREZID")
  matS$ENSEMBL <- ensembl_gene$ENSEMBL[match(matS$GeneID, ensembl_gene$ENTREZID)]
  
  GO_ID <- matS$GO
  GO_term<- AnnotationDbi::select(GO.db,keys=GO_ID, columns="TERM", keytype="GOID")
  matS$GOterm <- GO_term$TERM[match(matS$GO, GO_term$GOID)]
  
  GO_onto <- GO_term<- AnnotationDbi::select(GO.db,keys=GO_ID, columns="ONTOLOGY", keytype="GOID")
  matS$GOOntology <- GO_onto$ONTOLOGY[match(matS$GO, GO_term$GOID)]
  
  genotypes <- geno(vcf)
  genotype_data <- data.frame(genotypes$GT)
  sample_name <- names(genotype_data[1])
  matS$genotype <- genotype_data[[sample_name]][match(matS$Variant, rownames(genotype_data))]
  allele_data <- data.frame(genotypes$AF)
  allelic_depth <- data.frame(genotypes$AD)
  allelic_depth <- allelic_depth[1]
  col_name <- paste("Reads", sample_name, sep = "")
  allelic_depth[[col_name]] <- format(allelic_depth[[sample_name]], big.mark = ",")
  split_allelic <- separate(allelic_depth, col_name, into = c("ReadsRef", "ReadsAlt"), sep = ",")
  matS$ReadsRef <- split_allelic$ReadsRef[match(matS$Variant, rownames(split_allelic))]
  matS$ReadsAlt <- split_allelic$ReadsAlt[match(matS$Variant, rownames(split_allelic))]
  
  gene_ids <- matS$ENSEMBL
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  omim_description <- getBM(
    attributes = c("ensembl_gene_id", "mim_morbid_description"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = ensembl
  )
  
  omim_description_combined <- omim_description %>% group_by(ensembl_gene_id) %>% summarize(mim_morbid_description = paste(mim_morbid_description, collapse = "/ "))
  matS$DiseaseOMIM <- omim_description_combined$mim_morbid_description[match(matS$ENSEMBL, omim_description_combined$ensembl_gene_id)]
  
  getType <- function(ref, alt){
    if(nchar(ref) == nchar(alt)){
      if(nchar(ref) == 1){
        return("SNP")
      } else {
        return("Complex variant")
      }
    } else if (nchar(ref) > nchar(alt)) {
      return("Deletion")
    } else {
      return("Insertion")
    }
  }
  matS$Type <- mapply(getType, matS$ref_allele, matS$alt_allele)
  
  freq_data <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
  minorallfreq <- getBM(attributes = c("refsnp_id", "minor_allele", "minor_allele_freq"), filters = "snp_filter", values = matS$rsID, mart = freq_data)
  matS$MAF <- minorallfreq$minor_allele_freq[match(matS$rsID,minorallfreq$refsnp_id)]
  return(matS)
}