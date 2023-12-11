annotation_1 <- function(input_case){
  vcf <- readVcf(input_case, "hg19")
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  all <- VariantAnnotation::locateVariants(vcf, txdb, AllVariants())
  variants_vcf_chr <- ranges(all)
  location_variants <- mcols(all)$LOCATION
  variant_location_info <- data.frame(Variant_Ranges = variants_vcf_chr,
                                      Location_Info = location_variants)
  variant_location_info <- unique(variant_location_info)
  coding <- predictCoding(vcf, txdb, seqSource = Hsapiens)
  matA <- data.frame(Variant=names(coding),
                     chromosome=seqnames(coding),
                     start=start(coding),end=end(coding),
                     ref_allele=as.character(coding$REF),
                     alt_allele=unlist(lapply(lapply(
                       coding$ALT,`[[`,1),as.character)),
                     GeneID=coding$GENEID,
                     Protein_posi=unlist(lapply(lapply(
                       coding$PROTEINLOC,`[[`,1),as.integer)),
                     ref_AA=as.character(coding$REFAA),
                     alt_AA=as.character(coding$VARAA),
                     aaType = coding$CONSEQUENCE)
  for(i in 1:length(matA$aaType)){
    if (matA$aaType[i] != "frameshift"){
      matA$aaChange[i]<- paste0("p.",matA$ref_AA[i],matA$Protein_posi[i],matA$alt_AA[i])
    }else{
      matA$aaChange[i] <- NA
    }
  }
  matA <- matA %>% distinct(Variant, .keep_all = TRUE)
  return(list(matA = matA, variant_location_info = variant_location_info))
}


