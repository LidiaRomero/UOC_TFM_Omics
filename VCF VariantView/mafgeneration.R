mafgeneration <- function(filtered_matS){
  matS_maf <- filtered_matS[,c ("Gene", "chromosome", "start", "end","ref_allele", "alt_allele", "Variant_Classification", "Type", "SampleID")]
  colnames(matS_maf) <- c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode')
  matS_maf$Start_Position <- as.integer(matS_maf$Start_Position)
  matS_maf$End_Position <- as.integer(matS_maf$End_Position)
  matS_maf$Variant_Classification <- as.character(matS_maf$Variant_Classification)
  matS_maf$AAChange <- filtered_matS$aaChange[match(matS_maf$Hugo_Symbol, filtered_matS$Gene)]
  maf_object <- read.maf(maf = matS_maf)
  return(maf_object)
}