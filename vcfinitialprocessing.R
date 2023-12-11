vcfinitialprocessing <- function(input_case) {
  vcf <- readVcf(input_case, "hg19")
  rd_chr <- rowRanges(vcf)
  all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  tar_chr <- as.vector(seqnames(rd_chr)@values)
  tar_chr <- gsub("chr", "", tar_chr)
  tar_chr[grepl(tar_chr, pattern = "M")] <- "MT"
  tar_chr_keep <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")
  condition <- tar_chr %in% tar_chr_keep
  tar_chr <- tar_chr[condition]
  my_snps <- snpsBySeqname(all_snps, c(tar_chr))
  seqlevelsStyle(my_snps) <- "UCSC"
  genome(my_snps) <- "hg19"
  snp_ID <- data.frame(posIDX = paste0(seqnames(my_snps), ":", pos(my_snps)), rsID = my_snps$RefSNP_id, stringsAsFactors = FALSE)
  matV1 <- data.frame(Variant = names(rd_chr), stringsAsFactors = FALSE)
  matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)", "\\1", matV1$Variant)
  matV1$start <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
  matV1$end <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
  matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\3", matV1$Variant)
  matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\4", matV1$Variant)
  matV1$posIDX <- gsub("(.*)_(.*)", "\\1", matV1$Variant)
  matS <- merge(matV1, snp_ID, all.x = TRUE, by = "posIDX")
  matS <- dplyr::select(matS, -posIDX)
  return(matS)
}
