annotation_3 <- function(matS){
  matS$Type <- ifelse(matS$Type == "Deletion", "DEL",
                      ifelse(matS$Type == "Insertion", "INS", matS$Type))
  matS$chromosome <- sub("^chr", "", matS$chromosome)
  get_mutation_type <- function(ref, alt, var_type, var_class) {
    if (is.na(var_class)) {
      return("Missing_Classification") # Handle missing values in Variant_Classification
    } else if (var_type == "DEL" && var_class == "frameshift") {
      return("Frame_Shift_Del") # Frameshift deletion
    } else if (var_type == "INS" && var_class == "frameshift") {
      return("Frame_Shift_Ins") # Frameshift insertion
    } else if (var_class == "nonsynonymous") {
      if (is.na(ref) || is.na(alt)) {
        return("Missing_Allele") # Handle missing values in Allele columns
      } else if (nchar(ref) != nchar(alt)) {
        if (nchar(ref) > nchar(alt)) {
          return("In_Frame_Del") # In-frame deletion
        } else {
          return("In_Frame_Ins") # In-frame insertion
        }
      } else if (var_type == "SNP") {
        return("Missense_Mutation") # SNP classified as Missense Mutation
      }
    } else if (var_class == "synonymous") {
      return("Silent") # Synonymous mutations as Silent
    } else if (var_class == "nonsense") {
      return("Nonsense_Mutation") # Nonsense mutations
    }
    return("Other") # Other cases or unchanged rows
  }
  matS$Variant_Classification <- mapply(get_mutation_type, matS$ref_allele, matS$alt_allele, matS$Type, matS$aaType)
  matS <- matS[, -which(names(matS) == "aaType")]
  return(matS)
  }