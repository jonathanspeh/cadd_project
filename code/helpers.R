# Remove duplicates based on Consequence Score and smallest relative Position in gene  
remove_duplicates <- function(cadd_data, id = Ident, pos = relcDNApos) {
  cadd_data |> dplyr::group_by({{ id }}) |>
    dplyr::slice_max(ConsScore, with_ties = TRUE) |>
    dplyr::slice_min({{ pos }}, with_ties = TRUE) |>
    dplyr::arrange(Chrom, Pos) |> 
    dplyr::ungroup()
  }

remove_duplicates_dt <- function(cadd_data, id = "Ident") {
  require(data.table)
  cadd_data <-  as.data.table(cadd_data)
  cadd_data <- cadd_data[cadd_data[, .I[ConsScore == max(ConsScore)], by = id]$V1]
  cadd_data[cadd_data[, .I[relcDNApos == min(relcDNApos)], by = id]$V1]
}






# Process Overlaps with gnomad
process_caddscoring <- function(.data, cols = keep, autosomes_only = TRUE){
  processed <- .data |>
    dplyr::select(dplyr::all_of({{ cols }})) |>
    dplyr::mutate("ChromPos" = paste(Chrom, Pos, sep=":"), 
                  "Ident" = paste(ChromPos, Ref, Alt, sep = ":"))
  if(autosomes_only) processed <- dplyr::filter(processed, !Chrom %in% c("X", "Y"))
  return(processed)
}

process_caddscoring_dt <- function(.data, cols = keep, autosomes_only = TRUE){
  DT <- data.table(.data)
  colnames(DT)[which(colnames(DT)=="#Chrom")] <- "Chrom"
  colnames(DT)[which(colnames(DT)=="GenomAD-Exomes")] <- "GnomAD_Exomes"
  DT <- DT[, ..cols]
  DT[, Ident := paste0(Chrom, ":", Pos, Ref, ">", Alt)]
  if(autosomes_only) DT <- DT[!Chrom %in% c("X", "Y")]
  tibble(DT)
}


# Parse gnomad column
# Can be used to get all gnomad annotations, however is quite slow and heavy on memory
# Needs optimisation
parse_gnomad <- function(gnomad_data, column = GnomAD_Exomes){
  gnomad_data |> 
    tidyr::separate_longer_delim({{ column }}, ";") |>
    tidyr::separate_wider_delim({{ column }}, "=", names = c("key", "value"),
                                too_few = "align_start", too_many = "merge") |>
    tidyr::pivot_wider(names_from = key, values_from = value) |>
    utils::type.convert(as.is = TRUE)
  }

# Gets AC, AN and AF of gnomad data
# Faster and uses way less memory than parse_gnomad()
get_afs <- function(gnomad_data) {
  require(data.table)
  require(dplyr)
  gnomad_data <- as.data.table(gnomad_data)
  gnomad_data <- cbind(gnomad_data[, !"GnomAD_Exomes"], gnomad_data[, tstrsplit(GnomAD_Exomes, split = ";")[1:3]])
  gnomad_data <- mutate(gnomad_data, across(starts_with("V"), function(x) sub(".*=", "", x)))
  colnames(gnomad_data)[which(colnames(gnomad_data)=="V1")] <- "AC"
  colnames(gnomad_data)[which(colnames(gnomad_data)=="V2")] <- "AN"
  colnames(gnomad_data)[which(colnames(gnomad_data)=="V3")] <- "AF"
  gnomad_data <- utils::type.convert(gnomad_data, as.is = TRUE)
  as_tibble(gnomad_data)
}




# Print SNV information
print_snv_info <- function(data, row){
  data <- data <- data[row,]
  out <- data|>
    mutate("Position" = paste0(Chrom, ":", Pos)) |>
    dplyr::select(all_of(c("Position", "Ref", "Alt", "Allele Count" = "AC", 
                       "Allele Frequency" = "AF", "PHRED", "Gene" = "GeneName"))) |>
    t() |> data.frame() 
  
  colnames(out) <- " "
  out
  }
