# Remove duplicates based on Consequence Score and smallest relative Position in gene  
remove_duplicates <- function(cadd_data, id = Ident, pos = relcDNApos) {
  cadd_data |> dplyr::group_by({{ id }}) |>
    dplyr::slice_max(ConsScore, with_ties = TRUE) |>
    dplyr::slice_min({{ pos }}, with_ties = TRUE) |>
    dplyr::arrange(Chrom, Pos) |> 
    dplyr::ungroup()
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


# Parse gnomad column
parse_gnomad <- function(gnomad_data, column = GnomAD_Exomes){
  gnomad_data |> 
    tidyr::separate_longer_delim({{ column }}, ";") |>
    tidyr::separate_wider_delim({{ column }}, "=", names = c("key", "value"),
                                too_few = "align_start", too_many = "merge") |>
    tidyr::pivot_wider(names_from = key, values_from = value) |>
    utils::type.convert(as.is = TRUE)
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
