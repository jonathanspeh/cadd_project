# Remove duplicates based on Consequence Score and smallest relative Position in gene  
remove_duplicates <- function(cadd_data, id = Ident, pos = relcDNApos) {
  cadd_data |> dplyr::group_by({{ id }}) |>
    dplyr::slice_max(ConsScore, with_ties = TRUE) |>
    dplyr::slice_min({{ pos }}, with_ties = TRUE) |>
    dplyr::arrange(Chrom, Pos) |> 
    ungroup()
  }




#TODO - Process gnomad overlaps
# process_overlaps <- function(tsv, keep = keep){}


# Parse gnomad column
parse_gnomad <- function(gnomad_data, column = GnomAD_Exomes){
  gnomad_data |> 
    tidyr::separate_longer_delim({{ column }}, ";") |>
    tidyr::separate_wider_delim({{ column }}, "=", names = c("key", "value"),
                                too_few = "align_start", too_many = "merge") |>
    tidyr::pivot_wider(names_from = key, values_from = value) |>
    utils::type.convert(as.is = TRUE)
  }




