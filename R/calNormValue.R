#' Calculate PSI
#' 
#' The number of read support a novel junction divided by the maximum number of read support the annotated junction that share a same splice site with the novel junction.
#' 
#' @return number PSI
#' @export
calNormValue <-
function(gene_junc, index, docking_site){
  docking_pos = gene_junc[index, docking_site]
  count = gene_junc[index, 1]
  count_canonical = gene_junc[(gene_junc[[docking_site]] == docking_pos) & (gene_junc[["coding.transcript"]] == TRUE), ]
  if (nrow(count_canonical) == 0){
    normalized_val = 1 
  }else{
    normalized_val = count / (max(count_canonical[, 1]) + count) 
  }
  normalized_val
}
