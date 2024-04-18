#' Find annotated junction that share a same annotated splice site with a novel junction
#' 
#' @return coordinates of junctions concatenated by "|"
#' @export
annoJuncWithShareSites <- 
  function(gene_junc, index, docking_site){
  docking_pos = gene_junc[index, docking_site]
  count_canonical = gene_junc[(gene_junc[[docking_site]] == docking_pos) , ]
  
  if (nrow(count_canonical) == 0){
    NA
  }else{
    paste(rownames(count_canonical), collapse="|")
  }
}
