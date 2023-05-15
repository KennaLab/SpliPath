#' Calculate PSI of novel junctions within each gene
#' 
#' A function to connect calJuncPSI function and normCount_ function. If there are only novel junctions in a gene, the PSI will not be calculated (set to 0 by calJuncPSI function).
#' 
#' @return data.frame 
#' @export
normCount <-
function(junc_count){
  genes = unique(junc_count[(junc_count[, 1] > 0) & (junc_count$is.anno == FALSE) , 'gene.id'])
  for (gene in genes){
    gene_junc = junc_count[(junc_count$gene.id == gene) & (junc_count[, 1] > 0), ]
    gene_junc$norm_val = NA
    if (nrow(gene_junc) == 1){
      next
    }
    if (length(unique(gene_junc$is.anno)) > 1){
      gene_junc = normCount_(gene_junc)
      junc_count[rownames(gene_junc), 'norm_val'] = gene_junc$norm_val
    }
  }
  round(junc_count$norm_val, 4)
}
