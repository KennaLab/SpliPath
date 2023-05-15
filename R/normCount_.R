#' Calculate PSI of different splicing events in a gene
#'
#' The function to connect normCount function and calNormValue function.
#' 
#' @return data.frame
#' @export
normCount_ <-
function(gene_junc){
  novel_junc = rownames(gene_junc[gene_junc$is.anno == FALSE, ])
  for (i in novel_junc){
    if ((gene_junc[i, 'event'] == 'exon_skipping') | (gene_junc[i, 'event'] == 'alternative_intron')){
      normalized_val1 = calNormValue(gene_junc, i, 'start')
      normalized_val2 = calNormValue(gene_junc, i, 'end')
      gene_junc[i, 'norm_val'] = min(normalized_val1, normalized_val2)
    }else if(gene_junc[i, 'start.anno'] == FALSE){
      gene_junc[i, 'norm_val'] = calNormValue(gene_junc, i, 'end')
    }else if(gene_junc[i, 'end.anno'] == FALSE){
      gene_junc[i, 'norm_val'] = calNormValue(gene_junc, i, 'start')
    }
  }
  gene_junc
}
