#' Calculate percentage spliced in (PSI) for novel junctions 
#'
#' @param junc_file string or data.frame The output file of mergeJuncData function, which contain junction-sample read counts table. The first four columns should be the coordinates of each junction: "chr", "start", "end", "strand".
#' @param junc_anno string or data.frame The output file of annotateJunc function, which contain junction annotations.
#' @param annotation vector Annotation sources (should in colnames in junc_anno) that were combined to determine novel junction. A junction not annotated in any given source is treated as novel junction. Default: c("in.ensembl", "in.snaptron")
#' @param output_prefix string Prefix of output PSI file. The junction-sample PSI table will be written in a prefix_norm.txt.gz file.
#'
#' @return data.frame Seven-column dataframe. Contain coordinates, mapped gene, and splicing event of each novel junction.
#' @export

calJuncPSI <-
function(junc_file, junc_anno, output_prefix, annotation = c("in.ensembl", "in.snaptron")){
  if (is.character(junc_file)){
    junc_file = read.table(junc_file, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
  }
  if (is.character(junc_anno)){
    junc_anno = read.table(junc_anno, sep = "\t", header = T, stringsAsFactors = F)
  }
  samples = colnames(junc_file)[5:ncol(junc_file)]
  junc_anno_col = colnames(junc_anno)
  
  ### Add annotation columns to junction - sample dataframe
  junc_file$idx = paste(junc_file[["chr"]], junc_file[["start"]], junc_file[["end"]], junc_file[["strand"]], sep=":")
  junc_anno$idx = paste(junc_anno[["chr"]], junc_anno[["start"]], junc_anno[["end"]], junc_anno[["strand"]], sep=":")
  junc_file = junc_file[junc_file[['idx']] %in% junc_anno[['idx']], ]
  junc_file = subset(junc_file, select = -idx)
  junc_anno = subset(junc_anno, select = -idx)
  junc_file = dplyr::inner_join(junc_file, junc_anno, by=c('chr', 'start', 'end', 'strand'), multiple = "all")
  junc_file[, annotation] = sapply(junc_file[, annotation], FUN = as.logical)
  junc_file$is.anno = rowSums(junc_file[, annotation], na.rm = T) > 0
  
  ### Cal PSI for novel junctions
  junc_file$norm_val = 0
  norm_val_mtx = junc_file[, c('chr', 'start', 'end', 'strand', 'gene.id', 'gene.name', annotation, 'is.anno', 'event')] #
  for (sample in samples){
    norm_val_mtx[[sample]] = normCount(junc_file[, c(sample, 'norm_val', junc_anno_col, "is.anno")])
  }
  norm_val_mtx = norm_val_mtx[(norm_val_mtx[['is.anno']] == FALSE), ]
  
  psi_table_gz = gzfile(sprintf("%s_norm.txt.gz", output_prefix), "w")
  write.table(format(norm_val_mtx, scientific=FALSE), psi_table_gz, row.names = F, col.names = T, sep = "\t", quote = F)
  close(psi_table_gz)
  
  norm_val_mtx[, c("chr", "start", "end", "strand", "gene.id", "gene.name", "event")]
}
