#' Prepare data files for SpliPath data browser
#' 
#' Separate splice junction tables and sRV tables per gene. The separated files will be written in SpliPath_Browser/data/junction and SpliPath_Browser/data/sRV_candidate
#'
#' @param junc_count vector A list of output files of mergeJuncData funciton, which contain junction-sample read count table. The input file name list in the first four arguments should be matched.
#' @param junc_psi vector A list of output files of calJuncPSI function, which contain junction-sample PSI table.
#' @param junc_anno vector A list of output files of annotateJunc function, which contain junction annotations.
#' @param hotspot vector A list of output files of nominateHotspot function. 
#' @param output_dir string The path to data output directory.
#' 
#' @import foreach
#' @import ggplot2
#' @import stringr
#' @import ggVennDiagram

#' @export
browserData <-
function(junc_count, junc_psi, junc_anno, hotspot, min_read = 5, output_dir = "."){
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  
  total_novel_junc = data.frame()
  
  for (n_file in 1:length(junc_count)){
    
    junc_anno = read.table(junc_anno[n_file], header=T, stringsAsFactors = F, check.names = F)
    junc_count = read.table(junc_count[n_file], header=T, stringsAsFactors = F, check.names = F)
    rownames(junc_count) = paste(junc_count$chr, junc_count$start, junc_count$end, junc_count$strand, sep=":")
    junc_psi = read.table(junc_psi[n_file], header=T, stringsAsFactors = F, check.names = F)
    sRV = read.table(hotspot[n_file], header=1, stringsAsFactors = F, check.names = F)
    
    if (nrow(total_novel_junc) == 0){
      total_novel_junc = data.frame(colSums(junc_count[junc_anno[junc_anno$event %in% c("novel_acceptor", "novel_donor", "exon_skipping"), "pos"], 5:ncol(junc_count)] >= min_read))
      rownames(total_novel_junc) = colnames(junc_count)[5:ncol(junc_count)]
      colnames(total_novel_junc) = sprintf("nr_novel_junction_min_%s_reads", min_read)
    }else{
      total_novel_junc = total_novel_junc + rowSums(junc_count[junc_anno$pos, 5:ncol(junc_count)] >= 5)
    }
    
    gene_list = unique(junc_anno[, c("gene.id", "gene.name")])
    for (g in 1:nrow(gene_list)){
      gene_id = gene_list[g, "gene.id"]
      gene_name = toupper(gene_list[g, "gene.name"])
  
      junc_anno_gene = subset(junc_anno, gene.id == gene_id)
      
      gz_file = gzfile( paste(output_dir, sprintf("%s_%s_intron_count.txt.gz", gene_id, gene_name), sep=.Platform$file.sep), "w")
      write.table(cbind(junc_anno_gene[, c("pos", "gene.id", "gene.name", "event")], junc_count[junc_anno_gene$pos, ]), gz_file, row.names = F, col.names=T, sep='\t', quote=F)
      close(gz_file)
      
      gz_file = gzfile( paste(output_dir, sprintf("%s_%s_intron_psi.txt.gz", gene_id, gene_name), sep=.Platform$file.sep), "w")
      write.table(subset(junc_psi , gene.id == gene_id), gz_file, row.names = F, col.names=T, sep='\t', quote=F)
      close(gz_file)
      
      gz_file = gzfile( paste(output_dir, sprintf("%s_%s_sRV_candidate.txt.gz", gene_id, gene_name), sep=.Platform$file.sep), "w")
      write.table(subset(sRV, Gene_id == gene_id), gz_file, row.names = F, col.names=T, sep='\t', quote=F)
      close(gz_file)
    }
  }
  write.table(total_novel_junc, paste(output_dir, "number_of_novel_junction_per_subject.txt", sep=.Platform$file.sep), row.names = T, col.names=T, sep='\t', quote=F)
}
