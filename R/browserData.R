#' Prepare data files for SpliPath data browser
#' 
#' Make splice junction tables and csQTL tables for each gene. 
#' 
#' @param paired_data logical "TRUE" for paired data mode or "FALSE" for unpaired data mode. If TRUE, arguments tissues_leafcutter_pvals_file, junc_count_files, rna_meta, and wgs_meta are required. The function matches the paired samples by 'SubjectID' in the metadata.   
#' @param junc_anno_files vector A list of output files of annotateJunc function, which contain junction annotations.
#' @param cs_qtl_files string The output file of collapsedsQTL function. One row in the the data.frame is an annotated DNA variant that mapped to a novel junction in a subject's tissues. 
#' @param tissues_leafcutter_pvals_file list The file paths to the LeafCutter analysis P values of each tissue type. Each element is a LeafCutter file path named under the tissue type. The tissues should be the same with those in the RNAseq sample metadata file.
#' @param junc_count_files vector A list of output files of mergeJuncData funciton, which contain junction-sample read count table. The input file name list in the first four arguments should be matched.
#' @param output_dir string The path to data output directory.
#' 
#' @import foreach
#' @import ggplot2
#' @import stringr

#' @export
browserData<-
  function(paired_data, junc_anno_files, cs_qtl_files, tissues_leafcutter_pvals_file = NULL, junc_count_files = NULL, output_dir = "."){
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    
    # total_novel_junc = data.frame()
    
    if (paired_data){
      leafcutter_pvals_merged = NULL
      for (tissue in names(tissues_leafcutter_pvals_file)){
        
        leafcutter_pvals_file = tissues_leafcutter_pvals_file[[tissue]]
        leafcutter_pvals = read.table(leafcutter_pvals_file, header = T, sep="\t", row.names = 1, stringsAsFactors = F, check.names=F)
        colnames(leafcutter_pvals) = do.call(rbind, strsplit(colnames(leafcutter_pvals), split = ".Aligned"))[, 1]
        rename_row = do.call(rbind, strsplit(rownames(leafcutter_pvals), "\\:|\\_"))[, c(1:3,6)]
        colnames(rename_row) = c("chr", "start", "end", "strand")
        rename_row[, "end"] = as.integer(rename_row[, "end"]) - 1
        rownames(leafcutter_pvals) = apply(rename_row, 1, paste, collapse=":")
  
        if (is.null(leafcutter_pvals_merged)){
          leafcutter_pvals_merged = leafcutter_pvals
        }else{
          leafcutter_pvals_merged = merge(leafcutter_pvals_merged, leafcutter_pvals, by=0, all=TRUE)
          rownames(leafcutter_pvals_merged) = leafcutter_pvals_merged$Row.names
          leafcutter_pvals_merged = leafcutter_pvals_merged[, -1]
        }
        
      }
    }
    for (n_file in 1:length(junc_anno_files)){
      
      if (paired_data){
        junc_count = read.table(junc_count_files[n_file], header=T, stringsAsFactors = F, check.names = F)
        rownames(junc_count) = paste(junc_count$chr, junc_count$start, junc_count$end, junc_count$strand, sep=":")
      } 
      junc_anno = read.table(junc_anno_files[n_file], header=T, stringsAsFactors = F, check.names = F)
      sRV = read.table(cs_qtl_files[n_file], header=1, stringsAsFactors = F, check.names = F)
      sRV$SpliceAI_pred_match_junction = sRV$match_by_threshold
        

      gene_list = unique(junc_anno[, c("gene.id", "gene.name")])
      for (g in 1:nrow(gene_list)){
        gene_id = gene_list[g, "gene.id"]
        gene_name = toupper(gene_list[g, "gene.name"])
        
        junc_anno_gene = subset(junc_anno, gene.id == gene_id)
        
        if (paired_data){
          gz_file = gzfile( paste(output_dir, sprintf("%s_%s_intron_count.txt.gz", gene_id, gene_name), sep=.Platform$file.sep), "w")
          write.table(cbind(junc_anno_gene[, c("pos", "gene.id", "gene.name", "event")], junc_count[junc_anno_gene$pos, ]), gz_file, row.names = F, col.names=T, sep='\t', quote=F)
          close(gz_file)
        
          gz_file = gzfile( paste(output_dir, sprintf("%s_%s_leafcutter_pval.txt.gz", gene_id, gene_name), sep=.Platform$file.sep), "w")
          write.table(leafcutter_pvals_merged[(rownames(leafcutter_pvals_merged) %in% junc_anno_gene$pos), ], gz_file, row.names = T, col.names=T, sep='\t', quote=F)
          close(gz_file)
        }
        gz_file = gzfile( paste(output_dir, sprintf("%s_%s_csQTL.txt.gz", gene_id, gene_name), sep=.Platform$file.sep), "w")
        write.table(subset(sRV, Gene_id == gene_id), gz_file, row.names = F, col.names=T, sep='\t', quote=F)
        close(gz_file)
        
      }
    }

  }
