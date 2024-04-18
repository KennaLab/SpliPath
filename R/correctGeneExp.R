#' Correcting for gene expression level
#'
#' @param gene_exp string The file of gene expression table with gene IDs in row and sample IDs in column. The first line of the file should be the samples IDs and the first column should be gene IDs.
#' @param rna_meta string or data.frame The meta data for RNAseq samples. It should contain at least 4 columns: SubjectID, SampleID, Group, Tissue. The confounders that will be used to correct for gene expression level should be in the columns. 
#' @param confounder vector The confounders that should be corrected. Should be in the columns of the meta data.
#' @param output_prefix string Prefix of output file. The output files names will be "prefix_tissue_gene_exp_zscore.txt.gz"
#'
#' @export
correctGeneExp -> 
  function(gene_exp, rna_meta, confounder, output_prefix){
    
    gene_exp = read.table(gene_exp, header=T, row.names=1, stringsAsFactors = F, check.names = F)
    rna_meta = read.table(rna_meta, header=T, sep="\t", stringsAsFactors = F)
    
    if (sum(!rna_meta$SampleID %in% colnames(gene_exp)) > 0){
      stop(sprintf("Error: %s sample IDs in the metadata were not found in the gene expression tables"), sum(!rna_meta$SampleID %in% colnames(gene_exp)) )
    }
    
    tissueGeneExpCorr=function(tissue){ 
      
      sample_meta = rna_meta[rna_meta$Tissue == tissue, ]
      tissue_gene_exp = gene_exp[, sample_meta$SampleID]
      
      
      
      
      
      
      output_gz = gzfile(sprintf("%s_novel_junction_exp_%s.txt.gz", output_prefix, tissue), "w")
      write.table(tissue_map[, c("Subject", "Coordinates_of_novel_junc", "Read", "BG", "PSI", "P_binom", "Gene_expression_Z")], output_gz, sep="\t", quote=F, row.names = F)
      close(output_gz)
    }
    
    
    for (tissue in unique(rna_meta$Tissue)){
      tissueGeneExpCorr(tissue)
    }
    
    
  }