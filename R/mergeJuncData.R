#' Merge junction data of samples to a junction - sample data frame
#'
#' @param sample_path_file string Name of tab-delimited file which contain two columns: 1, unique sample ids; 2, paths to the junction files.
#' @param chrom string The junctions in this chromosome will be merged among samples.
#' @param output_prefix string Prefix of output file name. The junction-sample table will be written in prefix_chr*.txt.gz.
#' @import dplyr
#' 
#' @return A data frame of coordinates of unique junctions in four columns: chr, start, end, and strand.
#' @export
mergeJuncData <-
function(sample_path_file, chrom, output_prefix = ""){
  
  readJuncFile = function(sample_id, file_name, chrom){

    ### Read junctions in regtools BED format
    junc = read.table(text = gsub(",", "\t", readLines(file_name)), header = F, sep = "\t", stringsAsFactors = F)
    colnames(junc)[c(1,5,6)] = c("chr", sample_id, "strand")
    junc$start = junc$V2 + junc$V13
    junc$end = junc$V3 - junc$V14
    
    junc = junc[(junc$chr == chrom) & (junc$strand %in% c("+", "-")),]
    return(junc[,c('chr','start', 'end', 'strand', sample_id)])
  }
  
  sample_list = read.table(sample_path_file, sep="\t", header=T, stringsAsFactors = F)
  merged_junc = readJuncFile(sample_list[1, "SampleID"], sample_list[1, "Path"], chrom)
  if (nrow(sample_list) > 1){
    for (idx in 2:nrow(sample_list)){
      sample_junc = readJuncFile(sample_list[idx, "SampleID"], sample_list[idx, "Path"], chrom)
      merged_junc = dplyr::full_join(merged_junc, sample_junc, by = c('chr','start', 'end', 'strand'))
    }
  }
  merged_junc[is.na(merged_junc)] = 0
  merged_junc = merged_junc[order(as.numeric(merged_junc$start), as.numeric(merged_junc$end)), ]
  junc_table_gz = gzfile(sprintf("%s_chr%s.txt.gz", output_prefix, chrom), "w")
  write.table(format(merged_junc, scientific=FALSE), junc_table_gz, row.names = F, col.names = T, sep = "\t", quote = F)
  close(junc_table_gz)
  
  merged_junc[, c('chr', 'start', 'end', 'strand')]
}

