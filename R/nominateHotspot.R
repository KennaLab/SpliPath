#' Nominate hotspots of splice-altering variants from variants that mapped to novel junctions
#'
#' @param tissue_var2junc string or data.frame The output The output file or data.frame of mapVariantJunc function. One row in the the dataframe is an annotated DNA variant that mapped to a novel junction in a subject's tissues. 
#' @param min_psi number The minimum PSI required for observed novel junction.
#' @param min_read number The minimum number of reads required to support the observed novel junction. 
#' @param min_read_include number Also include junctions supported by more than min_read_include reads even of their PSI are under min_psi.
#' @param tissues vector The names of tissues.
#' @param min_nr_tissue number Require novel junction observation concordant across min_nr_tissue tissues.
#' @param min_splice_pred list The minimum prediction score of prediciton tools. Default: list(SpliceAI = 0.2, dbscSNV = 0.7).
#' @param max_allele_freq string The maximum allele frequency.
#' @param splicai_pred_match string or data.frame Optional. The output file or dataframe of matchPredNObs function, which shows whether the SpliceAI prediction of a variant match an observed novel junction. If provided, information of whether predictions match observations will be add in the output file.
#' @param blacklist string or NULL Optional. A genomic region blacklist. If provided, the regions will be removed.  Default: system.file("extdata/Reference/blacklist/hg38-blacklist-nochr.v2.bed.gz", package="SpliPath").
#' @param output_prefix string The nominated mutation hotspot will be written in output_prefix_hotspot.txt.gz. The "Hotspot" column shows whether the variants and novel junctions are passed the given thresholds and nominated as mutation hotspots.
#'
#' @import tidyr
#' @import dplyr
#' @import rtracklayer
#' @export

nominateHotspot <-
function(tissue_var2junc, 
         min_psi,
         min_read,
         min_read_include,
         tissues,
         min_nr_tissue,
         min_splice_pred = list(SpliceAI = 0.2, dbscSNV = 0.7),
         max_allele_freq = 0.05,
         splicai_pred_match = NULL,
         blacklist = system.file(paste("extdata", "Blacklist", "hg38-blacklist-nochr.v2.bed.gz", sep=.Platform$file.sep), package = "SpliPath"),
         output_prefix = ""){
  
  if (!is.data.frame(tissue_var2junc)){
    tissue_var2junc = read.table(tissue_var2junc, header=T, sep="\t", stringsAsFactors = F)
  }
  
  if (!is.null(blacklist)){
    blacklist = rtracklayer::import(blacklist, format="bed")
    tissue_var2junc = tidyr::separate(tissue_var2junc, col="Var_region", sep=":", into=c("chr", "start", "end"), remove = F)
    tissue_var2junc = GenomicRanges::makeGRangesFromDataFrame(tissue_var2junc, seqnames.field = "chr", start.field = "start", end.field = "end", ignore.strand = T, keep.extra.columns = TRUE)
    overlap = GenomicRanges::findOverlaps(blacklist, tissue_var2junc, ignore.strand=T)
    if(length(S4Vectors::subjectHits(overlap))>0){
      tissue_var2junc = tissue_var2junc[-unique(S4Vectors::subjectHits(overlap))]
    }
    tissue_var2junc = as.data.frame(tissue_var2junc, stringsAsFactors = F, row.names = NULL)
    tissue_var2junc = subset(tissue_var2junc, select=-c(seqnames, start, end, width, strand))
  }
  
  if (!is.null(splicai_pred_match)){
    if (!is.data.frame(splicai_pred_match)){
      splicai_pred_match = read.table(splicai_pred_match, sep="\t", header=T, stringsAsFactors = F)
    }
    tissue_var2junc$idx = paste(tissue_var2junc$VAR_id, tissue_var2junc$Coordinates_of_novel_junc, sep=":")
    tissue_var2junc$SpliceAI_pred_match_junction = splicai_pred_match[match(tissue_var2junc$idx, splicai_pred_match$idx), "match_by_threshold"]
    tissue_var2junc$SpliceAI_pred_cryptic_exon = splicai_pred_match[match(tissue_var2junc$idx, splicai_pred_match$idx), "pred_cryptic_exon"]
  }

  ### Nominate mutation hotspot and sRV by thresholds
  psi_colnames = paste0("PSI_", tissues)
  read_colnames = paste0("Read_", tissues)
  psi_pass1 = tissue_var2junc[, psi_colnames] >= min_psi
  read_pass1 = tissue_var2junc[, read_colnames] >= min_read

  read_pass2 = tissue_var2junc[, read_colnames] >= min_read_include
  
  PASS = rowSums((psi_pass1 & read_pass1) | read_pass2, na.rm = T) >= min_nr_tissue
  
  PASS = PASS & (tissue_var2junc$AF <= max_allele_freq)

  splice_pass = as.data.frame(matrix(NA, nrow=nrow(tissue_var2junc), ncol=length(names(min_splice_pred))), col.names = names(min_splice_pred))
  for (pred in names(min_splice_pred)){
    splice_pass[[pred]] = tissue_var2junc[[pred]] >= min_splice_pred[[pred]]
  }
  splice_pass = rowSums(splice_pass[, names(min_splice_pred)], na.rm = T) > 0
  
  PASS = PASS & splice_pass

  tissue_var2junc$Hotspot = F
  tissue_var2junc[PASS, "Hotspot"] = T
  
  output_file = gzfile(sprintf("%s_hotspot.txt.gz", output_prefix), "w")
  write.table(tissue_var2junc, output_file, sep="\t", quote=F, row.names = F)
  close(output_file)
}
