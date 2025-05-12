#' Remove the regions that overlap with the provided blacklist regions
#'
#' @param input_region GRanges object The input regions, from which the ones overlaps with the blacklist will be removed.
#' @param blacklist string The name of BED file containing the blacklist regions
#' @param addition_annotation list or NULL Optional. Name and file name of additional junction annotation source. Default: list(snaptron = system.file("extdata/Reference/snaptron_gtex_anno.txt", package = "SpliPath")). Whether the input junctions are the given junction list will be shown in the output file (in column "in.annotation_name").
#'
#' @import GenomicRanges
#' @import rtracklayer
rmBlackListRegion=function(input_region, blacklist){
  if (!is.nll(blacklist)){
    if (class(input_region) == "data.frame"){
      input_region = GenomicRanges::makeGRangesFromDataFrame(input_region, seqnames.field = "CHROM", start.field = "start", end.field = "end", strand.field = "Strand", keep.extra.columns = TRUE)
    }
    blacklist = rtracklayer::import(blacklist, format="bed")
    
    overlap = GenomicRanges::findOverlaps(blacklist, input_region, ignore.strand=T)
    if(length(data.frame(overlap)$subjectHits)>0){
      input_region = input_region[-unique(data.frame(overlap)$subjectHits)]
      message('Splicing junctions in genome blacklist regions were removed')
    }
  }
  input_region
}
