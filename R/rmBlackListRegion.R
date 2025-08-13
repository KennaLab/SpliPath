#' Remove the regions that overlap with the provided blacklist regions
#'
#' @param input_region GRanges object The input regions, from which the ones overlaps with the blacklist will be removed.
#' @param blacklist string The name of BED file containing the blacklist regions
#'
#' @import GenomicRanges
#' @import rtracklayer
rmBlackListRegion=function(input_region, blacklist){
  if (!is.null(blacklist)){
    blacklist = rtracklayer::import(blacklist, format="bed")
    overlap = GenomicRanges::findOverlaps(blacklist, input_region, ignore.strand=T)
    
    if(length(data.frame(overlap)$subjectHits)>0){
      input_region = input_region[-unique(data.frame(overlap)$subjectHits)]
      message('Splicing junctions in genome blacklist regions were removed')
    }
  }
  input_region
}
