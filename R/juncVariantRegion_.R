#' Find the intron disrupted by the unannotated junction and the two flanking exons
#' 
#' @export
juncVariantRegion_ <-
function(junc, novel_junc, exon){
  exon_gene = exon[exon$gene.id == novel_junc[junc, "gene.id"],]
  if (((novel_junc[junc, "strand"] == "+") & (novel_junc[junc, "event"] == "novel_acceptor")) | ((novel_junc[junc, "strand"] == "-") & (novel_junc[junc, "event"] == "novel_donor"))){
    upstream_exon = exon_gene[exon_gene$end == novel_junc[junc, "start"],]
    upstream_exon$dist = novel_junc[junc, "start"] - upstream_exon$start
    upstream_exon = upstream_exon[upstream_exon$dist == max(upstream_exon$dist),][1,]
    expand_start = upstream_exon$start
    
    downstream_exon = exon_gene[(exon_gene$start != novel_junc[junc, "end"]) & (exon_gene$end > novel_junc[junc, "end"]),]
    if (nrow(downstream_exon) == 0){downstream_exon = exon_gene[(exon_gene$start == novel_junc[junc, "end"]),]}
    downstream_exon$start_dist = downstream_exon$start - novel_junc[junc, "end"]
    downstream_exon$end_dist = downstream_exon$end - novel_junc[junc, "end"]
    downstream_exon = downstream_exon[abs(downstream_exon$start_dist) == min(abs(downstream_exon$start_dist)),]
    downstream_exon = downstream_exon[downstream_exon$end_dist == max(downstream_exon$end_dist),][1,]
    expand_end = downstream_exon$end
    if (is.na(expand_end)){expand_end = novel_junc[junc, "end"]}
    
  }else if(((novel_junc[junc, "strand"] == "+") & (novel_junc[junc, "event"] == "novel_donor")) | ((novel_junc[junc, "strand"] == "-") & (novel_junc[junc, "event"] == "novel_acceptor"))){
    upstream_exon = exon_gene[(exon_gene$end != novel_junc[junc, "start"]) & (exon_gene$start < novel_junc[junc, "start"]),]
    if (nrow(upstream_exon) == 0){upstream_exon = exon_gene[(exon_gene$end == novel_junc[junc, "start"]),]}
    upstream_exon$end_dist = novel_junc[junc, "start"] - upstream_exon$end 
    upstream_exon$start_dist = novel_junc[junc, "start"] - upstream_exon$start
    upstream_exon = upstream_exon[abs(upstream_exon$end_dist) == min(abs(upstream_exon$end_dist)),]
    upstream_exon = upstream_exon[upstream_exon$end_dist == max(upstream_exon$end_dist),][1,]
    expand_start = upstream_exon$start
    if (is.na(expand_start)){expand_start = novel_junc[junc, "start"]}
    
    downstream_exon = exon_gene[exon_gene$start == novel_junc[junc, "end"],]
    downstream_exon$dist = downstream_exon$end - novel_junc[junc, "end"]
    downstream_exon = downstream_exon[downstream_exon$dist == max(downstream_exon$dist),][1,]
    expand_end = downstream_exon$end
    
  }else if((novel_junc[junc, "event"] == "exon_skipping") | (novel_junc[junc, "event"] == "alternative_intron")){
    upstream_exon = exon_gene[exon_gene$end == novel_junc[junc, "start"],]
    upstream_exon$dist = novel_junc[junc, "start"] - upstream_exon$start
    upstream_exon = upstream_exon[upstream_exon$dist == max(upstream_exon$dist),][1,]
    expand_start = upstream_exon$start
    if (is.na(expand_start)){expand_start = novel_junc[junc, "start"]}
    
    downstream_exon = exon_gene[exon_gene$start == novel_junc[junc, "end"],]
    downstream_exon$dist = downstream_exon$end - novel_junc[junc, "end"]
    downstream_exon = downstream_exon[downstream_exon$dist == max(downstream_exon$dist),][1,]
    expand_end = downstream_exon$end
    if (is.na(expand_end)){expand_end = novel_junc[junc, "end"]}
    
  }
  c(unlist(novel_junc[junc,]), expand_start, expand_end)
}
