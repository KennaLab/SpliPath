#' Annotate junctions 
#'
#' @import dplyr
#' @import GenomicRanges
#' @return data.frame One row in the the data.frame is annotations of a unique junction. 
annotateJunc_ <-
  function(merged_junc, gene_ref, ensembl_pos, protein_coding_gene, intron_coding, intron_noncoding, addition_annotation){
    
    options(scipen = 999)

    merged_junc = GenomicRanges::makeGRangesFromDataFrame(merged_junc, keep.extra.columns = T, starts.in.df.are.0based=T)
    gene_ref = GenomicRanges::makeGRangesFromDataFrame(gene_ref, keep.extra.columns = T, starts.in.df.are.0based=T)
    overlaps = GenomicRanges::findOverlaps(merged_junc, gene_ref)
    
    merged_junc = data.frame(merged_junc[data.frame(overlaps)$queryHits])
    merged_junc[, c("seqnames", "strand")] = sapply(merged_junc[, c("seqnames", "strand")], as.character)
    merged_junc = merged_junc[, c("seqnames", "start", "end", "strand")]
    colnames(merged_junc) = c("chr",  "start", "end", "strand")
    
    gene_ref = data.frame(gene_ref[data.frame(overlaps)$subjectHits])
    gene_ref[, c("seqnames", "strand")] = sapply(gene_ref[, c("seqnames", "strand")], as.character)
    gene_ref = gene_ref[, c("seqnames", "start", "end", "gene.id", 'gene.name', "strand")]
    colnames(gene_ref) = c("chr.map", "start.map", "end.map", "gene.id", 'gene.name', "strand.map")
    
    junc_map = cbind(merged_junc, gene_ref)
    junc_map[, c("start", "start.map")] = junc_map[, c("start", "start.map")] - 1 ### Make coords 0-based
    rm(merged_junc, gene_ref)
    # End of new lines
    
    junc_map$pos = paste(junc_map$chr, junc_map$start, junc_map$end, junc_map$strand, sep=':')
    length(unique(junc_map$pos[duplicated(junc_map$pos)]))/length(unique(junc_map$pos))
    junc_map$start.pos = paste(junc_map$chr, junc_map$start, junc_map$strand, junc_map$gene.id, sep=':')
    junc_map$end.pos = paste(junc_map$chr, junc_map$end, junc_map$strand, junc_map$gene.id, sep=':')
    
    if (!is.null(addition_annotation)){
      for (addition_anno in names(addition_annotation)){
        addition_anno_pos = read.table(addition_annotation[[addition_anno]], header=F, stringsAsFactors = F)$V1
        anno_colname = paste0("in.", addition_anno)
        junc_map[[anno_colname]] = junc_map$pos %in% addition_anno_pos
      }
    }
    junc_map$in.ensembl = junc_map$pos %in% ensembl_pos
    junc_map$coding.gene = junc_map$gene.id %in% protein_coding_gene
    junc_map$coding.transcript = paste(junc_map$pos, junc_map$gene.id, sep= ':') %in% intron_coding$pos
    junc_map$noncoding.transcript = paste(junc_map$pos, junc_map$gene.id, sep=':') %in% intron_noncoding$pos
    
    junc_map$start.anno = junc_map$start.pos %in% intron_coding$start.pos
    junc_map$end.anno = junc_map$end.pos %in% intron_coding$end.pos
    junc_map$start.map = junc_map$start > junc_map$start.map
    junc_map$end.map = junc_map$end < junc_map$end.map
    
    junc_map$map[junc_map$start.map & junc_map$end.map] = "mapped"
    junc_map$map[!junc_map$start.map | !junc_map$end.map] = "d/a_site_unmapped"
    junc_map$map[!junc_map$start.map & !junc_map$end.map] = "unmapped"
    
    junc_map$event = 'annotated'
    junc_map$event[!junc_map$start.anno & junc_map$end.anno & (junc_map$strand == "+")] = "novel_donor"
    junc_map$event[!junc_map$start.anno & junc_map$end.anno & (junc_map$strand == "-")] = "novel_acceptor"
    junc_map$event[junc_map$start.anno & !junc_map$end.anno & (junc_map$strand == "+")] = "novel_acceptor"
    junc_map$event[junc_map$start.anno & !junc_map$end.anno & (junc_map$strand == "-")] = "novel_donor"
    junc_map$event[junc_map$start.anno & junc_map$end.anno & !junc_map$coding.transcript & !junc_map$noncoding.transcript] = "exon_skipping"
    junc_map$event[!junc_map$coding.transcript & junc_map$noncoding.transcript] = "noncoding_transcript"
    junc_map$event[!junc_map$start.anno & !junc_map$end.anno & !junc_map$noncoding.transcript & !junc_map$in.ensembl] = "novel_d&a_site"
    
    exon_skipping_idx = which(junc_map$event == "exon_skipping")
    exon_skipping_anno = lapply(exon_skipping_idx, FUN = function(idx, junc_map){
      ref_intron_gene = intron_coding[intron_coding$gene.id == junc_map[idx, "gene.id"],]
      start = junc_map[idx, "start"]
      end = junc_map[idx, "end"]
      start_transcript = ref_intron_gene[ref_intron_gene$start == start, "transcript.id"]
      end_transcript = ref_intron_gene[ref_intron_gene$end == end, "transcript.id"]
      span_exon_start = TRUE %in% ((ref_intron_gene$end > start) & (ref_intron_gene$end < end))
      span_exon_end = TRUE %in% ((ref_intron_gene$start > start) & (ref_intron_gene$start < end))
      
      if ((length(intersect(start_transcript, end_transcript)) > 0) & (span_exon_start & span_exon_end)){
        return("exon_skipping")
      }else{return("alternative_intron")}
    },
    junc_map = junc_map
    )
    exon_skipping_anno = unlist(exon_skipping_anno)
    junc_map[exon_skipping_idx, "event"] = exon_skipping_anno
    
    # Select junctions mapped within protein-coding genes
    junc_map = junc_map[junc_map$map == "mapped" & junc_map$coding.gene & junc_map$event != "novel_d&a_site",]
    
    # To remove the less likely map when junctions were mapped to multiple genes
    junc_map = junc_map[!(junc_map$coding.transcript == FALSE & junc_map$noncoding.transcript == FALSE & junc_map$in.ensembl == TRUE),]

    junc_map
  }

