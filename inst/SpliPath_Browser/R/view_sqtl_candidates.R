
draw_dbscSNV = function(data, tissue, min_read, min_psi, min_score){
  psi_col = paste0("PSI_", tissue)
  rc_col = paste0("Read_", tissue)
  colnames(data)[colnames(data) == psi_col] = "PSI"
  colnames(data)[colnames(data) == rc_col] = "Read"
  
  data = data[!is.na(data$dbscSNV) & (!is.na(data$PSI)),]
  if (nrow(data) == 0){
    empty_ggplot_with_info(info = sprintf("No splicing evidence in\n %s,\n Try other tissues", tissue))
  }else{
    ggplot(data, aes(x=PSI, y=dbscSNV)) + 
      geom_point(aes(colour = Read >= min_read, shape = AF < 0.05), size = 2) + 
      scale_colour_manual(name = sprintf('Read count >= %s', min_read), values=setNames(c('#0066CC','#CCCCCC'),c(T, F))) +
      scale_shape_manual(name = 'Annotated AF < 0.05', values=setNames(c(16, 17), c(F, T))) +
      geom_rect(aes(xmin = min_psi-0.01, xmax = 1.01, ymin = min_score-0.01, ymax = 1.01), fill = "transparent", color = "red") +
      xlim(0,1.01) + ylim(0, 1.01) +
      theme(panel.background = element_rect(fill = "white", colour = "grey50"))
  }
}

draw_Spliceai = function(data, tissue, min_read, min_psi, min_score){
  psi_col = paste0("PSI_", tissue)
  rc_col = paste0("Read_", tissue)
  colnames(data)[colnames(data) == psi_col] = "PSI"
  colnames(data)[colnames(data) == rc_col] = "Read"

  data = data[!is.na(data$SpliceAI) & (!is.na(data$PSI)),]
  if (nrow(data) == 0){
    empty_ggplot_with_info(info = sprintf("No splicing evidence in\n %s,\n Try other tissues", tissue))
  }else{
    ggplot(data, aes(x=PSI, y=SpliceAI)) + 
      geom_point(aes(colour = Read >= min_read, shape = AF < 0.05), size = 2) + 
      scale_colour_manual(name = sprintf('Read count >= %s', min_read), values = setNames(c('#0066CC','#CCCCCC'),c(T, F))) +
      scale_shape_manual(name = 'Annotated AF < 0.05', values=setNames(c(16, 17), c(F, T))) +
      geom_rect(aes(xmin = min_psi-0.01, xmax = 1.01, ymin = min_score-0.01, ymax = 1.01), fill = "transparent",  color = "red") +
      xlim(0,1.01) + ylim(0, 1.01) +
      theme(panel.background = element_rect(fill = "white", colour = "grey50"))
  }
}


draw_gene_red_box = function(gene_name, gene_table, dir_path, tissue, rna_meta, min_read, min_psi, min_spliceai, min_dbscsnv){
  red_box_file = sprintf("%s%s%s_%s_sRV_candidate.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene_name)[1], gene_name)

  if (length(red_box_file) > 0){
    sqtl_candidate = read.table(red_box_file, header=T, sep="\t", stringsAsFactors = F)
    gene_red_box = list(SpliceAI = draw_Spliceai(sqtl_candidate, tissue, min_read, min_psi, min_spliceai), dbscSNV = draw_dbscSNV(sqtl_candidate, tissue, min_read, min_psi, min_dbscsnv)) 
  }else{
    sqtl_candidate = data.frame()
    gene_red_box = list(SpliceAI = empty_ggplot_with_info(sprintf("No sRV candidates")), dbscSNV = empty_ggplot_with_info(sprintf("No sRV candidates")))
  }
  
  intron_sample = count_cryptic_intron_gene(tissue, rna_meta, sqtl_candidate, gene.id = name2id(gene_table, gene_name)[1], min_psi, min_read)
  
  list(plots=gene_red_box, sqtl_candidate = sqtl_candidate, intron_sample = intron_sample)
}

select_sqtl_candidate = function(sqtl_candidate, tissue, predictor, bursh_psi_min, bursh_psi_max, bursh_score_min, bursh_score_max){
  psi_colname = paste0("PSI_", tissue)
  sqtl_candidate[is.na(sqtl_candidate)] = 0
  sqtl_candidate = sqtl_candidate[sqtl_candidate[[psi_colname]] >= bursh_psi_min & sqtl_candidate[[psi_colname]] <= bursh_psi_max & sqtl_candidate[[predictor]] >= bursh_score_min & sqtl_candidate[[predictor]] <= bursh_score_max, ]
  #                                c("SubjectID", "Junc", "Strand", "Event", "DNA_variant", "Spliceai", "SpliceAI_pred_match_junction", "dbscSNV", "AF", "Gene", "Gene_id", "Expand_region", "SpliceAI_pred_cryptic_exon")]
  #colnames(sqtl_candidate) = c("Patient", "Coordinates_of_novel_junction", "Strand", "Event", "DNA_variant", "SpliceAI", "SpliceAI_pred_match_junction", "dbscSNV", "AF", "gene_name", "gene_id", "Expand_region", "SpliceAI_pred_cryptic_exon")
  sqtl_candidate = separate(sqtl_candidate, col="Coordinates_of_novel_junc", sep=":-|:\\+", into = c("Coordinates_of_novel_junc", NA))
  sqtl_candidate
}


count_cryptic_intron_gene = function(tissue, rna_meta, sqtl_candidate, gene.id, min_psi, min_read){
  
  subject_cryptic_tbl = as.data.frame(subset(rna_meta, Tissue == tissue))
  rvjunc_pass = getsRVJunc(sqtl_candidate, psi_threshold = min_psi, rc_threshold = min_read, tissue, nr_tissues_threshold = 1, splice_score = "SpliceAI", splice_threshold = 0.2)
  rvjunc_pass = unique(sqtl_candidate[rvjunc_pass[["PASS"]], c("SubjectID", "Coordinates_of_novel_junc")])
  rvjunc_pass = rvjunc_pass %>% dplyr::count(SubjectID) 
  
  
  subject_cryptic_tbl$Nr_novel_junction_overlap_sRV = rvjunc_pass[match(subject_cryptic_tbl$SubjectID, rvjunc_pass$SubjectID), "n"]
  subject_cryptic_tbl$Nr_novel_junction_overlap_sRV[is.na(subject_cryptic_tbl$Nr_novel_junction_overlap_sRV)] = 0
  subject_cryptic_tbl = subject_cryptic_tbl[, c("SubjectID", "Group", "Nr_novel_junction_overlap_sRV")]
  subject_cryptic_tbl
}



