library(ggplot2)
library(stringr)

name2id = function(gene_tbl, name){
  gene_tbl[match(name, gene_tbl$gene.name), "gene.id"]
}

empty_ggplot_with_info = function(info){
  g = ggplot() + 
    annotate("text", x = 1, y = 1, size = 6, label = info) + theme_void()
  g
}


draw_Spliceai = function(data, tissue, min_read, min_psi, min_score){
  lcp_col = paste0("Outlier_P_", tissue)
  rc_col = paste0("Reads_", tissue)
  colnames(data)[colnames(data) == lcp_col] = "LeafCutter_Pval"
  colnames(data)[colnames(data) == rc_col] = "Reads"
  
  data = data[!is.na(data$SpliceAI) & (!is.na(data$LeafCutter_Pval)),]
  if (nrow(data) == 0){
    empty_ggplot_with_info(info = sprintf("No splicing evidence in\n %s,\n Try other tissues", tissue))
  }else{
    ggplot(data, aes(x=-log10(LeafCutter_Pval), y=SpliceAI)) + 
      geom_point(aes(colour = Reads >= min_read, shape = AF < 0.05), size = 2) + 
      scale_colour_manual(name = sprintf('Read count >= %s', min_read), values = setNames(c('#0066CC','#CCCCCC'),c(T, F))) +
      scale_shape_manual(name = 'Allele Freq < 0.05', values=setNames(c(16, 17), c(F, T))) +
      geom_rect(aes(xmin = -log10(0.05), xmax = Inf, ymin = min_score-0.01, ymax = 1.01), fill = "transparent",  color = "red") +
      ylim(0, 1.01) + expand_limits(x=0) +
      theme(panel.background = element_rect(fill = "white", colour = "grey50"))
  }
}


draw_gene_red_box = function(gene_name, gene_table, dir_path, tissue, rna_meta, min_read, min_psi, min_spliceai){
  red_box_file = sprintf("%s%s%s_%s_csQTL.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene_name)[1], gene_name)

  if (length(red_box_file) > 0){
    sqtl_candidate = read.table(red_box_file, header=T, sep="\t", stringsAsFactors = F)
    gene_red_box = list(SpliceAI = draw_Spliceai(sqtl_candidate, tissue, min_read, min_psi, min_spliceai)) 
  }else{
    sqtl_candidate = data.frame()
    gene_red_box = list(SpliceAI = empty_ggplot_with_info(sprintf("No csQTL candidates")))
  }
  
  # intron_sample = count_cryptic_intron_gene(tissue, rna_meta, sqtl_candidate, gene.id = name2id(gene_table, gene_name)[1], min_psi, min_read)
  
  list(plots=gene_red_box, sqtl_candidate = sqtl_candidate) # intron_sample = intron_sample
}

select_sqtl_candidate = function(sqtl_candidate, tissue, predictor, bursh_psi_min, bursh_psi_max, bursh_score_min, bursh_score_max){
  psi_colname = paste0("Outlier_P_", tissue)
  sqtl_candidate[[psi_colname]][is.na(sqtl_candidate[[psi_colname]])] = 1
  sqtl_candidate[is.na(sqtl_candidate)] = 0 ### 
  sqtl_candidate = sqtl_candidate[-log10(sqtl_candidate[[psi_colname]]) >= bursh_psi_min & -log10(sqtl_candidate[[psi_colname]]) <= bursh_psi_max & sqtl_candidate[[predictor]] >= bursh_score_min & sqtl_candidate[[predictor]] <= bursh_score_max, ]

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


compare_splice = function(gg_x, gg_y, sashimi2gg_tbl, counts, psi, rna_meta, subject, specific_junction=FALSE){

  if ( specific_junction==FALSE ){
    sashimi2gg_tbl = subset(sashimi2gg_tbl, start <= gg_x & end >= gg_x & (gg_y * ytext > 0))
    rownames(sashimi2gg_tbl) = paste(sashimi2gg_tbl$chr, sashimi2gg_tbl$startv, sashimi2gg_tbl$endv, sashimi2gg_tbl$strand, sep=":")
    sashimi2gg_tbl = unique(sashimi2gg_tbl[, c("startv", "endv")])

    plot_junction_list = rownames(sashimi2gg_tbl)
  }else{
    plot_junction_list = specific_junction
  }

  compare = list()
  c_idx = 1
  tissues = as.character(rna_meta[rna_meta$SubjectID == subject, "Tissue"])
  for (tissue in tissues){
    sample_record = subset(rna_meta, Tissue == tissue & SampleID %in% colnames(counts))
    sample_record$Group = str_wrap(sample_record$Group , width = 15)

    for (s in plot_junction_list){

      xlabs=paste(levels(as.factor(sample_record$Group)),"\n(N=",table(sample_record$Group),")", sep="")

      ### Compare read
      sample_record$read = unname(unlist(counts[s, as.character(sample_record$SampleID), drop=T]))
      compare[[c_idx]] = ggplot(data=sample_record, aes(x=Group, y=read)) +
                                      geom_boxplot() +
                                      ggtitle(sprintf("%s Junction: %s", subject, s)) + ylab(sprintf("Reads in %s samples", tissue)) +
                                      scale_x_discrete(labels=xlabs) +
                                      theme_classic() +
                                      theme(axis.title.x=element_blank())

      s_idx = match(subject, as.character(sample_record$SubjectID))
      if (!is.na(s_idx)){
        compare[[c_idx]] = compare[[c_idx]] + geom_point( data = data.frame(x=sample_record$Group[s_idx], y=sample_record$read[s_idx]),
                                                           aes(x=x, y=y),
                                                           color="red")
      }
      c_idx = c_idx + 1

      if (s %in% rownames(psi)){
        sample_record$psi = unname(unlist(psi[s, as.character(sample_record$SampleID), drop=T]))
        compare[[c_idx]] = ggplot(data=sample_record, aes(x=Group, y=-log10(psi))) +
          geom_boxplot() +
          ggtitle(sprintf("%s Junction: %s", subject, s)) + ylab(sprintf("-log10(LeafCutter P value) in %s samples", tissue)) +
          scale_x_discrete(labels=xlabs) +
          theme_classic() +
          theme(axis.title.x=element_blank())

        s_idx = match(subject, as.character(sample_record$SubjectID))
        if (!is.na(s_idx)){
          compare[[c_idx]] = compare[[c_idx]] + geom_point( data = data.frame(x=sample_record$Group[s_idx], y=-log10(sample_record$psi[s_idx])),
                                                            aes(x=x, y=y),
                                                            color="red")
        }
      }else{
        sample_record$psi = NA
        compare[[c_idx]] = empty_ggplot_with_info(sprintf("No LeafCutter P value for junction %s", s))
      }





      c_idx = c_idx + 1
    }
  }
  print(length(compare))
  compare
}


