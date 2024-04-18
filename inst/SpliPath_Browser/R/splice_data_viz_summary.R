library(ggplot2)
library(stringr)
library(ggVennDiagram)

name2id = function(gene_tbl, name){
  gene_tbl[match(name, gene_tbl$gene.name), "gene.id"]
}

empty_ggplot_with_info = function(info){
  g = ggplot() + 
      annotate("text", x = 1, y = 1, size = 6, label = info) + theme_void()
  g
}

plot_number_cryptic_intron = function(gene_splice, rna_meta, color_by, itissue, other_filters){
  
  for (x in other_filters[["add_filters"]]){
    type = class(rna_meta[[x]])
    if (type %in% c("numeric", "integer")) {
      min_value_idx = sprintf("%s_tab1_main_min", x)
      max_value_idx = sprintf("%s_tab1_main_max", x)
      rna_meta = rna_meta[(rna_meta[[x]] >= other_filters[[min_value_idx]]) & (rna_meta[[x]] <= other_filters[[max_value_idx]]), ]
      
    } else if (type %in% c("character", "logical", "factor", "Rle")) {
      value_idx = sprintf("%s_tab1_main", x)
      rna_meta = rna_meta[rna_meta[[x]] == other_filters[[value_idx]], ]
    }
  }
  
  if (nrow(rna_meta) == 0){
    summary_plot = empty_ggplot_with_info("All the samples were filtered out,\n please change filter thresholds.")
  }else{
    if (itissue != "All_Tissues"){
      rna_meta = subset(rna_meta, Tissue %in% itissue & SampleID %in% rownames(gene_splice))
      
      rna_meta = rna_meta[order(rna_meta$RIN, decreasing = T), ]
      rna_meta = rna_meta[!duplicated(rna_meta$SubjectID), ]
      rna_meta$novel_junction = gene_splice[rna_meta$SampleID, 1]
      
      if (nrow(rna_meta) == 0){
        summary_plot = empty_ggplot_with_info("All the samples were filtered out,\n please change filter thresholds.")
      }else{
        xlabs=paste(levels(as.factor(rna_meta[[color_by]])),"\n(N=",table(rna_meta[[color_by]]),")",sep="") 
        summary_plot = ggplot(rna_meta, aes_string(x = color_by, y = "novel_junction", fill = color_by, color = color_by)) +
        geom_violin(trim=TRUE, scale="count", alpha=0.5, color="white") +
        geom_boxplot(position = position_dodge(width = 0.9), width=0.1, fill="white") +
        ylab("Nr. novel junction") + 
        scale_x_discrete(labels=xlabs) +
        ggtitle(sprintf("Number of novel junctions in %s samples", nrow(rna_meta))) +
        scale_fill_manual(values = c( "#FFC61E", "#009ADE", "#AF58BA", "#F28522", "#00CD6C")) + #
        scale_color_manual(values = c( "#FFC61E", "#009ADE", "#AF58BA", "#F28522", "#00CD6C")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) 
      }
    }else{
      rna_meta = subset(rna_meta, SampleID %in% rownames(gene_splice))
      
      rna_meta = rna_meta[order(rna_meta$RIN, decreasing = T), ]
      rna_meta = rna_meta[!duplicated(rna_meta[, c("Tissue", "SubjectID")]), ]
      
      rna_meta$novel_junction = gene_splice[rna_meta$SampleID, ]
      if (nrow(rna_meta) == 0){
        summary_plot = empty_ggplot_with_info("All the samples were filtered out,\n please change filter thresholds.")
      }else{
        xlabs=paste(levels(as.factor(rna_meta[["Tissue"]])),"\n(N=",table(rna_meta[["Tissue"]]),")",sep="") 
        summary_plot = ggplot(rna_meta, aes_string(x = "Tissue", y = "novel_junction", color = color_by, fill=color_by)) +
          geom_violin(trim=TRUE, scale="count") +
          ylab("Nr. novel junction") + 
          scale_x_discrete(labels=xlabs) +
          ggtitle(sprintf("Number of novel junctions in %s samples", nrow(rna_meta))) +
          scale_fill_manual(values = c( "#FFC61E", "#009ADE", "#AF58BA", "#F28522", "#00CD6C")) + #
          scale_color_manual(values = c( "#FFC61E", "#009ADE", "#AF58BA", "#F28522", "#00CD6C")) +
          theme(panel.background = element_rect(fill = "white", colour = "grey50")) 
      }
    }
  }
  list(summary_plot = summary_plot, meta_table = rna_meta)
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
      
      ### Compare PSI
      c_idx = c_idx + 1
      # if (s %in% rownames(psi)){
      #   sample_record$psi = unname(unlist(psi[s, as.character(sample_record$SampleID), drop=T]))
      #   compare[[c_idx]] = ggplot(data=sample_record, aes(x=Group, y=psi)) +
      #     geom_boxplot() +
      #     ggtitle(sprintf("%s Junction: %s", subject, s)) + ylab(sprintf("PSI in %s samples", tissue)) +
      #     scale_x_discrete(labels=xlabs) +
      #     theme_classic() +
      #     theme(axis.title.x=element_blank())
      #   
      #   s_idx = match(subject, as.character(sample_record$SubjectID))
      #   if (!is.na(s_idx)){
      #     compare[[c_idx]] = compare[[c_idx]] + geom_point( data = data.frame(x=sample_record$Group[s_idx], y=sample_record$psi[s_idx]), 
      #                                                                     aes(x=x, y=y), 
      #                                                                     color="red")
      #   }
      # }else{
      #   sample_record$psi = NA
      #   compare[[c_idx]] = empty_ggplot_with_info(sprintf("PSI not recorded for annotated junction %s", s))
      # }
      
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
  
  compare
}


# getsRVJunc_v1 = function(splice_var_summary, psi_threshold, rc_threshold, tissues, nr_tissues_threshold, splice_score, splice_threshold){
#   
#   psi_column = paste0("PSI_", tissues)
#   read_column = paste0("Read_", tissues)
# 
#   psi_pass1 = splice_var_summary[, psi_column] >= psi_threshold
#   rc_pass1 = splice_var_summary[, read_column] >= rc_threshold
#   #rc_pass2 = splice_var_summary[, read_column] >= 20
#   
#   splice_pass = (splice_var_summary[[splice_score]] >= splice_threshold)
#   
#   PASS_per_tissue = (psi_pass1 * rc_pass1) #|rc_pass1
#   PASS_per_tissue = PASS_per_tissue & splice_pass 
#   PASS_per_tissue[is.na(PASS_per_tissue)] = FALSE
#   #PASS = rowSums(PASS_per_tissue, na.rm = T) >= nr_tissues_threshold
#   
#   list(PASS = PASS, PASS_per_tissue = PASS_per_tissue)
# }

getsRVJunc = function(splice_var_summary, psi_threshold, rc_threshold, tissues, nr_tissues_threshold, splice_score, splice_threshold){
  
  read_column = paste0("Read_", tissues)
  
  rc_pass1 = splice_var_summary[, read_column, drop=F] >= rc_threshold

  splice_pass = (splice_var_summary[[splice_score]] >= splice_threshold)
  
  PASS_per_tissue = rc_pass1 & splice_pass 
  PASS_per_tissue[is.na(PASS_per_tissue)] = FALSE

  list(PASS = rowSums(PASS_per_tissue, na.rm = T) >= nr_tissues_threshold, PASS_per_tissue = PASS_per_tissue)
}


getNovelJuncPassThres = function(psi_tbl, count_tbl, psi_threshold = 0.2, count_threshold=5){
  psi_tbl[psi_tbl < psi_threshold] = 0
  psi_tbl[psi_tbl >= psi_threshold] = 1
  
  count_tbl[count_tbl < count_threshold] = 0
  count_tbl[count_tbl >= count_threshold] = 1
  
  psi_tbl * count_tbl
}
  
gene_tissue_nr_junc_venn = function(gene.list, tissues, subject.list, dir_path, gene_table, rna_meta, color_by = "Group", 
                                    psi_threshold = 0.2, count_threshold=5, nr_tissues_threshold =1, 
                                    junc_type = "ur-sQTL novel junctions", splice_score = "SpliceAI", splice_score_thres = 0.2){

  rna_meta = rna_meta[!duplicated(rna_meta[, c("Tissue", "SubjectID")]), ]
  rna_meta = subset(rna_meta, Tissue %in% tissues )
  
  ### If no valid subject id input, return plots with error message
  if (sum(subject.list %in% as.character(rna_meta$SubjectID)) == 0){
    
    sample_junc_plot = empty_ggplot_with_info('No valid input SubjectID,\nplease use SubjectIDs in "Meta data" in "Sample Overview" panel')
    srv_tbl_genes = data.frame(info = 'No valid input SubjectID,\nplease use SubjectIDs in "Meta data" in "Sample Overview" panel')
    venn_plot = empty_ggplot_with_info('No valid input SubjectID,\nplease use SubjectIDs in "Meta data" in "Sample Overview" panel')
    
  }else if(sum(gene.list %in% upload_file$view_gene_list) == 0){
    
    sample_junc_plot = empty_ggplot_with_info('No valid input genes')
    srv_tbl_genes = data.frame(info = 'No valid input genes')
    venn_plot = empty_ggplot_with_info('No valid input genes')
    
  }else{
    
    if (subject.list != ""){
      rna_meta = subset(rna_meta, SubjectID %in% subject.list )
    }
    sample.list = as.character(rna_meta$SampleID)
    
    srv_tbl_genes = c()
    if (junc_type == "ur-sQTL novel junctions"){
      gene.list_ = c()
      for (gene.name in gene.list){
        # if (file.exists(sprintf("%s%s%s_%s_sRV_candidate.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene.name)[1], gene.name))){
        #   srv_file = sprintf("%s%s%s_%s_sRV_candidate.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene.name)[1], gene.name)
        if (file.exists(sprintf("%s%s%s_%s_crossref.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene.name)[1], gene.name))){
          srv_file = sprintf("%s%s%s_%s_crossref.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene.name)[1], gene.name)
          srv_tbl_genes = rbind(srv_tbl_genes, read.table(srv_file, header=T, sep="\t", stringsAsFactors = F))
          gene.list_ = c(gene.list_, gene.name)
        }else{print(sprintf("No sRV info for %s", gene.name))}
      }
      if (subject.list != ""){
        srv_tbl_genes = subset(srv_tbl_genes, SubjectID %in% subject.list)
      }
      if (nrow(srv_tbl_genes) == 0){
        srv_tbl_genes[1:(length(gene.list_)*length(subject.list)), ] = 0
        srv_tbl_genes$Gene = rep(gene.list_, each = length(subject.list))
        srv_tbl_genes$SubjectID = rep(subject.list, time = length(gene.list_))
      }
      filters = getsRVJunc(srv_tbl_genes, psi_threshold, count_threshold, tissues, nr_tissues_threshold, splice_score, splice_score_thres)
      srv_tbl_genes[, tissues] = filters$PASS_per_tissue
      srv_tbl_genes = srv_tbl_genes[filters$PASS, c("Gene", "Coordinates_of_novel_junc", "SubjectID", tissues)]
      
      # Genes habor novel junctions in each tissues
      venn_data = lapply(tissues, FUN = function(tissues_pass){
        unique(srv_tbl_genes[["Gene"]][srv_tbl_genes[[tissues_pass]]])
      })
      names(venn_data) = tissues
      
      # Nr of novel junctions per sample
      sample_junc_tbl = data.frame(Gene = rep(srv_tbl_genes$Gene, times=length(tissues)),
                                   Coordinates_of_novel_junc = rep(srv_tbl_genes$Coordinates_of_novel_junc, times=length(tissues)),
                                   SubjectID =  rep(srv_tbl_genes$SubjectID, times=length(tissues)),
                                   Tissue = rep(tissues, each = nrow(srv_tbl_genes)),
                                   novel_junction = unname(unlist(srv_tbl_genes[, tissues])), stringsAsFactors = F)
      sample_junc_tbl = sample_junc_tbl %>% group_by(SubjectID, Tissue, Coordinates_of_novel_junc) %>% summarise(novel_junction = sum(novel_junction) > 0)
      sample_junc_tbl = sample_junc_tbl %>% group_by(SubjectID, Tissue) %>% summarise(across(novel_junction, sum))
      rna_meta = left_join(rna_meta, sample_junc_tbl, by=c("SubjectID", "Tissue"))
      rna_meta[is.na(rna_meta$novel_junction), "novel_junction"] = 0
      
      # Table: Gene in row; Tissue in col; elements are nr of novel junction
      srv_tbl_genes = srv_tbl_genes %>% group_by(Gene) %>% summarise(across(one_of(tissues), sum))
    }else{
      
      psi_tbl_genes = c()
      count_tbl_genes = c()
      for (gene.name in gene.list){
        if (file.exists(sprintf("%s%s%s_%s_intron_psi.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene.name)[1], gene.name))){
          psi_file = sprintf("%s%s%s_%s_intron_psi.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene.name)[1], gene.name)
          psi_tbl_genes = rbind(psi_tbl_genes, read.table(psi_file, header=T, sep="\t", stringsAsFactors = F, check.names = F))
          count_file = sprintf("%s%s%s_%s_intron_count.txt.gz", dir_path, .Platform$file.sep, name2id(gene_table, gene.name)[1], gene.name)
          count_tbl_genes = rbind(count_tbl_genes, read.table(count_file, header=T, sep="\t", stringsAsFactors = F, check.names = F))
        }else{print(sprintf("No splicing junction info for %s", gene.name))}
      }
      
      psi_tbl_genes = unique(psi_tbl_genes)
      rownames(psi_tbl_genes) = paste(psi_tbl_genes$chr, psi_tbl_genes$start, psi_tbl_genes$end, psi_tbl_genes$strand, sep=":")
      psi_tbl_genes_name = psi_tbl_genes$gene.name
      psi_tbl_genes = psi_tbl_genes[, colnames(psi_tbl_genes) %in% sample.list]
      
      count_tbl_genes = unique(count_tbl_genes)
      rownames(count_tbl_genes) = paste(count_tbl_genes$chr, count_tbl_genes$start, count_tbl_genes$end, count_tbl_genes$strand, sep=":")
      count_tbl_genes = count_tbl_genes[rownames(psi_tbl_genes), colnames(psi_tbl_genes)]
  
      srv_tbl_genes = getNovelJuncPassThres(psi_tbl_genes, count_tbl_genes, psi_threshold, count_threshold)
      
      # Nr of novel junctions per sample
      rna_meta$novel_junction = colSums(srv_tbl_genes[, as.character(rna_meta$SampleID)])
  
      # Flatten matrix
      srv_tbl_genes = data.frame(Gene = rep(psi_tbl_genes_name, times = ncol(srv_tbl_genes)),
                                 SampleID = rep(colnames(srv_tbl_genes), each = nrow(srv_tbl_genes)),
                                  PASS = unname(unlist(srv_tbl_genes)), stringsAsFactors = F)
      srv_tbl_genes$Tissue = rna_meta[match(srv_tbl_genes$SampleID, as.character(rna_meta$SampleID)), "Tissue"]
      
      # Genes habor novel junctions in each tissues
      venn_data = lapply(tissues, FUN = function(tissues_pass){
        unique(srv_tbl_genes[(srv_tbl_genes$Tissue == tissues_pass) & (srv_tbl_genes$PASS > 0), "Gene"])
      })
      names(venn_data) = tissues
  
      # Table: Gene in row; Tissue in col; elements are nr of novel junction
      srv_tbl_genes = srv_tbl_genes %>% group_by(Gene, Tissue) %>% summarise(across(PASS, sum))
      srv_tbl_genes = spread(srv_tbl_genes, key = Tissue, value = PASS)
    }
    
    gene_no_junc_ = gene.list[!gene.list %in% srv_tbl_genes$Gene]
    if (length(gene_no_junc_) > 0){
      gene_no_junc = data.frame(matrix(0, nrow=length(gene_no_junc_), ncol=ncol(srv_tbl_genes)))
      colnames(gene_no_junc) = colnames(srv_tbl_genes)
      gene_no_junc$Gene = gene_no_junc_
      print(colnames(gene_no_junc) %in% tissues)
      gene_no_junc[, colnames(gene_no_junc) %in% tissues] = 0
      srv_tbl_genes = rbind(srv_tbl_genes, gene_no_junc)
    }
    
    rna_meta$Tissue = factor(rna_meta$Tissue, levels = tissues)
    xlabs=paste(levels(rna_meta[["Tissue"]]),"\n(N=", table(rna_meta[["Tissue"]]),")",sep="") 
    sample_junc_plot = ggplot(rna_meta, aes_string(x = "Tissue", y = "novel_junction", color = color_by, fill=color_by)) +
      geom_violin(trim=TRUE, scale="count", alpha=0.5, color="white") +
      geom_boxplot(position = position_dodge(width = 0.9), width=0.1, fill="white") +
      ylab(sprintf("Nr. %s", junc_type)) + 
      scale_x_discrete(labels=xlabs) +
      ggtitle(sprintf("%s in %s query genes in %s individuals", junc_type, length(gene.list), length(unique(as.character(rna_meta$SubjectID))))) +
      scale_fill_manual(values = c( "#FFC61E", "#009ADE", "#AF58BA", "#F28522", "#00CD6C")) + #
      scale_color_manual(values = c( "#FFC61E", "#009ADE", "#AF58BA", "#F28522", "#00CD6C")) +
      theme(panel.background = element_rect(fill = "white", colour = "grey50")) 
    
    # venn_plot = ggvenn(venn_data, show_percentage = F, stroke_size = 0.5, aes(fill=count)) +
    #   ggtitle("Number of genes expressed novel junctions in each tissue")
      
    venn_plot = ggVennDiagram(venn_data, label_alpha = 0, label = "count") +
                  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")  +
                  guides(fill=guide_legend(title="Number of genes")) +
                  ggtitle("Number of genes expressed novel junctions in each tissue")
    }
  list(sample_junc_plot = sample_junc_plot, gene_tissue_tbl = data.frame(srv_tbl_genes), gene_tissue_venn_plot = venn_plot)
}




