#' Find the junctions that share a same annotated splice sites with the query novel junctions
#'
#' @param junc_file string or data.frame The output file of mergeJuncData function, which contain junction-sample read counts table. The first four columns should be the coordinates of each junction: "chr", "start", "end", "strand".
#' @param junc_anno string or data.frame The output file of annotateJunc function, which contain junction annotations.
#' @param annotation vector Annotation sources (should in colnames in junc_anno) that were combined to determine novel junction. A junction not annotated in any given source is treated as novel junction. Default: c("in.ensembl", "in.snaptron")
#' @param output_prefix string Prefix of output PSI file. The junction-sample PSI table will be written in a prefix_norm.txt.gz file.
#'
#' @return data.frame Seven-column dataframe. Contain coordinates, mapped gene, and splicing event of each novel junction.
#' @export
findBackGroundJunc <- 
  function(junc_file, junc_anno, output_prefix, annotation = c("in.ensembl", "in.snaptron"), return_value = T){
    if (is.character(junc_file)){
      junc_file = read.table(junc_file, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
    }
    if (is.character(junc_anno)){
      junc_anno = read.table(junc_anno, sep = "\t", header = T, stringsAsFactors = F)
    }
    samples = colnames(junc_file)[5:ncol(junc_file)]
    junc_anno_col = colnames(junc_anno)
    
    ### Add annotation columns to junction - sample dataframe
    junc_file$idx = paste(junc_file[["chr"]], junc_file[["start"]], junc_file[["end"]], junc_file[["strand"]], sep=":")
    junc_anno$idx = paste(junc_anno[["chr"]], junc_anno[["start"]], junc_anno[["end"]], junc_anno[["strand"]], sep=":")
    junc_file = junc_file[junc_file[['idx']] %in% junc_anno[['idx']], ]
    junc_file = subset(junc_file, select = -idx)
    junc_anno = subset(junc_anno, select = -idx)
    junc_file = dplyr::inner_join(junc_file, junc_anno, by=c('chr', 'start', 'end', 'strand'), multiple = "all")
    junc_file[, annotation] = sapply(junc_file[, annotation], FUN = as.logical)
    junc_file$is.anno = rowSums(junc_file[, annotation], na.rm = T) > 0
    rm(junc_anno)
    
    ### Find the corresponding annotated junction that share same sites with the novel junctions
    rownames(junc_file) = paste(junc_file[["chr"]], junc_file[["start"]], junc_file[["end"]], junc_file[["strand"]], junc_file[["gene.id"]], sep=":")
    junc_info = junc_file[, c('chr', 'start', 'end', 'strand', 'gene.id', 'gene.name', annotation, 'is.anno', 'variant', junc_anno_col)] #
    junc_info = junc_info[junc_info$variant != "alternative_intron", ] 
    
    junc_info$bg.junc1 = NA
    junc_info$bg.junc2 = NA
    
    genes = unique(junc_info[(junc_info$is.anno == FALSE) , 'gene.id'])
    for (gene in genes){
      #gene_junc = annoJuncWithShareSites(junc_info, gene) ### release annoJuncWithShareSite function
      ### START
      gene_junc = junc_info[(junc_info$gene.id == gene), ]
      if (length(unique(gene_junc$is.anno)) > 1){
        novel_junc = rownames(gene_junc[gene_junc$is.anno == FALSE, ])
        for (i in novel_junc){
          if ((gene_junc[i, 'variant'] == 'exon_skipping') | (gene_junc[i, 'variant'] == 'alternative_intron')){
            gene_junc[i, 'bg.junc1'] = annoJuncWithShareSites(gene_junc, i, 'start')
            gene_junc[i, 'bg.junc2'] = annoJuncWithShareSites(gene_junc, i, 'end')
          }else if(gene_junc[i, 'start.anno'] == FALSE){
            gene_junc[i, 'bg.junc1'] = annoJuncWithShareSites(gene_junc, i, 'end')
          }else if(gene_junc[i, 'end.anno'] == FALSE){
            gene_junc[i, 'bg.junc1'] = annoJuncWithShareSites(gene_junc, i, 'start')
          }
        }
      }
      ### END
      junc_info[rownames(gene_junc), c("bg.junc1", "bg.junc2")] = gene_junc[, c("bg.junc1", "bg.junc2")]
    }
    junc_info = junc_info[(junc_info[['is.anno']] == FALSE), ]
    junc_info = junc_info[!is.na(junc_info[['bg.junc1']]), ]
    
    ### Sum the reads of the corresponding annotation junctions
    bg_ls = unique(c(na.omit(junc_info[["bg.junc1"]]), na.omit(junc_info[["bg.junc2"]])))
    bg_2_read = list()
    bg_2_read = lapply(bg_ls, function(bg){
      bg_vec = unname(unlist(strsplit(bg, split="[|]")))
      unname(colSums(junc_file[bg_vec, samples]))
    })
    names(bg_2_read) = bg_ls
    rm(junc_file)
    
    ### Output the coordinates of the background junctions that mapped to each novel junctions
    output_gz = gzfile(sprintf("%s_background_info.txt.gz", output_prefix), "w")
    write.table(format(junc_info[, c('chr', 'start', 'end', 'strand', 'gene.id', 'gene.name', annotation, 'variant', "bg.junc1", "bg.junc2")], scientific=FALSE), output_gz, row.names = F, col.names = T, sep = "\t", quote = F)
    close(output_gz)
    
    junc_info[samples] = NA
    for (novel_idx in rownames(junc_info)){
      
      bg.junc1 = junc_info[novel_idx, "bg.junc1"]
      bg.junc2 = junc_info[novel_idx, "bg.junc2"]
      
      if (junc_info[novel_idx, "variant"] %in% c("exon_skipping", "alternative_intron")){
        
        if (!is.na(bg.junc1)){
          bg.junc1.read = bg_2_read[[bg.junc1]]
        }else{bg.junc1.read = rep(NA, length(samples))}
        
        if (!is.na(bg.junc2)){
          bg.junc2.read =  bg_2_read[[bg.junc2]]
        }else{bg.junc2.read = rep(NA, length(samples))}
        
        junc_info[novel_idx, samples] = pmax(bg.junc1.read, bg.junc2.read, na.rm = TRUE)
        
      }else{ ### novel donor or acceptor
        if (!is.na(bg.junc1)){
          junc_info[novel_idx, samples] = bg_2_read[[bg.junc1]]
        }
      }
      
    }
    junc_info = junc_info[, c('chr', 'start', 'end', 'strand', 'gene.id', 'gene.name', annotation, 'variant', "bg.junc1", "bg.junc2", samples)]
    
    output_gz = gzfile(sprintf("%s_background_junc.txt.gz", output_prefix), "w")
    write.table(format(junc_info, scientific=FALSE), output_gz, row.names = F, col.names = T, sep = "\t", quote = F)
    close(output_gz)
    
    if (return_value){
      junc_info
    }
  }
