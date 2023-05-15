library(foreach)
library(ggplot2)
library(rvat)
library(RSQLite)
library(SummarizedExperiment)

length_transform <- function(g){ g/100 }

make_variant_plot = function(gdb_path, subject, gene.name, chr, start, end, wgs_id, wgs_meta, expand_width, plot_xlim){
  chrom = paste0("chr", chr)
  
  ### Get variant Info
  gdb=rvat::gdb(gdb_path)
  spliceaiTerm = c("spliceaiDS_AG", "spliceaiDS_AL", "spliceaiDS_DG", "spliceaiDS_DL")
  dbscSNVTerm =  c("ada_score", "rf_score")
  
  vars <- RSQLite::dbGetQuery(gdb, 
                              statement = sprintf("select VAR_id, var.POS, var.CHROM, var.REF, var.ALT,
                              SpliceAI.spliceaiSYMBOL, SpliceAI.%s, dbscSNV.ada_score, dbscSNV.rf_score
                              from var
                              left join SpliceAI using ('VAR_id')
                              left join dbscSNV using ('VAR_id')
                              where var.CHROM in ('%s') and var.POS between %s and %s
                              ", paste(spliceaiTerm, collapse=", SpliceAI."), paste0(paste0("chr", chr), "' , '", sub("^chr", "", chr), "' , '", chr), start - expand_width, end + expand_width))
  
  rownames(vars) = vars$VAR_id
  vars$DNA_variant = paste0(vars$CHROM, ":g.", vars$POS, vars$REF, ">", vars$ALT)
  vars[,c(spliceaiTerm, dbscSNVTerm)] = sapply(vars[,c(spliceaiTerm, dbscSNVTerm)], FUN=as.numeric)
  vars[is.na(vars)] = 0

  vars$SpliceAI = do.call(pmax, c(vars[, spliceaiTerm], na.rm=T))
  vars$dbscSNV = do.call(pmax, c(vars[, dbscSNVTerm], na.rm=T))
  GT = rvat::getGT(gdb, VAR_id=vars$VAR_id, cohort="pheno")
  AF = data.frame(rvat::getAF(GT))
  colnames(AF) = "AF"
  vars$AF = AF[match(as.character(vars$VAR_id), rownames(AF)), "AF"]
  
  ### Get genotype of the query subject
  SM=RSQLite::dbGetQuery(gdb,sprintf("select * from %s", "pheno"))
  wgs_id = wgs_id[wgs_id %in% SM$IID]
  
  geno_plot = list()
  if (length(wgs_id) == 0){
    geno_plot[[1]] = empty_ggplot_with_info(info = sprintf("There is no WGS data for subject %s", subject))
  }else{
    for (i in 1:length(wgs_id)){
      wgs_id_ = wgs_id[i]
      vars$geno = assays(GT)$GT[as.character(vars$VAR_id), wgs_id_]
      vars = subset(vars, (!is.na(geno) & geno != 0))
      
      geno_theme_empty <- theme_bw(base_size = 11 )
      geno_theme_empty$panel.background = element_rect(fill="white", colour = "white")
      geno_theme_empty$axis.line = element_line()
      geno_theme_empty$line = element_blank()
      geno_theme_empty$rect <- element_blank()
      geno_theme_empty$axis.title.x=element_blank()
      geno_theme_empty$axis.text=element_blank()
      geno_theme_empty$axis.ticks.x=element_blank()
      geno_theme_empty$axis.line.y = element_line(arrow = grid::arrow(length = unit(0.2, "cm"), ends = "last"))
      geno_theme_empty$legend.position="bottom"
      geno_theme_empty$legend.justification = 'right'
      geno_theme_empty$legend.text=element_text(size=10)
      
      if (nrow(vars) > 0){
        vars$dist = vars$POS - start
        vars$xaxis = length_transform(vars$dist)
        
        geno_plot[[i]] = ggplot(vars, aes(x=xaxis, y=SpliceAI, color=substr(as.character(dbscSNV>=0.7),1,1), label = DNA_variant)) +
          geom_hline(aes(yintercept=0.5, linetype = "0.5"), color = "grey", size=1) +
          geom_hline(aes(yintercept=0.2, linetype = "0.2"), color = "grey", size=1) +
          scale_linetype_manual(name = "SpliceAI threshold", values = c(3, 1), guide = guide_legend(override.aes = list(linetype = c(3, 1)))) +
          geom_point(aes(shape=substr(as.character(AF<0.05),1,1)), size=2) + 
          geom_text(aes(label=ifelse(AF<0.05 & (SpliceAI >= 0.2 | dbscSNV>=0.7) , DNA_variant, '')), vjust=-0.5) + 
          scale_shape_manual(name = 'Allele frequency in GDB < 0.05', values=setNames(c(16, 17), c("F", "T"))) +
          scale_colour_manual(name = 'dbscSNV>=0.7', values = setNames(c('#0066CC', 'black'),c("T", "F"))) +
          xlim(plot_xlim) + ylim(c(0, 1)) +
          geno_theme_empty
      }else{
        geno_plot[[i]] = ggplot() +
          geom_hline(aes(yintercept=0.5, linetype = "0.5"), color = "grey", size=1) +
          geom_hline(aes(yintercept=0.2, linetype = "0.2"), color = "grey", size=1) +
          scale_linetype_manual(name = "SpliceAI threshold", values = c(3, 1), guide = guide_legend(override.aes = list(linetype = c(3, 1)))) +
          xlim(plot_xlim) + ylim(c(0, 1)) +
          geno_theme_empty
      }
      if ("Sample_Source" %in% colnames(wgs_meta)){
        geno_plot[[i]] = geno_plot[[i]] + 
                         ggtitle(sprintf("%s DNA variants - from %s sample", subject, as.character(wgs_meta[match(wgs_id_, as.character(wgs_meta$SampleID)), "Sample_Source"])))
      }else{
        geno_plot[[i]] = geno_plot[[i]] + ggtitle(sprintf("%s DNA variants", subject))
      }
    }
  }
  if (length(geno_plot) > 1){
    for (i in 1:(length(geno_plot)-1)){
      geno_plot[[i]] = geno_plot[[i]] + guides(color=FALSE, shape=FALSE, linetype=FALSE)
    }
  }
  geno_plot
}

make_sashimi_like_plot <- function(
    subject,
    exons_table = NULL,
    meta = NULL,
    counts = NULL, 
    intron_meta = NULL,
    snp_pos=NA){
  
  meta$group=as.character(meta$Tissue)

  y <- t(counts)
  x <- meta$group
  groups=sort(unique(x))
  
  length_transform <- function(g){ g/100 }
  summary_func <- colSums
  
  max_log=.5*ceiling(2*log10( 1+max( unlist( foreach (tis=groups) %do% { intron_meta$counts=summary_func(y[ tis==x,,drop=F]) } ) ) ))
  breaks=if (max_log <= 2.5) seq(0,max_log,by=0.5) else seq(0,ceiling(max_log),by=1)
  limits=c(0.0,max_log)
  
  ### Get unique intron start and end for docking
  intron_meta$id=1:nrow(intron_meta) # number each junction
  temp=intron_meta[,c("id","start","end")]
  m=reshape2::melt(temp, id.vars = "id") # melt to get list of all coordinates
  s=unique(m$value) # get unique start and end values
  s=sort(s)
  d=s[2:length(s)]-s[1:length(s)-1] # get the difference between pairs of coordinates
  trans_d <- length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33 # apply a trasnformation function - doesn't work on negative numbers!
  coords <- c(0,cumsum(trans_d))
  names(coords)=s
  
  total_length=sum(trans_d) # ==max(coords)
  my_xlim=c(-.02*total_length,1.02*total_length)
  expand_width = round(total_length * 100 * 0.02)
  #my_xlim=c(0, total_length)
  first_plot=T

  ##############
  # PLOT SETTINGS
  ###############
  
  new_theme_empty <- theme_bw(base_size = 11 )
  new_theme_empty$panel.background = element_rect(fill="white", colour = "white")
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
  main_title = paste0(subject, " - Chr", intron_meta$chr[1], ":", min(s)-expand_width, "-", max(s)+expand_width)
  #main_title = paste0("Chr", intron_meta$chr[1], ":", min(s), "-", max(s))
  legend_title <- "Read"
  junction_colour <- "red"
  cryptic_colour <- "grey"
  mainPalette <- c(junction_colour, cryptic_colour)
  names(mainPalette) = c("Annotated", "Novel")
  
  min_height=0
  max_height=0
  curv = 0.15
  min_exon_length <- 0.5
  maxcount=0
  mincount=1.0
  yFactor = 0.65   # originally set as 0.65
  yConstant = -0.25 # originally set as 0.5
  labelTextSize=4 # orignally set as 5
  curveMax = 10
  curveExponent = 2
  yOffset = 0
  centreLineWidth = 3 
  
  #summary_func=function(a) apply( sweep(a,1,rowSums(a),"/"),2, function(g) mean(g, na.rm=T) ) # Cal PSI in leafcuter
  maxcount=max(c(max(y), maxcount)) # find max and min values of summarised counts
  mincount=min(c(min(y), mincount))
  last_group=groups[length(groups)]
  
  plots <- list()
  for( fancyVar in 1:length(groups) ){
    
    intron_meta$counts=y[groups[fancyVar]==x,,drop=T]
    intron_meta$prop=intron_meta$counts/sum(intron_meta$counts)  # this has been changed

    # for each junction
    
    allEdges=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      if (i%%2==1) return(NULL)  # only care about the even numbered junctions?
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      l=end-start  # length of junction
      edge = data.frame(start, end)
      edge$chr = intron_meta$chr[i]
      edge$strand = intron_meta$strand[i]
      edge$startv <- intron_meta$start[i]
      edge$endv <- intron_meta$end[i]
      edge$log10counts=intron_meta$prop[i]+1
      
      # label each junction
      edge$label=intron_meta$counts[i]
      edge$clu<- intron_meta$clu[i]
      edge$Group <- i
      edge$xtext <-start+l/2
      edge$ytext <- -( (l^(yFactor) / 2) + yConstant)  # magic formula here
      edge$verdict <- ifelse( intron_meta$verdict[i] == "Annotated", yes = "Annotated", no ="Novel")
      edge
    })
    

    allEdgesP=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      if (i%%2==0) return(NULL)  # just the odd numbered junctions
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      l=end-start

      edge = data.frame(start, end)
      edge$chr = intron_meta$chr[i]
      edge$strand = intron_meta$strand[i]
      edge$startv <- intron_meta$start[i]
      edge$endv <- intron_meta$end[i]
      edge$log10counts=intron_meta$prop[i]+1
      edge$label=intron_meta$counts[i]
      #edge$clu<- intron_meta$clu[i]
      edge$Group <- i
      edge$xtext <-start+l/2
      edge$ytext <- (l^(yFactor) / 2)  + yConstant
      #edge$SIZE <- intron_meta$prop[i]+1
      edge$verdict <- ifelse( intron_meta$verdict[i] == "Annotated", yes = "Annotated", no ="Novel")
      edge
    })
    
    allEdges2gg = rbind(allEdgesP, allEdges)
    YLIMP = 0.5
    YLIMN = -0.5
    
    if (!is.null(allEdgesP)){
      allEdgesP = subset(allEdgesP, label > 0)
      if (any(nrow(allEdgesP))){
        YLIMP <- max(allEdgesP$ytext) + 0.25 * max(allEdgesP$ytext)
      }
    }
    if (!is.null(allEdges)){
      allEdges = subset(allEdges, label > 0)
      if (any(nrow(allEdges))){
          YLIMN <- min(allEdges$ytext) + 0.25 * min( allEdges$ytext)
      }
    }
    
    if ( all(is.na(main_title)) | !first_plot){
      new_theme_empty$plot.title <- element_blank()
    }
    first_plot <- FALSE
    
    g <- ggplot() 
    if (is.data.frame(allEdgesP) & any(nrow(allEdgesP))){
      g = g + 
        geom_curve(data=allEdgesP, aes(x = start, xend = xtext, y = yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                              angle=90, curvature=-curv,lineend="round") +
        geom_curve(data=allEdgesP, aes(x = xtext, xend = end, y = ytext, yend = yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                   angle=90, curvature=-curv,lineend="round") +
        geom_label(data=allEdgesP,aes(x=xtext,y=0.95*ytext,label=label), size = labelTextSize, label.size = NA, parse=TRUE, fill = "white",colour = "black", label.r = unit(0.3,"lines"), label.padding = unit(0.3,"lines") )
    }
    
    if (is.data.frame(allEdges) & any(nrow(allEdges))){
      g = g + 
        geom_curve(data=allEdges, aes(x = start, xend = xtext, y = -yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                   angle=90,curvature=curv,lineend="round") +
        geom_curve(data=allEdges, aes(x = xtext, xend = end, y = ytext, yend = -yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                   angle=90,curvature=curv,lineend="round") +
        geom_label(data=allEdges,aes(x=xtext,y=0.95*ytext,label=label), size= labelTextSize, label.size = NA, parse=TRUE, fill = "white", colour = "black", label.r = unit(0.3,"lines"), label.padding = unit(0.3,"lines") ) 
        
    }
    g = g + 
      new_theme_empty +
      # make the y axis label the group
      ylab(paste0(groups[fancyVar])) +
      xlab("") +
      xlim(my_xlim) +
      ggtitle(paste0(groups[fancyVar]) ) +
      
      # horizontal line - smooth out the ends of the curves
      geom_hline(yintercept=0, size = centreLineWidth, colour = "white") +
      geom_hline(yintercept=0,alpha=.9, size=1, colour="black") +
      
      ylim(YLIMN, YLIMP) +
      scale_size_continuous(limits=c(0,10),guide='none')
    
    plots[[fancyVar]] <- g  # return the plot
  }
    
  # ADDING EXON ANNOTATION
  df <- data.frame(x=coords, xend=total_length*(s-min(s))/(max(s)-min(s)), y=0, yend=min_height) # coords: the position of splice sites in the plot
  exons_chr <- exons_table[ exons_table$gene.name==intron_meta$gene.name[1], ] # subset by chr
  stopifnot( nrow(exons_chr) > 0)
  
  exons_here <- exons_chr[ ( min(s) <= exons_chr$start & exons_chr$start <= (max(s) + 1)) | ( min(s) <= exons_chr$end & exons_chr$end <= (max(s) + 1) ), ] # find exons
  
  # if exons found - remove exons that don't actually start or end with a junction
  # and any repeated exons or any exons larger than 500bp
  if( nrow(exons_here) > 0 ){
    exons_here <-  unique(
      exons_here[ ( exons_here$end %in% intron_meta$start | exons_here$start %in% intron_meta$end | exons_here$start %in% (intron_meta$end + 1) ), ] #
    )
  }
  
  # if any exons survive the cull
  if ( nrow( exons_here) > 0) {
    n_genes <- seq(1, length( unique(exons_here$gene.name) ) )
    gene_name_df <- data.frame( x= 0.2*total_length,#(n_genes * total_length) / (max(n_genes) + 1),
                                y=YLIMN,
                                label=rev(unique(exons_here$gene.name))
    )
    # count the occurences of each gene's exons and order by it
    gene_name_df$N <- table(exons_here$gene.name)[gene_name_df$label]
    gene_name_df <- gene_name_df[ order(gene_name_df$N, decreasing = TRUE), ]
    
    # for gene name colouring (black by default, other colours for multiple genes)
    cbbPalette <- "#000000" #, "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    names(cbbPalette) <- gene_name_df$label[1]
    
    # add junction colours
    mainPalette <- c(cbbPalette, mainPalette)
    #names(mainPalette)[ (length(mainPalette)-1):length(mainPalette) ] <- c("Annotated","Novel")
    
    # fit exons within the cluster scale
    invert_mapping=function(pos){
      if (pos %in% s) coords[as.character(pos)] else
        if (pos < min(s)) my_xlim[1] else
          if (pos > max(s)) my_xlim[2] else {
            w=which( pos < s[2:length(s)] & pos > s[1:(length(s)-1)] )
            stopifnot(length(w)==1)
            coords[w] + (coords[w+1]-coords[w])*(pos - s[w])/(s[w+1]-s[w])
          }
    }
    
    exon_df <- data.frame( x=sapply(exons_here$start,invert_mapping),
                           xend=sapply(exons_here$end,invert_mapping),
                           y=0,
                           yend=0,
                           label = exons_here$gene.name)
    
    # alter exon sizes to conform to a minimum exon length
    exon_df[ (exon_df$xend - exon_df$x) < min_exon_length, ]$xend <- exon_df[ (exon_df$xend - exon_df$x) < min_exon_length, ]$x + min_exon_length
    
    # remove exons that are now duplicates - overlapping exons that extend outside of the plotting space and are squished to same coordinates
    exon_df <- exon_df[ !duplicated( paste(exon_df$x, exon_df$xend) ),]
    
    ## STRANDING
    # pick principal gene - the one that contributes the most exons to the cluster will be at the top of the df
    principal_gene <- gene_name_df$label[1]
    gene_strand <- unique(exons_here[exons_here$gene.name == principal_gene,]$strand)
    if( length(gene_strand)  > 1 & !is.na(gene_strand) ){
      gene_strand <- NULL
    }else{
      # assign strand arrows based on lengths of intron
      exon_intervals <-  intervals::Intervals( matrix(data = c(exon_df$x, exon_df$xend), ncol = 2) )
      #exon_intervals <- intervals::interval_union( exon_intervals )
      intron_intervals <- intervals::interval_complement( exon_intervals )
      intron_intervals <- intron_intervals[ 2:(nrow(intron_intervals)-1),]
      # strand arrows should be placed between exons
      strand_df <- as.data.frame(intron_intervals)
      # remove strand arrows where the introns are too small
      strand_df <- strand_df[ (strand_df$V2 - strand_df$V1) > 0.025 * total_length,]
      strand_df$midpoint <- strand_df$V1 + (strand_df$V2 - strand_df$V1) / 2
      
      group <- c()
      for(i in 1:nrow(strand_df)){
        group <- c(group, rep(i,2))
      }
      # make strand_df
      if( gene_strand == "+"){
        strand_df <- data.frame( x = c(rbind(strand_df$V1, strand_df$midpoint)),
                                 group = group,
                                 y = 0)
      }
      if( gene_strand == "-"){
        strand_df <- data.frame( x = c(rbind(strand_df$midpoint, strand_df$V2)),
                                 group = group,
                                 y = 0)
      }
      
      # add strand position - which way should the arrows go?
      if( gene_strand == "+" ){
        strand_pos <- "last"
      }
      if( gene_strand == "-" ){
        strand_pos <- "first"
      }
      # add strand arrows to plot
      for (i in 1:length(plots) ){
        plots[[i]] <- plots[[i]] +
          geom_line( data = strand_df,
                     aes( x = x, y = y, group = group ), colour = "black", size=1,
                     arrow = arrow(ends = strand_pos, type = "open", angle = 30, length = unit(0.1, units = "inches" )))
      }
    }
    
    # add exons to plots
    for (i in 1:length(plots) ){
      plots[[i]] <- plots[[i]] +
        geom_segment( data=exon_df, aes(x=x,y=y,xend=xend,yend=yend, colour = label), alpha=1, size=6) +
        geom_segment( data = exon_df, aes(x = x, xend = x+0.01, y = y, yend = yend), colour = "white", size = 6, alpha = 1) +
        geom_segment( data = exon_df, aes(x = xend-0.01, xend=xend, y = y, yend = yend), colour = "white", size = 6, alpha = 1)
    }
  }
  # TITLE
  
  # give a main_title argument, ideally a vector of c(gene_name, cluster_name)
  if( all( !is.na(main_title) ) ){
    if( length(main_title) > 1){
      plots[[1]] <- plots[[1]] + ggtitle(main_title[1],subtitle = main_title[2]) +
        theme(plot.title = element_text(face="bold.italic", colour="black", size = 10, hjust = 0.45), # centre titles - account for xlabels
              plot.subtitle = element_text(hjust = 0.45))
    }else{
      plots[[1]] = plots[[1]] + ggtitle(main_title)
    }
  }
  
  # add colour palette - different depending on whether exons are included or not - hide in top plot
  plots[[1]] <- plots[[1]] +
    scale_colour_manual("", values = mainPalette ) + theme(legend.position="top", legend.justification = 'right', legend.text=element_text(size=10))
  if (length(plots) > 1){
    for (p in 2:(length(plots))){
      plots[[p]] <- plots[[p]] +
        scale_colour_manual("", values = mainPalette ) + guides(colour=FALSE) # don't show colour legend in top plot
    }
  }
  
  list(plots=plots, start = min(s), end = max(s), expand_width = expand_width, plot_xlim = my_xlim, edges=allEdges2gg)
}

plot_splicing =  function(dir_path, gdb_path, rna_meta, wgs_meta, subject, gene.name, srv_candidate_tbl, exon_table, gene_table, gene.id=NULL, coord=NULL, gene_wise_only = F){
  meta = rna_meta[as.character(rna_meta$SubjectID) == subject , c("SampleID", "Tissue", "Group")] #& Tissue %in% tissues, 
  
  wgs_id = as.character(wgs_meta[as.character(wgs_meta$SubjectID) == subject, "SampleID"])

  ### Load gene intron count table
  if (is.null(gene.id)){
    gene.id = name2id(gene_table, gene.name)[1]
  }
  
  counts_gene = read.table(sprintf("%s%s%s_%s_intron_count.txt.gz", dir_path, .Platform$file.sep, gene.id, gene.name), header=T, sep="\t", stringsAsFactors = F, check.names = F)
  rownames(counts_gene) = paste(counts_gene$chr, counts_gene$start, counts_gene$end, counts_gene$strand, sep=":")
  counts_gene = counts_gene[rowSums(counts_gene[, as.character(meta$SampleID), drop=F]) != 0, ]
  chr = counts_gene$chr[1]
  
  psi_gene = read.table(sprintf("%s%s%s_%s_intron_psi.txt.gz", dir_path, .Platform$file.sep,  gene.id, gene.name), header=T, sep="\t", stringsAsFactors = F, check.names = F)
  rownames(psi_gene) = paste(psi_gene$chr, psi_gene$start, psi_gene$end, psi_gene$strand, sep=":")
  psi_gene = psi_gene[rowSums(psi_gene[, as.character(meta$SampleID), drop=F]) != 0, ]
  
  counts_gene$verdict = "Annotated"
  counts_gene$verdict[counts_gene$pos %in% rownames(psi_gene[psi_gene$is.anno == F, ])] = "Novel"
  counts_gene[counts_gene$verdict == "Novel", "event"] = "annotated"
  
  intron_anno = counts_gene[, !colnames(counts_gene) %in% as.character(rna_meta$SampleID)]
  
  ### Pre define output
  intron_meta = data.frame()
  counts = data.frame()
  psi = data.frame()
  sashimi_plots = list(plots = list(ggplot(), ggplot()), edges=data.frame())
  
  ### Draw local sashimi plot
  if (gene_wise_only == FALSE){
    coord = gsub("chr", "", coord)
    coord = as.numeric(unlist(strsplit(coord, split="[:-]")))
    chr=coord[1]
    plot_start=coord[2]
    plot_end=coord[3]

    counts = counts_gene[counts_gene$start >= plot_start & counts_gene$end <= plot_end, c("chr", "start", "end", "strand", as.character(rna_meta$SampleID))]
    intron_meta = cbind(counts_gene[rownames(counts), !colnames(counts_gene) %in% as.character(rna_meta$SampleID)], counts[, as.character(meta$SampleID)])
    psi = psi_gene[psi_gene$start >= plot_start & psi_gene$end <= plot_end, ]
    

    if (nrow(intron_meta) != 0){
      ### Make RNA splice shashimi-like plots, and get the coordinates of actual plotting region
      sashimi_plots = make_sashimi_like_plot(subject, exons_table = exon_table, meta = meta, intron_meta = intron_meta, counts = counts[rownames(intron_meta), as.character(meta$SampleID)])
      
      local_var = make_variant_plot(gdb_path, subject, gene.name, chr, sashimi_plots$start, sashimi_plots$end, wgs_id, wgs_meta, sashimi_plots$expand_width, sashimi_plots$plot_xlim)
      sashimi_plots$plots = c(sashimi_plots$plots, local_var)

    }else{
      sashimi_plots = list()
      sashimi_plots[["plots"]] = list(empty_ggplot_with_info(sprintf("There is no splice data within %s in gene %s", coord, gene.name)), ggplots())
      sashimi_plots[["edges"]] = data.frame()
    }
    
    # Add junction match variant prediction info
    if (is.data.frame(srv_candidate_tbl)){
      if ("SpliceAI_pred_match_junction" %in% colnames(srv_candidate_tbl)){
        srv_candidate_tbl$Coordinates_of_novel_junction = paste(srv_candidate_tbl$Coordinates_of_novel_junc, srv_candidate_tbl$Strand, sep=":")
        variants_ls = unique(srv_candidate_tbl$DNA_variant)
        variants_ls_colnames = paste0("Match SpliceAI pred of ", variants_ls)
        for (v in variants_ls){
          match_colname = paste0("Match SpliceAI pred of ", v)
          match_junc = unique(subset(srv_candidate_tbl, DNA_variant == v & Coordinates_of_novel_junction %in% rownames(intron_meta), select = c(Coordinates_of_novel_junction, SpliceAI_pred_match_junction)))
          intron_meta[match_junc$Coordinates_of_novel_junction, match_colname] = match_junc$SpliceAI_pred_match_junction
        }
        intron_meta = intron_meta[,c("chr", "start", "end", "strand", "verdict", variants_ls_colnames, as.character(meta$SampleID))]
      }else{
        intron_meta = intron_meta[,c("chr", "start", "end", "strand", "verdict", as.character(meta$SampleID))]
      }
    }else{
      intron_meta = intron_meta[,c("chr", "start", "end", "strand", "verdict", as.character(meta$SampleID))]
    }
    intron_meta = intron_meta %>% dplyr::rename(setNames(as.character(meta$SampleID), as.character(meta$Tissue)))
  }
  
  ### Draw gene wise splice plot
  exons_gene = exon_table[exon_table$gene.name == gene.name & exon_table$chr == chr, ] 
  
  intron_gene = counts_gene[,c("chr", "start", "end", "strand", as.character(meta$SampleID))]
  intron_gene$read_counts = rowSums(intron_gene[, as.character(meta$SampleID), drop=F])
  intron_gene = cbind(intron_gene, intron_anno[rownames(intron_gene), ])
  gene_wise = make_gene_wise_plot(subject, exons_gene, counts_gene, intron_gene)
  
  gene_wise_var = make_variant_plot(gdb_path, subject, gene.name, chr, gene_wise$start, gene_wise$end, wgs_id, wgs_meta, 0, gene_wise$plot_xlim)
  
  gene_wise_plots = list()
  gene_wise_plots[[1]] = gene_wise$plots
  gene_wise_plots = c(gene_wise_plots, gene_wise_var)

  list(local_plots=sashimi_plots$plots, local_table=intron_meta, local_counts=counts, local_psi = psi, local_sashimi2gg = sashimi_plots$edges,
       gene_wise_plots = gene_wise_plots, gene_wise_sashimi2gg = gene_wise$edges)
}

