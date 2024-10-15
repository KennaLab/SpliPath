make_gene_wise_plot <- function(
    exons_gene = NULL,
    intron_meta = NULL
    ){
  
  summary_func <- colSums
  
  ### Get unique intron start and end for docking
  intron_meta$id=1:nrow(intron_meta) # number each junction
  temp=intron_meta[,c("id","start","end")]
  m=reshape2::melt(temp, id.vars = "id") # melt to get list of all coordinates
  s=unique(c(m$value, exons_gene$start, exons_gene$end)) # get unique start and end values, plot from the first exon to the last exon
  s=sort(s)
  d=s[2:length(s)]-s[1:length(s)-1] # get the difference between pairs of coordinates
  trans_d <- length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33 # apply a trasnformation function - doesn't work on negative numbers!
  coords <- c(0,cumsum(trans_d))
  names(coords)=s
  
  total_length=sum(trans_d) 
  my_xlim=c(0, total_length)

  # PLOT SETTINGS
  
  new_theme_empty <- theme_bw(base_size = 11 )
  new_theme_empty$panel.background = element_rect(fill="white", colour = "white")
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
  main_title = paste0("Unannotated junctions in gene ", exons_gene$gene.name[1])
  legend_title <- "Read"
  junction_colour <- "red"
  cryptic_colour <- "grey"
  mainPalette <- c(junction_colour, cryptic_colour)
  names(mainPalette) = c("Annotated", "Unannotated")
  
  min_height=0
  max_height=0
  curv = 0.15
  min_exon_length = 0.5
  maxcount=0
  mincount=1.0
  yFactor = 0.65   # originally set as 0.65
  yConstant = -0.35 # originally set as 0.5
  labelTextSize= 4 # orignally set as 5
  curveMax = 10
  curveExponent = 2
  yOffset = 0
  centreLineWidth = 3 
  
  maxcount=max(c(max(intron_meta$read_counts), maxcount)) # find max and min values of summarised counts
  mincount=min(c(min(intron_meta$read_counts), mincount))

  ### Define junction curve parameters and draw gene wise plot
  intron_meta$prop=intron_meta$read_counts/sum(intron_meta$read_counts)  # this has been changed
  
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
    edge$clu<- intron_meta$clu[i]
    edge$Group <- i
    edge$xtext <-start+l/2
    edge$ytext <- -( (l^(yFactor) / 2) + yConstant)  # magic formula here
    edge$verdict <- ifelse( intron_meta$verdict[i] == "Annotated", yes = "Annotated", no ="Unannotated")
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
    edge$Group <- i
    edge$xtext <-start+l/2
    edge$ytext <- (l^(yFactor) / 2)  + yConstant
    edge$verdict <- ifelse( intron_meta$verdict[i] == "Annotated", yes = "Annotated", no ="Unannotated")
    edge
  })
  
  allEdges2gg = rbind(allEdgesP, allEdges)
  allEdges2gg$curv = curv
  allEdges2gg$yOffset = yOffset
  allEdges2gg$curveMax = curveMax

  MAXcounts <- max(c(max(with(allEdgesP,1)),max(with(allEdges,1))))
  
  YLIMP <- max( allEdgesP$ytext) + 0.25 * max( allEdgesP$ytext)
  YLIMN <- min( allEdges$ytext) + 0.25 * min( allEdges$ytext)
  print(YLIMP); print(YLIMN)
  gene_wise_plot <- ggplot() +
    geom_curve(data=allEdgesP, aes(x = start, xend = xtext, y = yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)*10 ),
               angle=90, curvature=-curv,lineend="round") +
    geom_curve(data=allEdgesP, aes(x = xtext, xend = end, y = ytext, yend = yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)*10 ),
               angle=90, curvature=-curv,lineend="round") +
    geom_curve(data=allEdges, aes(x = start, xend = xtext, y = -yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)*10 ),
               angle=90,curvature=curv,lineend="round") +
    geom_curve(data=allEdges, aes(x = xtext, xend = end, y = ytext, yend = -yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)*10 ),
               angle=90,curvature=curv,lineend="round") +
    
    new_theme_empty +
    ylab("") +
    xlab("") +
    xlim(my_xlim) +
  
    # horizontal line - smooth out the ends of the curves
    geom_hline(yintercept=0, size = centreLineWidth, colour = "white") +
    geom_hline(yintercept=0,alpha=.9, size=1) +
    
    # label the junctions
    ylim(YLIMN,YLIMP) + 
    guides(size=FALSE)
    # scale_size_continuous(limits=c(0,10),guide='none')
    
    
  # ADDING EXON ANNOTATION
  df <- data.frame(x=coords, xend=total_length*(s-min(s))/(max(s)-min(s)), y=0, yend=min_height) # coords: the position of splice sites in the plot

  exons_here <- exons_gene # find exons
  
  # if exons found - remove exons that don't actually start or end with a junction
  # if( nrow(exons_here) > 0 ){
  #   exons_here <-  unique(
  #     exons_here[ ( exons_here$end %in% intron_meta$start |
  #                     exons_here$start %in% intron_meta$end ), ]
  #   )
  # }

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
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # sort out colour pallette - most common gene goes first!
    cbbPalette <- cbbPalette[ 1:length(gene_name_df$label) ]
    names(cbbPalette) <- gene_name_df$label
    # add junction colours
    mainPalette <- c(cbbPalette, mainPalette)

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
    
    # add exons to plots
    gene_wise_plot = gene_wise_plot +
      geom_segment( data = exon_df, aes(x=x,y=y,xend=xend,yend=yend, colour = label), alpha=1, size=6) +
      geom_segment( data = exon_df, aes(x = x, xend = x+0.01, y = y, yend = yend), colour = "white", size = 6, alpha = 1) +
      geom_segment( data = exon_df, aes(x = xend-0.01, xend=xend, y = y, yend = yend), colour = "white", size = 6, alpha = 1)
    
  }


  # TITLE
  # give a main_title argument, ideally a vector of c(gene_name, cluster_name)
  gene_wise_plot = gene_wise_plot + ggtitle(main_title)

  # add colour palette - different depending on whether exons are included or not - hide in top plot
  gene_wise_plot = gene_wise_plot +
    scale_colour_manual("", values = mainPalette ) + theme(legend.position="top", legend.justification = 'right', legend.text=element_text(size=10))

  list(plots=gene_wise_plot, start = min(s), end = max(s), plot_xlim = my_xlim, edges=allEdges2gg)
}

