library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

print(getwd())
ui = fluidPage(
  
  navbarPage(
    title = ("SpliPath"),
    id = "navbar",
    
    ### 0th Penal: Home
    tabPanel(
      title = "Home",
      fluidRow(
        column(8, offset = 2, align="center", img(src='logo.png', align="center")),
        column(8, offset = 2, align="center", HTML(text = paste(readLines("./www/Introduction.html"), collapse=" ")))
        # column(8, offset = 2, align="center", img(src='Guide.png', width=1000, align="center"))
      )
   ),
    

    ### 3nd Penal: red box
    tabPanel(
      title = "Select csQTL",
      sidebarLayout(
        sidebarPanel(
          selectInput(inputId = "sqtl_gene", 
                    label = "Search csQTL candidates in gene:",
                    choices = upload_file$view_gene_list),
          
          uiOutput("update_red_box_tissues"),
          #selectInput(inputId = "sqtl_tissue", 
          #            label = "Tissue", 
          #            choices = upload_file$tissues), 
          actionButton(inputId = "plot_red_box", 
           label = "Search")
        ),
        mainPanel(
          titlePanel("csQTL candidates"),
          p('Drag cursor to select csQTL candidates in the plots.'),
          splitLayout(
            plotOutput(outputId = "gene_SpliceAI", brush = "gene_spliceai_brush"),
            # plotOutput(outputId = "gene_dbscSNV", brush = "gene_dbscsnv_brush")
          ),
          p('\n\n\nClick row in table to visualize the splicing event and the the DNA variant. \n'),
          DT::dataTableOutput(outputId = "sqtl_candidate_tbl"),
          
          # titlePanel("Other novel junctions in subjects"),
          # p('\n\n\nClick row in table to visualize gene splicing of the subject. \n'),
          # DT::dataTableOutput(outputId = "intron_sample_tbl")
        )
      )
    ),
    
    ### 4rd Penal: Sashimi Plot
    tabPanel(
      title = "sQTL - Sashimi Plot",
      value = "sashimiPage",
      verticalLayout(
        fluidRow(
          column(8, offset = 2, align="left", titlePanel("sQTL splicing"))
        ),
        fluidRow(
          column(8, offset = 2, align="left", p('\nDouble click the junction in the first splicing plot to visualize junction expression in subject groups. \n'))
        ),
        fluidRow(
          column(8, offset = 2, align="center", plotOutput(outputId = "sashimi1st", dblclick = "sashimi_click", height="200px"))
        ),
        fluidRow(
          column(8, offset = 2, align="center", uiOutput("sashimi.ui", height="200px"))
        ), 
        fluidRow(
          column(8, offset = 2, align="center", DT::dataTableOutput(outputId = "splice_table"))
        ),
        
        titlePanel("Gene splicing"),
        fluidRow(
          column(12, align="left", p('\nZoom in by selecting a region in the whole gene splicing plot.'))
        ),
        fluidRow(
          column(12, align="center", plotOutput(outputId = "gene_wise_splice", brush = "zoom_in", height="200px"))
        ),
        fluidRow(
          column(12, align="center", plotOutput(outputId = "gene_wise_var", height="200px"))
        )
      )
    ),
    
    ### 5th Penal: Statistics
    tabPanel(
      title = "sQTL - Expression Plot",
      value = "expPage",
      fluidRow(
        column(8, offset = 2, align="left", titlePanel("Expression of the selected junction among groups")),
        column(8, offset = 2, align="left", p("*Red points refer to junction expression in the query subject's tissues \n"))
      ),
      fluidRow(
        column(8, offset = 2, align="center", uiOutput("compare.ui", height="250px"))
      )
    ),
    tabPanel(
      title = "Help",
      fluidRow(
        column(8, offset = 2, align="center", img(src='Guide.png', align="center"))
     )
   )
   
  )
)

server = function(input, output, session) {
  

  rv = reactiveValues(summary = ggplot(),
                      meta_info = "",
                      meta_table = data.frame(),
                      other_filters = list(),
                      input_overview_genes = c(),
                      input_overview_subject = c(),
                      sample_junc_plot = ggplot(),
                      gene_tissue_tbl = data.frame(),
                      gene_tissue_venn_plot = ggplot(),
                      gene_SpliceAI = ggplot(),
                      min_read = 5,
                      sqtl_candidate_tbl = data.frame(),
                      # intron_sample_tbl = data.frame(),
                      gene_wise_splice = ggplot(),
                      gene_wise_var = list(ggplot()),
                      sqtl_exist = F,
                      sashimi = list(ggplot(), ggplot(), ggplot()),
                      splice_table=data.frame(),
                      splice_counts = data.frame(),
                      sashimi2gg = data.frame())
  
    
  ### 3nd panel
  output$update_red_box_tissues = renderUI({
    selectInput(inputId = "sqtl_tissue", 
                label = "Tissue", 
                choices = upload_file$tissues)
    
  })

  # Select gene and tissue, show red box
  observeEvent(input$plot_red_box, { 
    red_box = draw_gene_red_box(toupper(input$sqtl_gene), upload_file$gene_table, upload_file$browser_data_dir, input$sqtl_tissue, upload_file$rna_meta, 5, "rv$min_psi", min_spliceai = 0.2) 
    rv$gene_SpliceAI = red_box$plots$SpliceAI
    # rv$gene_dbscSNV = red_box$plots$dbscSNV
    rv$sqtl_candidate_tbl = red_box$sqtl_candidate
    rv$select_sqtl_candidate = data.frame()
    # rv$intron_sample_tbl = red_box$intron_sample
  })

  output$gene_SpliceAI = renderPlot({
    rv$gene_SpliceAI
  })
  # output$gene_dbscSNV = renderPlot({
  #   rv$gene_dbscSNV
  # })
  
  # Click gene red box plot, show sqtl candidates table
  observeEvent(input$gene_spliceai_brush, {
    session$resetBrush("gene_spliceai_brush")
    rv$select_sqtl_candidate = select_sqtl_candidate(rv$sqtl_candidate_tbl, input$sqtl_tissue, "SpliceAI", 
                                                     input$gene_spliceai_brush$xmin, input$gene_spliceai_brush$xmax, input$gene_spliceai_brush$ymin, input$gene_spliceai_brush$ymax)
    
  })

  output$sqtl_candidate_tbl = DT::renderDataTable({
    rv$select_sqtl_candidate[, colnames(rv$select_sqtl_candidate) %in%  c("SubjectID", "Coordinates_of_unannotated_junc", "Event", "DNA_variant", "SpliceAI", "SpliceAI_pred_match_junction", "AF", "SpliceAI_pred_cryptic_exon")]
  }, server=T, selection = "single", rownames= FALSE, options = list(pageLength = 10, scrollX=TRUE))
  
  
  # Click sqtl table, plot local sashimi and gene splice
  observeEvent(input$sqtl_candidate_tbl_rows_selected,{
    row_select = as.numeric(input$sqtl_candidate_tbl_rows_selected)
    splice = plot_splicing(upload_file$browser_data_dir, upload_file$gdb_path, upload_file$rna_meta, upload_file$wgs_meta, 
                           subject=rv$select_sqtl_candidate[row_select, "SubjectID"], 
                           coord=rv$select_sqtl_candidate[row_select, "Var_region"], 
                           gene.name=rv$select_sqtl_candidate[row_select, "Gene"], 
                           rv$select_sqtl_candidate,
                           exon_table = upload_file$exon_table, 
                           gene_table = upload_file$gene_table, 
                           gene.id=rv$select_sqtl_candidate[row_select, "Gene_id"]) 
    rv$sqtl_gene = rv$select_sqtl_candidate[row_select, "Gene"]
    rv$sashimi = splice$local_plots
    rv$splice_table = splice$local_table
    rv$splice_counts = splice$local_counts
    rv$splice_psi = splice$local_psi
    rv$sashimi2gg = splice$local_sashimi2gg
    rv$subject_id = rv$select_sqtl_candidate[row_select, "SubjectID"]
    
    rv$gene_wise_splice = splice$gene_wise_plots[[1]]
    rv$gene_wise_var = splice$gene_wise_plots[-1]
    rv$gene_wise_sashimi2gg = splice$gene_wise_sashimi2gg

    updateNavbarPage(session, 'navbar', selected = 'sashimiPage')
  })
  
  # # Show cryptic intron count per subject
  # output$intron_sample_tbl = DT::renderDataTable({
  #   rv$intron_sample_tbl
  # }, server=T, selection = "single", rownames= FALSE, options = list(pageLength = 10, scrollX=TRUE))
  # 
  # # Click intron per sample table, plot only gene splice 
  # observeEvent(input$intron_sample_tbl_rows_selected,{
  #   row_select = as.numeric(input$intron_sample_tbl_rows_selected)
  # 
  #   splice = plot_splicing(upload_file$browser_data_dir, upload_file$gdb_path, upload_file$rna_meta, upload_file$wgs_meta, 
  #                          subject=rv$intron_sample_tbl[row_select, "SubjectID"], 
  #                          gene.name=toupper(input$sqtl_gene), 
  #                          srv_candidate_tbl = NA,
  #                          exon_table = upload_file$exon_table, 
  #                          gene_table = upload_file$gene_table, 
  #                          gene_wise_only = T) 
  #   rv$sqtl_gene = input$sqtl_gene
  #   rv$sashimi = splice$local_plots
  #   rv$splice_table = splice$local_table
  #   rv$splice_counts = splice$local_counts
  #   rv$splice_psi = splice$local_psi
  #   rv$sashimi2gg = splice$local_sashimi2gg
  #   rv$subject_id = rv$intron_sample_tbl[row_select, "SubjectID"]
  #   
  #   rv$gene_wise_splice = splice$gene_wise_plots[[1]]
  #   rv$gene_wise_var = splice$gene_wise_plots[-1]
  #   rv$gene_wise_sashimi2gg = splice$gene_wise_sashimi2gg
  #   
  #   updateNavbarPage(session, 'navbar', selected = 'sashimiPage')
  # })
  
  ### 4rd panel
  
  # Gene wise plot
  gene_var_height = reactive({
    (length(rv$gene_wise_var))*200
  })
  
  output$gene_wise_splice = renderPlot({
    rv$gene_wise_splice
  }, height=200)
  output$gene_wise_var = renderPlot({
    gridExtra::grid.arrange(grobs = rv$gene_wise_var, ncol=1)
  })
  
  # Zoom in 
  observeEvent(input$zoom_in, { 
    
    locate_zoom_in=function(sashimi2gg, zoom_xmin, zoom_xmax){
      zoom_xmin = round((zoom_xmin - sashimi2gg$start[1])*100 + sashimi2gg$startv[1])
      zoom_xmax = round((zoom_xmax - sashimi2gg$start[1])*100 + sashimi2gg$startv[1])
      paste(sashimi2gg$chr[1], zoom_xmin, zoom_xmax, sep=":")
    }

    splice = plot_splicing(upload_file$browser_data_dir, upload_file$gdb_path, upload_file$rna_meta, upload_file$wgs_meta, 
                           subject=rv$subject_id, 
                           gene.name=toupper(rv$sqtl_gene),
                           srv_candidate_tbl = rv$sqtl_candidate_tbl,
                           exon_table = upload_file$exon_table, 
                           gene_table = upload_file$gene_table, 
                           coord=locate_zoom_in(rv$gene_wise_sashimi2gg, input$zoom_in$xmin, input$zoom_in$xmax))
    rv$sashimi = splice$local_plots
    rv$splice_table = splice$local_table
    rv$splice_counts = splice$local_counts
    rv$splice_psi = splice$local_psi
    rv$sashimi2gg = splice$local_sashimi2gg
    
    session$resetBrush("zoom_in")
  })
  
  # Local splicing plot
  sashimi_height = reactive({
                    (length(rv$sashimi) - 1)*150
                  })
  
  output$sashimi1st = renderPlot({
    rv$sashimi[[1]]
  }, height=200)
  output$sashimi = renderPlot({
    gridExtra::grid.arrange(grobs = rv$sashimi[-1], ncol=1)
  })
  output$sashimi.ui = renderUI({
    plotOutput("sashimi", height = sashimi_height())
  })
  output$splice_table = DT::renderDataTable({
    rv$splice_table
  }, server=T, selection = "single", rownames= FALSE, options = list(pageLength = 10, scrollX=TRUE)
  )

  ### 5th panel
  rv_3rd = reactiveValues(compare_plots = list(ggplot(), ggplot()))
  
  compare_height = reactive({
    (ceiling(length(rv_3rd$compare_plots)/4))*600
  })
  
  observeEvent(input$sashimi_click, { 
    rv_3rd$compare_plots = compare_splice(input$sashimi_click$x, input$sashimi_click$y, rv$sashimi2gg, rv$splice_counts, rv$splice_psi, upload_file$rna_meta, rv$subject_id)
    updateNavbarPage(session, 'navbar', selected = 'expPage')
  })
  
  observeEvent(input$splice_table_rows_selected,{
    row_select = as.numeric(input$splice_table_rows_selected)
    specific_junction = paste(rv$splice_table$chr, rv$splice_table$start, rv$splice_table$end, rv$splice_table$strand, sep=":")[row_select]
    rv_3rd$compare_plots = compare_splice(NA, NA, NA, rv$splice_counts, rv$splice_psi, upload_file$rna_meta, rv$subject_id, specific_junction=specific_junction)
    updateNavbarPage(session, 'navbar', selected = 'expPage')
  })
    
  output$compare = renderPlot({
   gridExtra::grid.arrange(grobs = rv_3rd$compare_plots, ncol=2)
  })
  output$compare.ui = renderUI({
    plotOutput("compare", height = compare_height())
  })

}

# shiny::runApp(shinyApp(ui = ui, server = server), display.mode="showcase")
shinyApp(ui = ui, server = server)

