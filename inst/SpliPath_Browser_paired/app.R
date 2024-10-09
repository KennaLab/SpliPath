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
        column(8, offset = 2, align="center", HTML(text = paste(readLines("./www/Introduction.html"), collapse=" "))),
        # column(8, offset = 2, align="center", img(src='Guide.png', width=1000, align="center"))
      )
   ),
    
    # ### 1st Penal: Sample Overview
    # tabPanel(
    #   title = "Sample Overview",
    #   sidebarLayout(
    #     sidebarPanel(
    # 
    #       uiOutput("add_variabbles_from_meta"),
    #     ),
    #     
    #     mainPanel(
    #       titlePanel("Novel junctions in genome"),
    #       plotOutput(outputId = "summary"),
    #       titlePanel("Meta data"),
    #       textOutput(outputId = "meta_info"),
    #       DT::dataTableOutput(outputId = "meta_table")
    #     )
    #   )
    # ),
    # 
    # ### 2nd Penal: Gene overview
    # tabPanel(
    #   title = "Gene Overview",
    #   sidebarLayout(
    #     sidebarPanel(
    # 
    #       p("View novel junctions in query gene(s) or a pre-defined gene set."),
    #       textAreaInput(inputId = "overview_genes", 
    #                 label = "Genes",
    #                 placeholder = "(One entry per line)"),
    #       fileInput(inputId = "overview_geneset", 
    #                 label = "Upload gene set"), 
    #       textAreaInput(inputId = "overview_subject", 
    #                 label = "SubjectIDs",
    #                 placeholder = "(One entry per line)"),
    #       p("*If SubjectIDs are not provided, all the individuals in the meta data will be included."),
    #       radioButtons(inputId = "junction_type_gene", 
    #                    label = "Novel junctions in genes", 
    #                    choices = c("ur-sQTL novel junctions", "All novel junctions"), 
    #                    selected = "ur-sQTL novel junctions"),
    #       
    #       selectInput(inputId = "prediction_tool", 
    #                   label = "Variants annotation methods", 
    #                   choices = c("SpliceAI", "dbscSNV")),
    #       numericInput(inputId = "splice_score_thres", 
    #                    label = "Min prediction score",
    #                    value = 0.2,
    #                    min = 0, max = 1, step=0.1),
    # 
    #       # numericInput(inputId = "psi_thres", 
    #       #              label = "Min junction PSI",
    #       #              value = 0.2,
    #       #              min = 0, max = 1, step=0.1),
    #       
    #       numericInput(inputId = "rc_thres", 
    #                    label = "Min junction read count",
    #                    value = 5,
    #                    min = 0, max = 100, step=1),
    #       numericInput(inputId = "nr_tissues", 
    #                    label = "Min number of concordant tissues",
    #                    value = 1,
    #                    min = 1, max = 4, step=1),
    #       actionButton(inputId = "search_gene", 
    #                      label = "Update")
    #     ),
    #     mainPanel(
    #       titlePanel("Novel junctions in genes"),
    #       DT::dataTableOutput(outputId = "gene_tissue_tbl"),
    #       plotOutput(outputId = "sample_junc_plot"),
    #       titlePanel("Concordance Across Tissues"),
    #       plotOutput(outputId = "gene_tissue_venn_plot")
    #     )
    #   )
    # ),
    
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
    )
  )
)

server = function(input, output, session) {
  
  # ### 1st panel
  # 
  # output$add_variabbles_from_meta = renderUI({
  #   if (nrow(upload_file$wgs_meta) > 0){  #input[["upload_data"]]){
  #     fluidRow(width = 12,
  #              column(12, h5("Plot number of novel junctions by:")),
  #              column(12,
  #                     selectInput(inputId = "tissue", 
  #                          label = "Tissue", 
  #                          choices = c("All_Tissues", upload_file$tissues),
  #                          selected = "All_Tissues")),
  #              column(12,
  #                     selectInput(inputId = "color_by", 
  #                          label = "Color by", 
  #                          choices = setdiff(colnames(upload_file$rna_meta), c("SampleID", "SubjectID"))[1])),
  #              column(12, h5("Filter samples by other info:")),
  #              column( 8,
  #                      selectInput(
  #                        "addvariables_select_tab1_main",
  #                        label = NULL,
  #                        choices = setdiff(colnames(upload_file$rna_meta), c("SampleID", "SubjectID")),
  #                      )
  #              ),
  #              column( 4, actionButton("addvariables_add_tab1_main", label = "Add") ## adds it to filtertab1_m,ain 
  #              ),
  #              column(12,uiOutput("addvariables_ui_tab1_main")),
  #              column(4, align="left", actionButton(inputId = "search_meta", label = "Update"))
  #             )
  #   }
  # })
  # 
  # # Add filters:
  # filters_tab1_main = reactiveValues(
  #   filters = list(),
  #   add_filters = list()
  # )
  # 
  # observeEvent(input[["addvariables_add_tab1_main"]], {
  #                filters_tab1_main[["add_filters"]] = c(filters_tab1_main[["add_filters"]], input[["addvariables_select_tab1_main"]]) 
  #                updateSelectInput(session,
  #                                  inputId = "addvariables_select_tab1_main",
  #                                  choices = setdiff(colnames(upload_file$rna_meta), c("SampleID", "SubjectID")))
  #              })
  # 
  # output$addvariables_ui_tab1_main = renderUI({
  #   if (length(filters_tab1_main[["add_filters"]]) > 0) {
  #     lapply(
  #       filters_tab1_main[["add_filters"]], # user-provided fields that the user wants to filter on 
  #       FUN = function(x) {
  #         type = class(upload_file$rna_meta[[x]])
  #         dat = isolate(filters_tab1_main[["filters"]])
  #         if (type %in% c("numeric", "integer")) {
  #           return(
  #             fluidRow(column(6,
  #               numericInput(
  #                 inputId = sprintf("%s_tab1_main_min", x),
  #                 value = min(upload_file$rna_meta[[x]], na.rm = T),
  #                 label = sprintf("%s_min", x),
  #                 min = min(upload_file$rna_meta[[x]], na.rm = T),
  #                 max = max(upload_file$rna_meta[[x]], na.rm = T)
  #               )
  #             ),
  #             column(6,
  #               numericInput(
  #                 inputId = sprintf("%s_tab1_main_max", x),
  #                 value = max(upload_file$rna_meta[[x]], na.rm = T),
  #                 label = sprintf("%s_max", x),
  #                 min = min(upload_file$rna_meta[[x]], na.rm = T),
  #                 max = max(upload_file$rna_meta[[x]], na.rm = T)
  #               )
  #           )))
  #         } else if (type %in% c("character", "logical", "factor", "Rle")) {
  #             choices = unique(upload_file$rna_meta[[x]])
  #           return(
  #             selectInput(
  #               inputId = sprintf("%s_tab1_main", x),
  #               label = x,
  #               choices = choices,
  #               selected = choices[1]
  #             )
  #           )
  #         }
  #       }
  #     )
  #   }
  # })  
  # 
  # rv = reactiveValues(summary = ggplot(), 
  #                     meta_info = "",
  #                     meta_table = data.frame(),
  #                     other_filters = list(),
  #                     input_overview_genes = c(),
  #                     input_overview_subject = c(),
  #                     sample_junc_plot = ggplot(),
  #                     gene_tissue_tbl = data.frame(),
  #                     gene_tissue_venn_plot = ggplot(),  
  #                     gene_SpliceAI = ggplot(),
  #                     # gene_dbscSNV = ggplot(),
  #                     min_read = 5,
  #                     # min_psi = 0.2,
  #                     sqtl_candidate_tbl = data.frame(),
  #                     intron_sample_tbl = data.frame(),
  #                     gene_wise_splice = ggplot(),
  #                     gene_wise_var = list(ggplot()),
  #                     sqtl_exist = F,
  #                     sashimi = list(ggplot(), ggplot(), ggplot()),
  #                     splice_table=data.frame(),
  #                     splice_counts = data.frame(),
  #                     sashimi2gg = data.frame())
  # 
  # ### Default Show
  # cryptic_intron = plot_number_cryptic_intron(upload_file$gene_splice, upload_file$rna_meta, setdiff(colnames(upload_file$rna_meta), c("SampleID", "SubjectID"))[1], "All_Tissues", other_filters=list()) 
  # rv$summary = cryptic_intron$summary_plot
  # rv$meta_table = cryptic_intron$meta_table
  # 
  # observeEvent(input$search_meta, { 
  #   # Put selected values for filters into `other_filters`
  #   rv$other_filters[["add_filters"]] = filters_tab1_main[["add_filters"]]
  #   for (x in rv$other_filters[["add_filters"]]){
  #     type = class(upload_file$rna_meta[[x]])
  #     if (type %in% c("numeric", "integer")) {
  #       value_idx = sprintf("%s_tab1_main_min", x)
  #       rv$other_filters[[value_idx]] = input[[value_idx]]
  #       value_idx = sprintf("%s_tab1_main_max", x)
  #       rv$other_filters[[value_idx]] = input[[value_idx]]
  #       
  #     } else if (type %in% c("character", "logical", "factor", "Rle")) {
  #       value_idx = sprintf("%s_tab1_main", x)
  #       rv$other_filters[[value_idx]] = input[[value_idx]]
  #     }
  #   }
  #   cryptic_intron = plot_number_cryptic_intron(upload_file$gene_splice, upload_file$rna_meta, input$color_by, input$tissue, rv$other_filters) 
  #   rv$summary = cryptic_intron$summary_plot
  #   rv$meta_table = cryptic_intron$meta_table
  # })
  # output$summary = renderPlot({
  #   rv$summary
  # })
  # 
  # p("\n\n\nMeta data of the selected samples.\n")
  # output$meta_table = DT::renderDataTable({
  #   rv$meta_table
  # }, options = list(pageLength = 20, select = TRUE, scrollX=TRUE), rownames= FALSE
  #    #extensions =  c("Select"), selection = list(target = "column"), 
  # )
  # 
  # ### 2nd panel
  # 
  # observeEvent(input$search_gene, {
  #   rv$min_read = input$rc_thres
  #   # rv$min_psi = input$psi_thres
  #   
  #   rv$input_overview_genes = c()
  #   if (! "" %in% input$overview_genes){
  #     rv$input_overview_genes = c(rv$input_overview_genes, setdiff(unlist(strsplit(input$overview_genes, split="\n")), ""))
  #   }
  #   if (!is.null(input$overview_geneset)){
  #     rv$input_overview_genes = c(rv$input_overview_genes, readLines(input$overview_geneset$datapath))
  #   }
  # 
  #   rv$input_overview_subjects = input$overview_subject
  #   if (! "" %in% input$overview_subject ){
  #     rv$input_overview_subjects = setdiff(unlist(strsplit(input$overview_subject, split="\n")), "")
  #   }
  #   
  #   gene_overview = gene_tissue_nr_junc_venn(toupper(rv$input_overview_genes), upload_file$tissues, rv$input_overview_subjects, dir_path = upload_file$browser_data_dir, gene_table = upload_file$gene_table, rna_meta = upload_file$rna_meta, "Group", 
  #                                            "input$psi_thres", input$rc_thres, input$nr_tissues, input$junction_type_gene, input$prediction_tool, input$splice_score_thres)
  #   rv$sample_junc_plot = gene_overview$sample_junc_plot
  #   rv$gene_tissue_tbl = gene_overview$gene_tissue_tbl
  #   rv$gene_tissue_venn_plot = gene_overview$gene_tissue_venn_plot  
  #     
  # })
  # 
  # output$sample_junc_plot = renderPlot({
  #   rv$sample_junc_plot
  # })
  # output$gene_tissue_tbl = DT::renderDataTable({
  #   rv$gene_tissue_tbl
  # }, rownames= FALSE)
  # output$gene_tissue_venn_plot = renderPlot({
  #   rv$gene_tissue_venn_plot
  # })
  
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
  # observeEvent(input$gene_dbscsnv_brush, {
  #   session$resetBrush("gene_spliceai_brush")
  #   rv$select_sqtl_candidate = select_sqtl_candidate(rv$sqtl_candidate_tbl, input$sqtl_tissue, "dbscSNV", 
  #                                                    input$gene_dbscsnv_brush$xmin, input$gene_dbscsnv_brush$xmax, input$gene_dbscsnv_brush$ymin, input$gene_dbscsnv_brush$ymax)
  #   
  # })
  
  output$sqtl_candidate_tbl = DT::renderDataTable({
    rv$select_sqtl_candidate[, colnames(rv$select_sqtl_candidate) %in%  c("SubjectID", "Coordinates_of_novel_junc", "Event", "DNA_variant", "SpliceAI", "SpliceAI_pred_match_junction", "AF", "SpliceAI_pred_cryptic_exon", "csQTL_candidate")]
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

shiny::runApp(shinyApp(ui = ui, server = server), display.mode="showcase")
# shinyApp(ui = ui, server = server)

