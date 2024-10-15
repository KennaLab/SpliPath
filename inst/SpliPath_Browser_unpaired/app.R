library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

print(getwd())
ui = fluidPage(
 
   
  navbarPage(
    title = ("SpliPath - csQTL analysis using unpaired genomic and transcriptomic data "),
    tabPanel(
      title = "Gene csQTL",
      value = "sashimiPage",
      
      sidebarLayout(
        sidebarPanel(
          selectInput(inputId = "sqtl_gene", 
                      label = "Search csQTL candidates in gene:",
                      choices = upload_file$view_gene_list),
          
          actionButton(inputId = "plot_gene_splice", 
                       label = "Search")
        ),
        mainPanel(
          verticalLayout(
            titlePanel("Gene csQTL"),
            fluidRow(
              column(12, align="center", plotOutput(outputId = "gene_wise_splice", brush = "zoom_in", height="200px"))
            ),
            fluidRow(
              column(12, align="center", plotOutput(outputId = "gene_wise_var", height="200px"))
            ),
            fluidRow(
              column(12, align="center", DT::dataTableOutput(outputId = "splice_table"))
            )
          )
        )
      )
    ),
    tabPanel(
      title = "Help",
      fluidRow(
        column(8, offset = 2, align="center", img(src='Guide.png', width=1100, align="center"))
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
                      gene_dbscSNV = ggplot(),
                      min_read = 5,
                      # min_psi = 0.2,
                      sqtl_candidate_tbl = data.frame(),
                      intron_sample_tbl = data.frame(),
                      gene_wise_splice = ggplot(),
                      gene_wise_var = list(ggplot()),
                      sqtl_exist = F,
                      sashimi = list(ggplot(), ggplot(), ggplot()),
                      splice_table=data.frame(),
                      splice_counts = data.frame(),
                      sashimi2gg = data.frame())
  

  # 
  observeEvent(input$plot_gene_splice, { 
    splice = plot_splicing(upload_file$browser_data_dir, 
                           upload_file$gdb_path,
                           gene.name=toupper(input$sqtl_gene),
                           exon_table = upload_file$exon_table, 
                           gene_table = upload_file$gene_table,
                           gene_wise_only = T)
    rv$gene_wise_splice = splice$gene_wise_plots[[1]]
    rv$gene_wise_var = splice$gene_wise_plots[-1]
    rv$gene_wise_sashimi2gg = splice$gene_wise_sashimi2gg
    rv$splice_table = splice$local_table
  })
  
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
  output$splice_table = DT::renderDataTable({
    rv$splice_table
  }, server=T, rownames= FALSE, options = list(pageLength = 10, scrollX=TRUE)
  )
  

}

# shiny::runApp(shinyApp(ui = ui, server = server), display.mode="showcase")
shinyApp(ui = ui, server = server)

