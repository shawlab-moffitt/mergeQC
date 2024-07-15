type_id <- paste0("v1.0.20240711")

library("shiny")
library("DT")
library("data.table")
library("RecordLinkage")
library("dplyr")
library("reshape2")
library("ggplot2")
library("shinyjqui")
library("shinycssloaders")
library("plotly")
library("svglite")

options(shiny.maxRequestSize = 10000*1024^2)

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab, na.rm = T)]
}

# UI ---------------------------------------------------------------------------

# Define UI for application that draws a histogram

ui <- navbarPage("{ MergeQC }",
                 ## Data Preview -----------------------------------------------
                 tabPanel("Data Preview",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("MatFileInput", label = "Matrix Files:",multiple = T, placeholder = "CTRL + left click files"),
                              fileInput("MetFileInput", label = "Meta Files:",multiple = T, placeholder = "CTRL + left click files")
                            ),
                            
                            # Show a plot of the generated distribution
                            mainPanel(
                              verbatimTextOutput("FileCheckAlerts"),
                              fluidRow(
                                column(2, #style = 'padding-right:0px',
                                       uiOutput("rendDnldMatrix")
                                ),
                                column(2, #style = 'padding-left:0px;padding-right:0px',
                                       uiOutput("rendDnldMeta")
                                ),
                                column(2, #style = 'padding-left:0px;padding-right:0px',
                                       uiOutput("rendDnldMatrixFull")
                                ),
                                column(2, #style = 'padding-left:0px;padding-right:0px',
                                       uiOutput("rendDnldMetaFull")
                                ),
                                column(2, #style = 'padding-left:0px',
                                       uiOutput("renddownload_notes")
                                )
                              ),
                              p(),
                              tabsetPanel(
                                id = "tabs"
                              )
                            )
                          )
                          ),
                 ## Data QC ----------------------------------------------------
                 tabPanel("Data QC",
                          sidebarLayout(
                            sidebarPanel(
                              width = 4,
                              conditionalPanel(condition = "input.QCtab == '1'",
                                               tabsetPanel(
                                                 tabPanel("Data Input",
                                                          p(),
                                                          selectizeInput("LogDataOpts","Select datasets to log2+1:", choices = NULL, selected = NULL, multiple = T),
                                                          div(DT::dataTableOutput("DataQCtab"), style = "font-size:12px"),
                                                          downloadButton("dnldDataQCtab","Download Table")
                                                 ),
                                                 tabPanel("Figure Settings",
                                                 )
                                               )
                                               ),
                              conditionalPanel(condition = "input.QCtab == '2'",
                                               tabsetPanel(
                                                 tabPanel("Data Input",
                                                          p(),
                                                          selectizeInput("AvgExprData1","Dataset 1:",choices = NULL, selected = 1),
                                                          selectizeInput("AvgExprData2","Dataset 2:",choices = NULL, selected = 1),
                                                          selectizeInput("AvgExprGene","Highligh Genes:",choices = NULL, selected = NULL, multiple = T)
                                                          ),
                                                 tabPanel("Figure Settings",
                                                 )
                                               )
                                               )
                            ),
                            mainPanel(
                              tabsetPanel(id = "QCtab",
                                          tabPanel("Expression Density",
                                                   p(),
                                                   shinyjqui::jqui_resizable(plotOutput("DensityPlot", width = "100%", height = "500px")),
                                                   downloadButton("dnldDensityPlot_SVG","SVG"),
                                                   p(),
                                                   div(DT::dataTableOutput("AvgExprDens"), style = "font-size:12px"),
                                                   downloadButton("dnldAvgExprDens","Download Table"),
                                                   value = 1),
                                          tabPanel("Average Expression",
                                                   p(),
                                                   shinyjqui::jqui_resizable(plotlyOutput("AvgExpr_Plot", width = "100%", height = "500px")),
                                                   downloadButton("dnldAvgExpr_Plot_SVG","SVG"),
                                                   p(),
                                                   div(DT::dataTableOutput("AvgExpr_df"), style = "font-size:12px"),
                                                   downloadButton("dnldAvgExpr_df","Download Table"),
                                                   value = 2)
                                          )
                            )
                          )
                          )
                 )



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  FileCheckAlerts_react <- reactiveVal()
  
  output$FileCheckAlerts <- renderPrint({
    
    req(FileCheckAlerts_react())
    text <- paste(FileCheckAlerts_react(), collapse = "\n")
    cat(text)
    
  })
  output$renddownload_notes <- renderUI({
    if (length(FileCheckAlerts_react()) > 0) {
      downloadButton("download_notes","Download data processing notes")
    }
  })
  output$download_notes <- downloadHandler(
    filename = function() {
      paste("Merge_Data_Processing_Notes_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      df <- FileCheckAlerts_react()
      writeLines(df, file)
    }
  )
  
  # Load in files --------------------------------------------------------------
  
  observeEvent(input$MatFileInput, {
    
    prependTab(inputId = "tabs",
              tabPanel("Inner Merged Matrix",
                       uiOutput("rendMatrixHeader"),
                       uiOutput("rendMatHead"),
                       dataTableOutput("Matrix_Merge_Preview")
                       ),
              select = T
              )
    
  })
  
  observe({
    
    req(matrix_files_react())
    appendTab(inputId = "tabs",
              tabPanel("Full Merged Matrix",
                       uiOutput("rendFullMatrixHeader"),
                       uiOutput("rendFullMatHead"),
                       dataTableOutput("FullMatrix_Merge_Preview")
              )
    )
    
  })
  
  matrix_files_df_react <- reactive({
    
    file_df <- input$MatFileInput
    file_df
    
  })
  observe({
    req(matrix_files_df_react())
    updateSelectizeInput(session,"LogDataOpts",choices = tools::file_path_sans_ext(matrix_files_df_react()[,"name"]), selected = NULL)
  })
  
  matrix_files_react <- reactive({
    
    req(matrix_files_df_react())
    file_df <- matrix_files_df_react()
    logData <- input$LogDataOpts
    if (nrow(file_df) == 1) {
      if (tools::file_ext(file_df$datapath[1]) %in% c("zip","ZIP","gz","GZ")) {
        filelist <- unzip(file_df$datapath, list = T)
        df_list <- lapply(filelist$Name, function(x) as.data.frame(fread(x)))
        df_list
      }
    } else if (nrow(file_df) > 1) {
      df_list <- list()
      for (row in seq(nrow(file_df))) {
        df <- as.data.frame(fread(file_df[row,"datapath"]))
        colnames(df)[1] <- "Gene"
        tabName <- tools::file_path_sans_ext(file_df[row,"name"])
        if (tabName %in% logData) {
          df[,-1] <- log2(as.matrix(df[,-1]+1))
        }
        df_list[[tabName]] <- df
      }
      df_list
    }
    
  })
  
  
  observeEvent(input$MetFileInput, {
    
    req(matrix_files_react())
    appendTab(inputId = "tabs",
              tabPanel("Full Merged Meta",
                       uiOutput("rendMetaHeader"),
                       uiOutput("rendMetHead"),
                       dataTableOutput("Meta_Merge_Preview")
                       )
              )
    
  })
  
  
  
  meta_files_df_react <- reactive({
    
    file_df <- input$MetFileInput
    file_df
    
  })
  
  meta_files_react <- reactive({
    
    req(meta_files_df_react())
    file_df <- meta_files_df_react()
    if (nrow(file_df) == 1) {
      if (tools::file_ext(file_df$datapath[1]) %in% c("zip","ZIP","gz","GZ")) {
        filelist <- unzip(file_df$datapath, list = T)
        df_list <- lapply(filelist$Name, function(x) as.data.frame(fread(x)))
        df_list
      }
    } else if (nrow(file_df) > 1) {
      df_list <- list()
      for (row in seq(nrow(file_df))) {
        df <- as.data.frame(fread(file_df[row,"datapath"]))
        colnames(df)[1] <- "OrigMetaID"
        df$File_Name <- file_df[row,"name"]
        df_list[[file_df[row,"name"]]] <- df
      }
      df_list
    }
    
  })
  # Matrix Merge ---------------------------------------------------------------
  
  output$rendDnldMatrix <- renderUI({
    req(matrix_merged_inn())
    downloadButton("dnldMatrix","Inner Merged Matrix")
  })
  
  output$rendDnldMatrixFull <- renderUI({
    req(matrix_merged_full())
    downloadButton("dnldMatrixFull","Full Merged Matrix")
  })
  
  matrix_merged_inn <- reactive({
    
    df_list <- matrix_files_react()
    file_df <- matrix_files_df_react()
    
    masterDF <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2),
                       df_list)
    masterDF
    
  })
  
  matrix_merged_full <- reactive({
    
    df_list <- matrix_files_react()
    file_df <- matrix_files_df_react()
    
    masterDF <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = TRUE),
                       df_list)
    masterDF
    
  })

  output$rendFullMatrixHeader <- renderUI({
    req(matrix_merged_full())
    h3("Full Merged Matrix Preview")
  })
  output$rendFullMatHead <- renderUI({
    req(matrix_merged_full())
    if (ncol(matrix_merged_full()) > 301) {
      radioButtons("FullMatHead",NULL, choices = c("View table head","View entire table"), inline = T)
    } else {
      radioButtons("FullMatHead",NULL, choices = c("View table head","View entire table"), inline = T, selected = "View entire table")
    }
  })
  output$FullMatrix_Merge_Preview <- renderDataTable({
    req(matrix_merged_full())
    req(input$FullMatHead)
    mat <- as.data.frame(matrix_merged_full())
    if (input$FullMatHead == "View table head") {
      mat <- head(mat,c(100,100))
    }
    datatable(mat,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  
  output$rendMatrixHeader <- renderUI({
    req(matrix_merged_inn())
    h3("Inner Merged Matrix Preview")
  })
  output$rendMatHead <- renderUI({
    req(matrix_merged_inn())
    if (ncol(matrix_merged_inn()) > 301) {
      radioButtons("MatHead",NULL, choices = c("View table head","View entire table"), inline = T)
    } else {
      radioButtons("MatHead",NULL, choices = c("View table head","View entire table"), inline = T, selected = "View entire table")
    }
  })
  output$Matrix_Merge_Preview <- renderDataTable({
    req(matrix_merged_inn())
    req(input$MatHead)
    mat <- as.data.frame(matrix_merged_inn())
    if (input$MatHead == "View table head") {
      mat <- head(mat,c(100,100))
    }
    datatable(mat,
              options = list(lengthMenu = c(5,10, 20, 100, 1000),
                             pageLength = 10,
                             scrollX = T),
              rownames = F)
  })
  
  output$dnldMatrixFull <- downloadHandler(
    filename = function() {
      paste0("Full_Merged_Matrix_",Sys.Date(),".txt")
    },
    content = function(file) {
      df <- matrix_merged_full()
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  output$dnldMatrix <- downloadHandler(
    filename = function() {
      paste0("Inner_Merged_Matrix_",Sys.Date(),".txt")
    },
    content = function(file) {
      df <- matrix_merged_inn()
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  
  
  # Meta Merge -----------------------------------------------------------------
  output$rendDnldMeta <- renderUI({
    
    req(meta_merged())
    downloadButton("dnldMeta","Inner Merged Meta")
    
  })
  
  meta_merged <- reactive({
    
    req(matrix_merged_inn())
    df_list <- meta_files_react()
    file_df <- meta_files_df_react()
    matrix_merged_inn <- matrix_merged_inn()
    
    merged_meta <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = TRUE),
                       df_list)

    matrix_colnames <- colnames(matrix_merged_inn)[-1]
    match_scores <- sapply(merged_meta$OrigMetaID, function (x) sapply(matrix_colnames, function(y) levenshteinSim(y, x)))
    merged_meta$BestMatchID <- rownames(match_scores)[apply(match_scores, 2, which.max)]
    merged_meta$BestMatchScore <- apply(match_scores, 2, max)
    merged_meta <- merged_meta %>% dplyr::relocate(File_Name)
    merged_meta <- merged_meta |> dplyr::relocate(BestMatchID)
    merged_meta <- merged_meta |> dplyr::relocate(BestMatchScore, .before = OrigMetaID)
    merged_meta <- merged_meta[match(matrix_colnames, merged_meta$BestMatchID),]
    merged_meta
    
  })
  
  output$rendMetaHeader <- renderUI({
    req(meta_merged())
    h3("Merged Meta Preview")
  })
  output$rendMetHead <- renderUI({
    req(meta_merged())
    if (ncol(meta_merged()) > 301) {
      radioButtons("MetHead",NULL, choices = c("View table head","View entire table"), inline = T)
    } else {
      radioButtons("MetHead",NULL, choices = c("View table head","View entire table"), inline = T, selected = "View entire table")
    }
  })
  output$Meta_Merge_Preview <- renderDataTable({
    req(meta_merged())
    req(input$MetHead)
    met <- as.data.frame(meta_merged())
    if (input$MetHead == "View table head") {
      met <- head(met,c(100,100))
    }
    datatable(met,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  
  output$dnldMeta <- downloadHandler(
    filename = function() {
      paste0("Merged_Meta_",Sys.Date(),".txt")
    },
    content = function(file) {
      df <- meta_merged()
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  
  # Data Processing ------------------------------------------------------------
  
  observe({
    
    SampsDiffInMet <- setdiff(meta_merged()[,1],intersect(meta_merged()[,1], colnames(matrix_merged_full())[-1]))
    SampsDiffInMat <- setdiff(colnames(matrix_merged_full())[-1],intersect(meta_merged()[,1], colnames(matrix_merged_full())[-1]))
    FileCheckAlerts_list <- c()
    FileCheckAlerts_list <- c(FileCheckAlerts_list,
                              paste0("Number of matrix files uploaded: ",length(matrix_files_react())),
                              paste0("Number of meta files uploaded: ",length(meta_files_react())),
                              paste0("Total number of samples from matrix data: ",ncol(matrix_merged_full()[,-1])),
                              paste0("Total number of samples from meta data: ",nrow(meta_merged())),
                              paste0("Total number of features consistent in all matrix data: ",nrow(matrix_merged_inn())))
    if (length(SampsDiffInMat) > 0) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                paste0("Number of samples in matrix data not found in meta data: ",SampsDiffInMat))
    }
    if (length(SampsDiffInMet) > 0) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                paste0("Number of samples in meta data not found in matrix data: ",SampsDiffInMet))
    }
    FileCheckAlerts_react(FileCheckAlerts_list)
    
  })
  
  # QC stats -------------------------------------------------------------------
  
  QC_react <- reactive({
    
    df_list <- matrix_files_react()
    
    df_list_qc <- lapply(df_list,function(x){
      quants <- c()
      quants <- c(quants,min = min(as.matrix(x[,-1]), na.rm = T))
      quants <- c(quants,quantile(as.matrix(x[,-1]),probs = c(0.25, 0.5, 0.75, 0.9, 0.99), na.rm = T))
      quants <- c(quants,max = max(as.matrix(x[,-1]), na.rm = T))
      quants <- c(quants,mode = Modes(as.matrix(x[,-1]))[1])
      return(quants)
    })
    
    df_list_qc_df <- as.data.frame(do.call(rbind,df_list_qc))
    
    df_list_qc_df <- cbind(data.frame(File = rownames(df_list_qc_df)),
                                      df_list_qc_df)
    df_list_qc_df

  })
  
  output$DataQCtab <- DT::renderDataTable({
    
    df <- QC_react()
    df <- df[,-c(6,7)]
    datatable(df,
              options = list(lengthMenu = c(5,10, 20, 100, 1000),
                             pageLength = 10,
                             scrollX = T),
              selection = list(selected = 1),
              rownames = F) %>%
      formatRound(columns = colnames(df)[-1], digits = 3)
    
  })
  output$dnldDataQCtab <- downloadHandler(
    filename = function() {
      paste0("Merged_Expression_QC",Sys.Date(),".txt")
    },
    content = function(file) {
      df <- QC_react()
      df <- df[,-c(6,7)]
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  
  
  # Density plot ---------------------------------------------------------------
  
  DensityPlot_df <- reactive({
    
    df_list <- matrix_files_react()
    
    df_list_avgExpr <- lapply(seq_along(df_list),function(x) {
      df <- df_list[[x]]
      df_name <- names(df_list)[[x]]
      print(head(df,c(5,5)))
      x_mean <- rowMeans(as.matrix(df[,-1]))
      x_mean_df <- data.frame(Gene = df[,1],
                              Avg_Expression = x_mean)
      colnames(x_mean_df)[2] <- df_name
      return(x_mean_df)
    })
    
    df_list_avgExpr_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = T),
                                 df_list_avgExpr)
    df_list_avgExpr_df
  })
  
  output$AvgExprDens <- DT::renderDataTable({
    
    df <- DensityPlot_df()
    datatable(df,
              options = list(lengthMenu = c(5,10, 20, 100, 1000),
                             pageLength = 10,
                             scrollX = T),
              rownames = F) %>%
      formatRound(columns = colnames(df)[-1], digits = 3)
    
    
  })
  output$dnldAvgExprDens <- downloadHandler(
    filename = function() {
      paste0("Merged_Average_Expression_",Sys.Date(),".txt")
    },
    content = function(file) {
      df <- DensityPlot_df()
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  
  DensityPlot_react <- reactive({
    
    if (length(input$DataQCtab_rows_selected) > 0) {
      df_list_avgExpr_df <- DensityPlot_df()
      print(input$DataQCtab_rows_selected)
      files_select <- QC_react()[input$DataQCtab_rows_selected,1]
      xMax <- max(QC_react()[input$DataQCtab_rows_selected,6])
      df_list_avgExpr_df_select <- df_list_avgExpr_df[,c(colnames(df_list_avgExpr_df)[1],files_select)]
      df_list_avgExpr_df_select_melt <- melt(df_list_avgExpr_df_select)
      
      p <- ggplot(df_list_avgExpr_df_select_melt, aes(x=value, fill=variable)) +
        geom_density(alpha=0.4) +
        theme_minimal() +
        xlim(c(min(df_list_avgExpr_df_select_melt[,3]),xMax))+
        labs(x = "Averge Expression", y = "Density",
             color = colnames(df_list_avgExpr_df)[1]) +
        theme(legend.position = "top",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14))
        
      
      p
    }
    
  })
  
  output$DensityPlot <- renderPlot({
    
    p <- DensityPlot_react()
    p
    
  })
  output$dnldDensityPlot_SVG <- downloadHandler(
    filename = function() {
      paste("Merged_Density_Plot_",Sys.Date(),".svg", sep = '')
    },
    content = function(file) {
      p <- DensityPlot_react()
      ggsave(file,p, width = 10, height = 8)
    }
  )
  
  # Average Expr ---------------------------------------------------------------
  
  observe({
    
    req(matrix_files_react())
    datasets <- names(matrix_files_react())
    
    updateSelectizeInput(session,"AvgExprData1", choices = datasets, selected = datasets[1])
    updateSelectizeInput(session,"AvgExprData2", choices = datasets, selected = datasets[2])
    
  })
  
  observe({
    
    genes <- AvgExpr_Plot_df()[,1]
    updateSelectizeInput(session,"AvgExprGene", choices = genes, selected = NULL, server = T)
    
  })
  
  AvgExpr_Plot_df <- reactive({
    
    req(DensityPlot_df())
    req(input$AvgExprData1)
    req(input$AvgExprData2)
    set1 <- input$AvgExprData1
    set2 <- input$AvgExprData2
    df <- DensityPlot_df()
    
    df <- df[,c(colnames(df)[1],set1,set2)]
    df <- df %>%
      mutate(ColorCol = case_when(
        df[,set2] > df[,set1] ~ paste0(set2," Expr > ",set1,' Expr'),
        df[,set1] > df[,set2] ~ paste0(set1," Expr > ",set2,' Expr')
      ))
    
    df
    
    
  })
  
  AvgExpr_Plot_react <- reactive({
    
    req(AvgExpr_Plot_df())
    
    set1 <- input$AvgExprData1
    set2 <- input$AvgExprData2
    plot_df <- AvgExpr_Plot_df()
    genes <- input$AvgExprGene
    
    
    p <- ggplot(plot_df, aes(x = !!sym(set1), y = !!sym(set2),
                             color = ColorCol,
                             text = paste0("</br>Gene: ",Gene,
                                           "</br>",set1,": ",round(!!sym(set1),3),
                                           "</br>",set2,": ",round(!!sym(set2),3)))) +
      geom_point(size = 1) +
      theme_minimal() +
      theme(legend.position="none",
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 20)) +
      xlab(paste0("Average Expression: ",set1))+
      ylab(paste0("Average Expression: ",set2)) +
      labs(title = paste0("Average Expression: ",set1," vs. ",set2)) +
      scale_color_manual(values = c("cadetblue3","lightcoral"))
    
    if (length(genes) > 0) {
      plot_df_genes <- plot_df[which(plot_df[,1] %in% genes),]
      p <- p + geom_point(data = plot_df_genes,
                          aes(x = !!sym(set1), y = !!sym(set2)),
                          pch = 1,
                          color = "black",
                          size = 1,
                          stroke = .3)
    }
    
    p
    
  })
  output$dnldAvgExpr_Plot_SVG <- downloadHandler(
    filename = function() {
      set1 <- input$AvgExprData1
      set2 <- input$AvgExprData2
      paste0(set1,"_",set2,"_Avg_Expression_plot_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- AvgExpr_Plot_react()
      ggsave(file,p, width = 10, height = 8)
    }
  )
  
  output$AvgExpr_Plot <- renderPlotly({
    
    p <- AvgExpr_Plot_react()
    set1 <- input$AvgExprData1
    set2 <- input$AvgExprData2
    plot_df <- AvgExpr_Plot_df()
    ply <- ggplotly(p, tooltip = "text")
    genes <- input$AvgExprGene
    if (length(genes) > 0) {
      plot_df_genes <- plot_df[which(plot_df[,1] %in% genes),]
      ply <- ply %>%
        add_annotations(x = plot_df_genes[,set1],
                        y = plot_df_genes[,set2],
                        text = plot_df_genes[,1],
                        showarrow = TRUE,
                        arrowhead = 4,
                        arrowsize = .5)
    }
    ply <- ply %>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          height = 800,
          width = 1000,
          filename = paste0(set1,"_",set2,"_Avg_Expression_",Sys.Date())
        )
      )
    ply
    
  })
  
  output$AvgExpr_df <- DT::renderDataTable({
    
    df <- AvgExpr_Plot_df()[,-4]
    datatable(df,
              options = list(lengthMenu = c(5,10, 20, 100, 1000),
                             pageLength = 10,
                             scrollX = T),
              rownames = F)
  })
  output$dnldAvgExpr_df <- downloadHandler(
    filename = function() {
      set1 <- input$AvgExprData1
      set2 <- input$AvgExprData2
      paste0(set1,"_",set2,"_Avg_Expression_",Sys.Date(),".txt")
    },
    content = function(file) {
      df <- AvgExpr_Plot_df()
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

