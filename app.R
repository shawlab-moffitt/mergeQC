type_id <- paste0("v1.0.20240711")

library("shiny")
library("shinycssloaders")
library("shinyjqui")
library("svglite")
library("ggplot2")
library("dplyr")
library("reshape2")
library("plotly")
library("DT")
library("data.table")
library("RecordLinkage")
library("janitor")

options(shiny.maxRequestSize = 10000*1024^2)

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab, na.rm = T)]
}

filter_expression_matrix <- function(data,criteria,proportion) {
  # Extract gene symbols
  gene_symbols <- data[, 1]
  # Extract the numerical data part of the matrix
  numeric_data <- data[, -1]
  # Calculate the number of columns with numerical data
  num_cols <- ncol(numeric_data)
  # Determine the threshold for the number of columns with values < 1
  threshold <- (proportion/100) * num_cols
  # Filter the rows based on the threshold
  keep_rows <- apply(numeric_data, 1, function(row) {
    sum(row < criteria) <= threshold
  })
  # Combine the gene symbols with the filtered numeric data
  filtered_matrix <- cbind(gene_symbols[keep_rows], numeric_data[keep_rows, ])
  colnames(filtered_matrix)[1] <- colnames(data)[1]
  return(filtered_matrix)
}
shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
  }
  inputs
}

# UI ---------------------------------------------------------------------------

# Define UI for application that draws a histogram

ui <- navbarPage("{ MergeQC }",
                 ## Data Preview -----------------------------------------------
                 tabPanel("Data Preview",
                          sidebarLayout(
                            sidebarPanel(
                              h3("File Upload"),
                              fileInput("MatFileInput", label = "Matrix Files:",multiple = T, placeholder = "CTRL + left click or Command + left click files"),
                              fileInput("MetFileInput", label = "Meta Files:",multiple = T, placeholder = "CTRL + left click or Command + left click files"),
                              fluidRow(
                                column(8, style = "padding-right:2px",
                                       uiOutput("rendCleanSampleNames"),
                                       uiOutput("rendCleanMetaNames")
                                       ),
                                column(4, style = "margin-top:10px;padding-left:2px",
                                       tags$a(href="https://github.com/shawlab-moffitt/mergeQC/raw/main/Example_Data/Example_Data.zip", "Download example data", target='_blank')
                                       )
                              )
                              
                            ),
                            
                            # Show a plot of the generated distribution
                            mainPanel(
                              verbatimTextOutput("FileCheckAlerts"),
                              fluidRow(
                                column(2,
                                       uiOutput("rendDnldMatrix")
                                ),
                                column(2,
                                       uiOutput("rendDnldMeta")
                                ),
                                column(2,
                                       uiOutput("rendDnldMatrixFull")
                                ),
                                column(2,
                                       uiOutput("rendDnldMetaFull")
                                ),
                                column(2,
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
                                                          fluidRow(
                                                            column(6,
                                                                   selectizeInput("LogDataOpts","Select datasets to log:", choices = NULL, selected = NULL, multiple = T,
                                                                                  options = list(placeholder = 'Select datasets to log'))
                                                                   ),
                                                            column(3,
                                                                   selectInput("LogMthod","Log Base:", choices = c("Log2","Log10","Log"))
                                                                   ),
                                                            column(3,
                                                                   numericInput("LogPseudo","Pseudocount:", value = 1)
                                                                   )
                                                          ),
                                                          fluidRow(
                                                            column(6,
                                                                   selectizeInput("ExpDataOpts","Select datasets to exponentiate:", choices = NULL, selected = NULL, multiple = T,
                                                                                  options = list(placeholder = 'Select datasets to exponentiate'))
                                                            ),
                                                            column(6,
                                                                   numericInput("ExpNum","Exponential Base:", value = 2)
                                                            )
                                                          ),
                                                          div(DT::dataTableOutput("DataQCtab"), style = "font-size:12px"),
                                                          downloadButton("dnldDataQCtab","Download Table"),
                                                          conditionalPanel(condition = "output.MatFileInputUploaded == true",
                                                                           hr(),
                                                                           h4("Filter out features with expression value of __ in __% of samples"),
                                                                           fluidRow(
                                                                             column(6,
                                                                                    numericInput("matFilterNum","Filter Value:",NULL)
                                                                             ),
                                                                             column(6,
                                                                                    numericInput("matFilterProp","Filter Proportion (%):",NULL)
                                                                             )
                                                                           )
                                                          )
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
                                                          selectizeInput("AvgExprGene","Highlight Genes:",choices = NULL, selected = NULL, multiple = T,
                                                                         options = list(placeholder = 'Select gene(s) to highlight')),
                                                          fluidRow(
                                                            column(4, style = "padding-right:2px;margin-top:10px",
                                                                   checkboxInput("ShowSDresidLines","Show SD Residual Lines", value = T),
                                                                   ),
                                                            column(3, style = "padding-right:4px",
                                                                   numericInput("SDresidLine","+/- SD of Residuals", value = 3, min = 0, step = 1)
                                                                   ),
                                                            column(5, style = "margin-top:25px",
                                                                   actionButton("RemoveOutGenes","Remove genes outside SD")
                                                            )
                                                          ),
                                                          selectizeInput("RemovedGenesDF","Remove Genes:",choices = NULL, selected = NULL, multiple = T),
                                                          #uiOutput("rendRemovedGenesDF"),
                                                          tags$head(
                                                            tags$style(HTML(
                                                              '.selectize-input {
                                                              max-height: 302px;
                                                              overflow-y: auto;
                                                              }'
                                                              )
                                                              )
                                                          )
                                                          ),
                                                 tabPanel("Figure Settings",
                                                 )
                                               )
                                               )
                            ),
                            mainPanel(
                              verbatimTextOutput("FileCheckAlerts2"),
                              fluidRow(
                                column(2,
                                       uiOutput("rendDnldMatrix2")
                                ),
                                column(2,
                                       uiOutput("rendDnldMeta2")
                                ),
                                column(2,
                                       uiOutput("rendDnldMatrixFull2")
                                ),
                                column(2,
                                       uiOutput("rendDnldMetaFull2")
                                ),
                                column(2,
                                       uiOutput("renddownload_notes2")
                                )
                              ),
                              p(),
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


# Server -----------------------------------------------------------------------

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  output$rendCleanSampleNames <- renderUI({
    req(input$MatFileInput)
    checkboxInput("CleanSampleNames","Clean Sample Names", value = FALSE)
  })
  output$rendCleanMetaNames <- renderUI({
    req(input$MetFileInput)
    checkboxInput("CleanMetaNames","Clean Meta Column Names", value = FALSE)
  })
  
  output$MatFileInputUploaded <- reactive({
    return(!is.null(matrix_files_react()))
  })
  outputOptions(output, 'MatFileInputUploaded', suspendWhenHidden=FALSE)
  
  FileCheckAlerts_react <- reactiveVal()
  GenesToRemove <- reactiveVal(NULL)
  
  output$FileCheckAlerts <- renderPrint({
    
    req(FileCheckAlerts_react())
    text <- paste(FileCheckAlerts_react(), collapse = "\n")
    cat(text)
    
  })
  output$FileCheckAlerts2 <- renderPrint({
    
    req(FileCheckAlerts_react())
    text <- paste(FileCheckAlerts_react(), collapse = "\n")
    cat(text)
    
  })
  output$renddownload_notes <- renderUI({
    if (length(FileCheckAlerts_react()) > 0) {
      downloadButton("download_notes","Download data processing notes")
    }
  })
  output$renddownload_notes2 <- renderUI({
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
  output$download_notes2 <- downloadHandler(
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
  
  MatFileIn <- reactive({
    
    if (isTruthy(matrix_files_react())) {
      TRUE
    } else { FALSE }
    
  })
  
  observeEvent(input$MatFileInput,{
    
    #req(matrix_files_react())
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
    updateSelectizeInput(session,"ExpDataOpts",choices = tools::file_path_sans_ext(matrix_files_df_react()[,"name"]), selected = NULL)
  })
  
  matrix_files_in_react <- reactive({
    
    req(matrix_files_df_react())
    file_df <- matrix_files_df_react()
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
        df_list[[tabName]] <- df
      }
      df_list
    }
    
  })
  
  
  
  matrix_files_react <- reactive({
    
    req(matrix_files_in_react())
    df_list <- matrix_files_in_react()
    logData <- input$LogDataOpts
    logMethod <- input$LogMthod
    logPseudo <- input$LogPseudo
    expData <- input$ExpDataOpts
    ExpNum <- input$ExpNum
    FiltNum <- input$matFilterNum
    FiltProp <- input$matFilterProp
    for (tabName in names(df_list)) {
      df <- df_list[[tabName]]
      if (tabName %in% logData) {
        if (logMethod == "Log2") {
          df[,-1] <- log2(as.matrix(df[,-1]+logPseudo))
        } else if (logMethod == "Log10") {
          df[,-1] <- log10(as.matrix(df[,-1]+logPseudo))
        } else if (logMethod == "Log") {
          df[,-1] <- log(as.matrix(df[,-1]+logPseudo))
        }
      }
      if (tabName %in% expData) {
        df[,-1] <- ExpNum^(as.matrix(df[,-1]))
      }
      if (isTruthy(FiltNum) & isTruthy(FiltProp)) {
        if (FiltNum > 0 & FiltProp > 0) {
          df <- filter_expression_matrix(df,FiltNum,FiltProp)
        }
      }
      df_list[[tabName]] <- df
    }
    df_list
    
  })
  
  
  observeEvent(input$MatFileInput, {
    
    #req(matrix_files_react())
    appendTab(inputId = "tabs",
              tabPanel("Merged Meta",
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
    #req(input$CleanMetaNames)
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
  output$rendDnldMatrix2 <- renderUI({
    req(matrix_merged_inn())
    downloadButton("dnldMatrix2","Inner Merged Matrix")
  })
  
  output$rendDnldMatrixFull <- renderUI({
    req(matrix_merged_full())
    downloadButton("dnldMatrixFull","Full Merged Matrix")
  })
  output$rendDnldMatrixFull2 <- renderUI({
    req(matrix_merged_full())
    downloadButton("dnldMatrixFull2","Full Merged Matrix")
  })
  
  matrix_merged_inn <- reactive({
    
    df_list <- matrix_files_react()
    file_df <- matrix_files_df_react()
    
    if (!is.null(input$CleanSampleNames)) {
      if (input$CleanSampleNames) {
        df_list <- lapply(df_list,function(x) {
          x <- x %>% janitor::clean_names(case = "none")
          return(x)
        })
      }
    }
    
    masterDF <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2),
                       df_list)
    
    if (isTruthy(input$RemoveOutGenes)) {
      if (length(input$RemovedGenesDF) > 0) {
        masterDF <- masterDF[which(!masterDF[,1] %in% input$RemovedGenesDF),]
      }
    }
    
    masterDF
    
  })
  
  matrix_merged_full <- reactive({
    
    df_list <- matrix_files_react()
    file_df <- matrix_files_df_react()
    
    if (!is.null(input$CleanSampleNames)) {
      if (input$CleanSampleNames) {
        df_list <- lapply(df_list,function(x) {
          x <- x %>% janitor::clean_names(case = "none")
          return(x)
        })
      }
    }
    
    masterDF <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = TRUE),
                       df_list)
    
    if (isTruthy(input$RemoveOutGenes)) {
      if (length(input$RemovedGenesDF) > 0) {
        masterDF <- masterDF[which(!masterDF[,1] %in% input$RemovedGenesDF),]
      }
    }
    
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
  output$dnldMatrixFull2 <- downloadHandler(
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
  output$dnldMatrix2 <- downloadHandler(
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
  output$rendDnldMeta2 <- renderUI({
    req(meta_merged())
    downloadButton("dnldMeta2","Inner Merged Meta")
  })
  
  meta_merged <- reactive({
    
    if (isTruthy(meta_files_df_react())) {
      req(matrix_merged_inn())
      df_list <- meta_files_react()
      file_df <- meta_files_df_react()
      matrix_merged_inn <- matrix_merged_inn()
      
      if (!is.null(input$CleanMetaNames)) {
        if (input$CleanMetaNames) {
          df_list <- lapply(df_list,function(x) {
            x <- x %>% janitor::clean_names(case = "none")
            return(x)
          })
        }
      }
      
      if (!is.null(input$CleanSampleNames)) {
        if (input$CleanSampleNames) {
          df_list <- lapply(df_list,function(x) {
            x <- x %>% janitor::clean_names(case = "none")
            return(x)
          })
        }
      }
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
    } else {
      req(matrix_merged_inn())
      matrix_merged_inn <- matrix_merged_inn()
      df_list <- matrix_files_react()
      merged_meta <- data.frame(SampleName = colnames(matrix_merged_inn)[-1],
                                File_Name = unname(unlist(sapply(names(df_list),function(x) {
                                  rep(x,ncol(df_list[[x]])-1)
                                }))))
      merged_meta
    }
    
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
  output$dnldMeta2 <- downloadHandler(
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
    
    FileCheckAlerts_list <- c()
    if (isTruthy(input$MatFileInput) & isTruthy(input$MetFileInput)) {
      SampsDiffInMet <- setdiff(meta_merged()[,1],intersect(meta_merged()[,1], colnames(matrix_merged_full())[-1]))
      SampsDiffInMat <- setdiff(colnames(matrix_merged_full())[-1],intersect(meta_merged()[,1], colnames(matrix_merged_full())[-1]))
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
    }
    if (isTruthy(input$MatFileInput) & !isTruthy(input$MetFileInput)) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                paste0("Number of matrix files uploaded: ",length(matrix_files_react())),
                                paste0("Total number of samples from matrix data: ",ncol(matrix_merged_full()[,-1])),
                                paste0("Total number of features consistent in all matrix data: ",nrow(matrix_merged_inn())))
    }
    
    if (length(input$LogDataOpts) > 0) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,paste0("Datasets transformed by ",input$LogMthod,"+",input$LogPseudo," :"))
      FileCheckAlerts_list <- c(FileCheckAlerts_list,paste0("  ",input$LogDataOpts))
    }
    if (length(input$ExpDataOpts) > 0) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,paste0("Datasets transformed exponentially with base",input$ExpNum," :"))
      FileCheckAlerts_list <- c(FileCheckAlerts_list,paste0("  ",input$ExpDataOpts))
    }
    if (isTruthy(input$matFilterNum) & isTruthy(input$matFilterProp)) {
      if ((input$matFilterNum > 0) & (input$matFilterProp > 0)) {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Features with an expression level less than ",input$matFilterNum," in ",input$matFilterProp,"% of samples filtered out"))
      }
    }
    if (isTruthy(input$RemoveOutGenes)) {
      if (length(input$RemovedGenesDF) > 0) {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste0("When comparing average expression of ",input$AvgExprData1," and ",input$AvgExprData2,", genes that were expressed\n",
                                         "  outside ",input$SDresidLine," standard deviations of the residual were removed from the merged data.\n  Number of genes removed: ",length(input$RemovedGenesDF)))
      }
    }
    
    FileCheckAlerts_react(FileCheckAlerts_list)
    
  })
  
  # QC stats -------------------------------------------------------------------
  
  QC_react <- reactive({
    
    df_list <- matrix_files_react()
    
    df_list_qc <- lapply(df_list,function(x){
      
      if (isTruthy(AvgExpr_Plot_df())) {
        if (isTruthy(input$RemoveOutGenes)) {
          if (length(input$RemovedGenesDF) > 0){
          removeGenes <- input$RemovedGenesDF
          #removeGenes <- AvgExpr_Plot_df()[which(AvgExpr_Plot_df()[,"Color_Col"] == "gray84"),1]
          x <- x[which(!x[,1] %in% removeGenes),]
          }
        }
      }
      
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
      x_mean <- rowMeans(as.matrix(df[,-1]))
      x_mean_df <- data.frame(Gene = df[,1],
                              Avg_Expression = x_mean)
      colnames(x_mean_df)[2] <- df_name
      return(x_mean_df)
    })
    
    df_list_avgExpr_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = T),
                                 df_list_avgExpr)
    
    
    #if (input$RemoveOutGenes) {
      #if (length(GenesToRemove()) > 0) {
      #  df_list_avgExpr_df <- df_list_avgExpr_df[which(!df_list_avgExpr_df[,1] %in% GenesToRemove()),]
      #}
    #}
    
    df_list_avgExpr_df
  })
  
  output$AvgExprDens <- DT::renderDataTable({
    
    df <- DensityPlot_df()
    if (isTruthy(input$RemoveOutGenes)) {
      if (length(input$RemovedGenesDF) > 0) {
        df <- df[which(!df[,1] %in% input$RemovedGenesDF),]
      }
    }
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
      if (length(GenesToRemove()) > 0) {
        df <- df[which(!df[,1] %in% GenesToRemove()),]
      }
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  
  DensityPlot_react <- reactive({
    
    if (length(input$DataQCtab_rows_selected) > 0) {
      df_list_avgExpr_df <- DensityPlot_df()
      if (isTruthy(input$RemoveOutGenes)) {
        if (length(input$RemovedGenesDF) > 0) {
          df_list_avgExpr_df <- df_list_avgExpr_df[which(!df_list_avgExpr_df[,1] %in% input$RemovedGenesDF),]
        }
      }
      files_select <- QC_react()[input$DataQCtab_rows_selected,1]
      xMax <- max(QC_react()[input$DataQCtab_rows_selected,6])
      df_list_avgExpr_df_select <- df_list_avgExpr_df[,c(colnames(df_list_avgExpr_df)[1],files_select)]
      df_list_avgExpr_df_select_melt <- melt(data.table(df_list_avgExpr_df_select), id.vars = "Gene")
      
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
    SDfromResid <- input$SDresidLine
    
    
    if (length(df[,set1]) > 0 & length(df[,set2]) > 0) {
      df <- df[,c(colnames(df)[1],set1,set2)]
      
      df <- df[complete.cases(df),]
      
      reg = lm(df[,set2] ~ df[,set1])
      R2 = summary(reg)$r.squared
      
      df[,"Regression Fit"] <- reg %>% predict(df[,set1, drop = F])
      
      #save(list = ls(), file = "Avg_Expr_Resid.RData", envir = environment())
      
      df$Residuals <- residuals(reg, na.rm = T)
      sd_residuals <- sd(df$Residuals, na.rm = T)
      
      df[,"Upper Residual Line"] <- df[,"Regression Fit"] + SDfromResid * sd_residuals
      df[,"Lower Residual Line"] <- df[,"Regression Fit"] - SDfromResid * sd_residuals
      
      
      df[,"Outside Upper Residual Line"] <- ifelse(df[,set2] > df[,"Upper Residual Line"],
                                                  TRUE,FALSE)
      df[,"Outside Lower Residual Line"] <- ifelse(df[,set2] < df[,"Lower Residual Line"],
                                                  TRUE,FALSE)
      
      df[,"Color_Col"] <- ifelse(df[,"Outside Upper Residual Line"] == TRUE | df[,"Outside Lower Residual Line"] == TRUE,"gray84","cadetblue3")
      
    }
    
    df
    
    
  })
  
  #RemoveGeneDF <- reactiveVal()
  
  observeEvent(input$RemoveOutGenes, {
    df <- AvgExpr_Plot_df()
    removeGenes <- c(GenesToRemove())
    removeGenes <- c(removeGenes,AvgExpr_Plot_df()[which(AvgExpr_Plot_df()[,"Color_Col"] == "gray84"),1])
    removeGenes <- unique(removeGenes)
    GenesToRemove(removeGenes)
  })
  
  observe({
    #GenesWereRemoved(input$RemovedGenesDF)
    updateSelectizeInput(session,"RemovedGenesDF", choices = DensityPlot_df()[,1], selected = GenesToRemove(), server = T)
  })
  #observeEvent(input$RemoveOutGenes, {
  #  output$rendRemovedGenesDF <- renderUI({
  #    selectizeInput("RemovedGenesDF","Genes Removed:",choices = DensityPlot_df()[,1], selected = NULL, multiple = T)
  #  })
  #})
  
  
  
  
  #js <- c(
  #  "$('[id^=checkb]').on('click', function(){",
  #  "  var id = this.getAttribute('id');",
  #  "  var i = parseInt(/checkb(\\d+)/.exec(id)[1]);",
  #  "  var value = $(this).prop('checked');",
  #  "  var info = [{row: i, col: 2, value: value}];",
  #  "  Shiny.setInputValue('RemovedGenesDF_cell_edit:DT.cellInfo', info);",
  #  "})"
  #)
  #output$RemovedGenesDF <- DT::renderDataTable({
  #  
  #  df <- RemoveGeneDF()
  #  datatable(df,
  #            options = list(paging = F),
  #            escape = FALSE,
  #            editable = list(target = "cell", disable = list(columns = 2)),
  #            selection = "none",
  #            callback = JS(js),
  #            rownames = F)
  #  
  #})
  #
  #observeEvent(input[["RemovedGenesDF_cell_edit"]], { 
  #  info <- input[["RemovedGenesDF_cell_edit"]] # this input contains the info of the edit
  #  print(head(info))
  #  RemoveGeneDF(editData(RemoveGeneDF(), info))
  #})
  #
  #observe({
  #  
  #  print(head(input[["RemovedGenesDF_cell_edit"]]))
  #  
  #})
  
  AvgExpr_Plot_react <- reactive({
    
    req(AvgExpr_Plot_df())
    
    set1 <- input$AvgExprData1
    set2 <- input$AvgExprData2
    plot_df1 <- AvgExpr_Plot_df()
    genes <- input$AvgExprGene
    RemoveGenes <- input$RemovedGenesDF
    plot_df <- plot_df1
    
    
    #if (input$RemoveOutGenes) {
      if (length(RemoveGenes) > 0) {
        plot_df <- plot_df[which(!plot_df[,1] %in% RemoveGenes),]
      }
    #}
    
    
    p <- ggplot(plot_df, aes(x = !!sym(set1), y = !!sym(set2),
                             color = Color_Col,
                             text = paste0("</br>Gene: ",Gene,
                                           "</br>",set1,": ",round(!!sym(set1),3),
                                           "</br>",set2,": ",round(!!sym(set2),3)))) +
      geom_point(size = 1) +
      #geom_point(size = 1,color = "cadetblue") +
      theme_minimal() +
      theme(legend.position="none",
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 20),
            plot.margin = margin(1, 0, 0, 0, "cm")) +
      xlab(paste0("Average Expression: ",set1))+
      ylab(paste0("Average Expression: ",set2)) +
      labs(title = paste0("Average Expression: ",set1," vs. ",set2)) +
      scale_color_manual(values = c("cadetblue3","gray84"))
    
    
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
    ScatterTitle_in <- paste0("Average Expression: ",set1," vs. ",set2)
    showRSD <- input$ShowSDresidLines
    SDfromResid <- input$SDresidLine
    ply <- ggplotly(p, tooltip = "text")
    genes <- input$AvgExprGene
    if (length(genes) > 0) {
      plot_df_genes <- plot_df[which(plot_df[,1] %in% genes),]
      ply <- ply %>%
        add_annotations(x = plot_df_genes[,set1],
                        y = plot_df_genes[,set2],
                        text = paste0("<b>",plot_df_genes[,1],"</b>"),
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
    if (length(plot_df[,set1]) > 0 & length(plot_df[,set2]) > 0) {
      reg = lm(plot_df[,set2] ~ plot_df[,set1])
      R2 = summary(reg)$r.squared
      # Add title and subtitle
      coef <- paste0("Coefficients: y = ",round(coefficients(reg)[2],3),"x"," + ",round(coefficients(reg)[1],3))
      rSqu <- paste0("R-Squared: ",round(R2,3))
      ply <- ply %>% layout(title = list(text = paste0(ScatterTitle_in,
                                                   '<br>',
                                                   '<sup>',
                                                   coef,
                                                   '<br>',
                                                   rSqu,
                                                   '</sup>'),
                                     x = 0,
                                     xref='paper',
                                     yref='paper',
                                     align = "left"
      )
      )
      colnames(plot_df)[which(colnames(plot_df) %in% c(set1,set2))] <- c("set1","set2")
      ply <- ply %>%
        add_lines(data = plot_df, x = ~set1, y = ~`Regression Fit`, name = "Regression Fit",
                  line = list(color = "black", width=2, dash="dash"))
      if (showRSD) {
        ply <- ply %>%
          add_lines(data = plot_df, x = ~set1, y = ~`Upper Residual Line`, name = paste("Residual +",SDfromResid,"SD"),
                    line = list(color = "red", width=2, dash="dash")) %>%
          add_lines(data = plot_df, x = ~set1, y = ~`Lower Residual Line`, name = paste("Residual -",SDfromResid,"SD"),
                    line = list(color = "red", width=2, dash="dash"))
      }
      ply
    }

  })
  
  output$AvgExpr_df <- DT::renderDataTable({
    
    df <- AvgExpr_Plot_df()
    if (isTruthy(input$RemoveOutGenes)) {
      if (length(input$RemovedGenesDF) > 0) {
        df <- df[which(!df[,1] %in% input$RemovedGenesDF),]
      }
    }
    df <- df[,-ncol(AvgExpr_Plot_df())]
    datatable(df,
              options = list(lengthMenu = c(5,10, 20, 100, 1000),
                             pageLength = 10,
                             scrollX = T),
              rownames = F) %>%
      formatRound(columns = colnames(df)[c(2:7)], digits = 4)
  })
  output$dnldAvgExpr_df <- downloadHandler(
    filename = function() {
      set1 <- input$AvgExprData1
      set2 <- input$AvgExprData2
      paste0(set1,"_",set2,"_Avg_Expression_",Sys.Date(),".txt")
    },
    content = function(file) {
      df <- AvgExpr_Plot_df()[,-ncol(AvgExpr_Plot_df())]
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

