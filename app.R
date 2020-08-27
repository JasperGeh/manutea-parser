# Libraries ####
rm(list = ls())
library(shiny)
library(dplyr) # for quickly adding split genes as column to df
library(container) # includes dict data structure
library(data.table) # for fast transposing of dfs
library(stringr) # for string replacement


# UI ####
ui <- fluidPage(
  h1(style="text-align: center;", "Manutea's GO Parser"),
  p(style="text-align: center;", "Hey Manutea! Load the CSV and switch to the 'Gene Dictionary'-tab to select a gene and see the sorted ontologies. The 'Compiled Table'-tab contains the top-3 ontologies for each term. Let me know if you need more or something doesn't work. —Jasper, 2020-08-27"),
  fluidRow(
    column(10, 
           align='center', 
           offset = 1, 
           fileInput('datafile', 'Choose CSV file',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')))),
  tabsetPanel(
    tabPanel("Input Table",
             br(),
             dataTableOutput(outputId = 'inTab')),
    tabPanel("Gene Dictionary",
             br(),
             wellPanel(uiOutput("gene")),
             dataTableOutput(outputId = 'dictTab')),
    tabPanel("Compiled Table",
             br(),
             fluidRow(
               column(
                 10, 
                 align='center', 
                 offset=1, 
                 downloadButton("dlComp", "Download Compiled Table"))),
             dataTableOutput(outputId = 'resTab'))
    )
)


# Server #####
server <- function(input, output) {

  filedata <- reactive({
    infile <- input$datafile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    df <- read.csv(infile$datapath, stringsAsFactors = FALSE)
    df <- df %>%
      mutate(genes = str_replace_all(genes, "\\|", ", ")) %>%
      mutate(genes = strsplit(genes, ", ")) # create new column Genelist which splits Gene-column at commas
    df
  })
  dictio <- reactive({
    df <- filedata()
    genedict <- dict() # Create empty dict from container library
    
    for (i in 1:nrow(df)) { # Iterate through rows of the DF
      l = df[[i,3]] # save list of split genes as "l"
      for (e in l) { # Iterate through list of genes
        if (!genedict$has(e)) { # check if dict already has gene in it and if not
          new_entry <- list(c(df[[i,1]], df[[i,4]], df[[i,2]]))
          genedict$set(e, new_entry, add=TRUE) # add GO with p-Value to dict as list
        } else { # check if pValue of current row is smaller than previously saved and if yes
          curr <- genedict$peek(e) # save current dict content for that gene
          new_entry <- append(curr, list(c(df[[i,1]], df[[i,4]], df[[i,2]])))
          genedict$set(e, new_entry, add=TRUE) # append current GO with p-Value to previous content and save
        }
      }
    }
    
    results <- list() # create empty list 
    for (k in genedict$keys()) { # iterate through genes in dict
      dfk = transpose(data.frame(genedict$get(k))) # extract dict tuples into transposed df
      dfk$V3 = as.numeric(dfk$V3) # convert results to num for ordering
      names(dfk) = c("Gene Ontology", "Term", "FDR-Value")
      results = append(results, list(dfk[order(dfk$"FDR-Value")[1:nrow(dfk)], ])) # [1:n] to save rows 1-n of ordered p-Values
    }
    names(results) = genedict$keys() # Name dfs in result list
    results
  })
  
  res <- reactive({
    results <- dictio()
    compTab <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c("Gene", "GO1", "Term1", "FDR1", "GO2", "Term2", "FDR2", "GO3", "Term3", "FDR3"))
    for (e in 1:length(results)) {
      sub <- results[[e]][1:3,]
      compTab[e,] <- c(names(results[e]), sub[1,1], sub[1,2], sub[1,3], sub[2,1], sub[2,2], sub[2,3], sub[3,1], sub[3,2], sub[3,3])
    }
    compTab
  })
  
  
  
  
  
  #The following function output content
  output$gene <- renderUI({
    df <- dictio()
    if (is.null(df)) return(NULL)
    
    items=sort(names(df))
    names(items)=items
    selectInput("gene", "Gene", items)
  })
  
  output$inTab <- renderDataTable({
    filedata()
  })
  
  output$dictTab <- renderDataTable({
    dictio()[[input$gene]]
  })
  
  output$resTab <- renderDataTable({
    res()
  })
  
  output$dlComp <- downloadHandler(
    filename = "CompiledTable.csv",
    content = function(file) {
      write.csv(res(), file, row.names = FALSE)
    }
  )
}



shinyApp(ui = ui, server = server)
