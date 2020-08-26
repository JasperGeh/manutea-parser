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
  p(style="text-align: center;", "Hey Manutea! Load the CSV, hit 'Analyse', and switch to the 'Analysis'-tab to select a gene and see the sorted ontologies.
    Let me know if you need more or something doesn't work. —Jasper"),
  fluidRow(
    column(10, 
           align='center', 
           offset = 1, 
           fileInput('datafile', 'Choose CSV file',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
           actionButton("go", "Analyse"))),
  tabsetPanel(
    tabPanel("Input Table",
             dataTableOutput(outputId = 'inTab')),
    tabPanel("Analysis",
             br(),
             wellPanel(uiOutput("gene")),
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
  res <- eventReactive(input$go, {
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
  
  #The following set of functions populate the column selectors
  output$gene <- renderUI({
    df <-res()
    if (is.null(df)) return(NULL)
    
    items=sort(names(df))
    names(items)=items
    selectInput("gene", "Gene", items)
  })
  
  output$inTab <- renderDataTable({
    filedata()
  })
  
  output$resTab <- renderDataTable({
    res()[[input$gene]]
  })
}



shinyApp(ui = ui, server = server)
