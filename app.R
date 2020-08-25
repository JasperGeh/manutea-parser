rm(list = ls())
library(shiny)
library(dplyr) # for quickly adding split genes as column to df
library(container) # includes dict data structure
library(data.table) # for fast transposing of dfs


ui <- fluidPage(
  numericInput(inputId = "Gene", 
              label = "Select Gene", 
              value = 1, min = 1, max = 5),
  tableOutput(outputId = "tabl")
)


server <- function(input, output) {
  df <- data.frame(Ontology = c("GO1", "GO2", "GO3", "GO4", "GO5", "GO6", "GO7"), 
                   Genes = c("A, B", "A, D", "C", "A, D, F", "D, C, F", "C, F", "A, B, F"), 
                   pValue = c(0.001, 0.06, 0.000002, 0.5, 0.00015, 0.0004, 0.01), stringsAsFactors=FALSE)
  
  df <- df %>%
    mutate(Genelist = strsplit(Genes, ", ")) # create new column Genelist which splits Gene-column at commas
  
  genedict <- dict() # Create empty dict from container library
  
  for (i in 1:nrow(df)) { # Iterate through rows of the DF
    l = df[[i,4]] # save list of split genes as "l"
    for (e in l) { # Iterate through list of genes
      if (!genedict$has(e)) { # check if dict already has gene in it and if not
        genedict$set(e, list(c(df[[i,1]], df[[i,3]])), add=TRUE) # add GO with p-Value to dict as list
      } else { # check if pValue of current row is smaller than previously saved and if yes
        curr <- genedict$peek(e) # save current dict content for that gene
        genedict$set(e, append(curr, list(c(df[[i,1]], df[[i,3]]))), add=TRUE) # append current GO with p-Value to previous content and save
      }
    }
  }
  
  results <- list() # create empty list 
  for (k in genedict$keys()) { # iterate through genes in dict
    dfk = transpose(data.frame(genedict$get(k))) # extract dict tuples into transposed df
    dfk$V2 = as.numeric(dfk$V2) # convert results to num for ordering
    results = append(results, list(dfk[order(dfk$V2)[1:3], ])) # save rows 1-3 of ordered p-Values
  }
  names(results) = genedict$keys() # Name dfs in result list
  rm(list = c("e", "i", "k", "l", "curr", "dfk", "genedict")) # Remove intermediate variables
  
  
  output$tabl <- renderDataTable({
    results[[input$Gene]]
  })
}


shinyApp(ui = ui, server = server)

