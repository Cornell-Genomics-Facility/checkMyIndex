library(shiny)

options(shiny.sanitize.errors = FALSE,   # to display informative error messages
        shiny.maxRequestSize = 5*1024^2) # limit size for the input file to upload (5Mo here)

shinyServer(function(input, output, session) {
  
  # reset input parameters when pressing the reset button
  observeEvent(input$reset, {
    shinyjs::reset("allParameters")
    # previous inputFile disapears but is still sent to the server
    # known issue: https://github.com/daattali/shinyjs/issues/104
    rv$inputFile <- NULL
    rv$inputFile2 <- NULL
    rv$testdata <- "none"
  })
  # automatically go to the first tab when pressing "reset" or providing any input file
  observeEvent(c(input$inputFile, input$inputFile2, input$testdata, input$reset), ignoreInit=TRUE, {
    updateTabsetPanel(session, "mainPanel", selected = "inputIndexes")
    shinyjs::hide("proposedSolution")
    shinyjs::hide("visualization")
    shinyjs::hide("colorBalancing")
  })
  # hide solution and heatmap when changing one of these parameters
  observeEvent(c(input$nbSamples, input$multiplexingRate, input$chemistry), ignoreInit=TRUE, {
    shinyjs::hide("proposedSolution")
    shinyjs::hide("visualization")
    shinyjs::hide("colorBalancing")
  })
  # piece of code to fix the input file reset problem (trick: pass it into a reactive value)
  rv <- reactiveValues(inputFile=NULL, inputFile2=NULL, testdata="none")
  observeEvent(input$inputFile, {rv$inputFile <- input$inputFile})
  observeEvent(input$inputFile2, {rv$inputFile2 <- input$inputFile2})
  observeEvent(input$testdata, {rv$testdata <- input$testdata})
  # get input files
  fileInput <- reactive({rv$inputFile})
  fileInput2 <- reactive({rv$inputFile2})
  testData <- reactive({rv$testdata})
  # tell UI if inputFiles are present
  output$inputFileProvided <- reactive({!is.null(rv$inputFile)})
  output$inputFile2Provided <- reactive({!is.null(rv$inputFile2)})
  output$testdataProvided <- reactive({rv$testdata != "none"})
  outputOptions(output, "inputFileProvided", suspendWhenHidden=FALSE)
  outputOptions(output, "inputFile2Provided", suspendWhenHidden=FALSE)
  outputOptions(output, "testdataProvided", suspendWhenHidden=FALSE)
  # automatically go to the proposed solution when pressing "search for a solution" and show solution and heatmap
  observeEvent(input$go, {
    if (is.null(tryCatch({displaySolution()$solution}, error = function(e) NULL))){
      shinyjs::show("proposedSolution")
      shinyjs::show("visualization")
      shinyjs::show("colorBalancing")
    }
    updateTabsetPanel(session, "mainPanel", selected = "proposedSolution")
  })
  
  # list of input index 1
  inputIndex <- reactive({
    # print('Input index 1 entered')
    if (testData() %in% c("simple", "dual")){
      file <- ifelse(testData() == "simple", "www/inputIndexesExample.txt", "www/index24-i7.txt")
    } else{
      if (!is.null(fileInput())) file <- fileInput()$datapath else return(NULL)
    }
    index <- tryCatch({readIndexesFileWithWeights(file)$index}, 
                      error = function(e) stop("An error occured when loading index 1 file, please check its structure."))
    index <- addColors(index, input$chemistry)
    # print(paste0('index$color: ', index$color))
    index$score <- scores(index$sequence)
    return(index)
  })
  output$indexUploaded <- reactive({!is.null(inputIndex())})
  outputOptions(output, "indexUploaded", suspendWhenHidden=FALSE)
  output$inputIndex <- DT::renderDT({inputIndex()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  output$textIndex <- renderText({tryCatch({
    index <- inputIndex()
    if (is.null(index)){
      "No index file loaded yet, use the left panel to select an input file."
    } else{
      paste0("The table below shows the ", nrow(index), " indices 1 (i7) uploaded with the colors corresponding
             to the chosen Illumina chemistry and the minimum number of mismatches with the other indices.
             Note that the smallest number of mismatches between two indices of this list is ", min(index$score), ".")
    }
  }, error = function(e) NULL)})
  
  # list of input indices 2
  inputIndex2 <- reactive({
    # print('Input index 2 entered')
    inputIndex()
    if (testData() == "simple") return(NULL)
    if (testData() == "dual"){
        file2 <- "www/index24-i5.txt"
    } else{
        if (!is.null(fileInput2())){
            if (is.null(inputIndex())) stop("Please load index 1 (i7) first.")
            file2 <- fileInput2()$datapath
        } else {
            if (!is.null(inputIndex())){
                file2 <- fileInput()$datapath
            } else {
                return(NULL)
            }
        }
    }

    # index2 <- tryCatch({readIndexesFileWithWeights(file2)$index}, 
    #                    error = function(e) stop("An error occured when loading index 2 file, please check its structure."))
    if (testData() == "dual" | !is.null(fileInput2())){
       index2 <- tryCatch({readIndexesFileWithWeights(file2)$index}, 
                      error = function(e) stop("An error occured when loading index 2 file, please check its structure."))
    } else {
       result <- tryCatch({readIndexesFileWithWeights(file2)}, 
                        error = function(e) stop("An error occured when loading index 2 file, please check its structure."))
       if (!is.null(result$index2)) {
           index2 <- result$index2
       } else {
           return(NULL)
       }
    }
    index2 <- addColors(index2, input$chemistry)
    # print(paste0('index2$color: ', index2$color))
    index2$score <- scores(index2$sequence)
    return(index2)
  })
  output$indexUploaded2 <- reactive({!is.null(inputIndex2())})
  outputOptions(output, "indexUploaded2", suspendWhenHidden=FALSE)
  output$inputIndex2 <- DT::renderDT({inputIndex2()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  output$textIndex2 <- renderText({tryCatch({
    index2 <- inputIndex2()
    if (is.null(index2)){
      ""
    } else{
      paste0("The table below shows the ", nrow(index2), " index 2 (i5) uploaded with the colors corresponding
             to the chosen Illumina chemistry and the minimum number of mismatches with the other indices.
             Note that the smallest number of mismatches between two indices of this list is ", min(index2$score), ".")
    }
  }, error = function(e) NULL)})
  
  # i7 and i5 pairing
  output$i7i5sameLength <- reactive({
    !is.null(inputIndex()) & !is.null(inputIndex2()) && nrow(inputIndex())==nrow(inputIndex2())
  })
  outputOptions(output, "i7i5sameLength", suspendWhenHidden=FALSE)
  output$i7i5pairing <- renderUI({checkboxInput("i7i5pairing", "i7 and i5 indices pairing (Illumina UDIs)", value=FALSE)})
  output$textPairingTable <- renderText({tryCatch({
    index <- inputIndex()
    index2 <- inputIndex2()
    if (is.null(index) | is.null(index2) || !input$i7i5pairing){
      ""
    } else{
      "The table below shows the pairing between the i7 and i5 input indices"
    }
  }, error = function(e) NULL)})
  output$pairingTable <- DT::renderDT({
    index <- inputIndex()
    index2 <- inputIndex2()
    if (is.null(index) | is.null(index2) || !input$i7i5pairing){
      return(NULL)
    } else{
      out <- cbind(index[, 1:2], index2[, 1:2])
      names(out) <- paste(names(out), c("i7", "i7", "i5", "i5"), sep=".")
      return(out)
    }
  }, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  
  # propose both the possible nb samples and multiplexing rates according to the input list of indices
  output$nbSamples <- renderUI({
    index <- tryCatch({inputIndex()}, error = function(e) NULL)
    if (is.null(index)){
      return("")
    } else{
      nr <- ifelse(input$chemistry == "2", nrow(index[substr(index$sequence, 1, 2) != "GG",]), nrow(index))
      index2 <- tryCatch({inputIndex2()}, error = function(e) NULL)
      if (is.null(index2)){
        nbSamples <- nr
      } else{
        if (input$i7i5pairing){
          nbSamples <- ifelse(input$chemistry == "2",
                              sum(substr(index$sequence, 1, 2) != "GG" & substr(index2$sequence, 1, 2) != "GG"),
                              nr)
        } else{
          nr2 <- ifelse(input$chemistry == "2", nrow(index2[substr(index2$sequence, 1, 2) != "GG",]), nrow(index2))
          nbSamples <- nr*nr2
        }
      }
      numericInput("nbSamples", label="Total number of samples in the experiment", value=nbSamples, min=2, step=1)
    }
  })
  output$multiplexingRate <- renderUI({
    index <- tryCatch({inputIndex()}, error = function(e) NULL)
    if (is.null(index) | is.null(input$nbSamples)){
      return("")
    } else{
      nbSamples <- as.numeric(input$nbSamples)
      if (is.na(nbSamples) || nbSamples %% 1 != 0) stop("Number of samples must be an integer.")
      if (nbSamples <= 1) stop("Number of samples must be greater than 1.")
      mr <- 1:nbSamples
      choices <- mr[sapply(mr, function(x) nbSamples %% x == 0)]
      nr <- ifelse(input$chemistry == "2", nrow(index[substr(index$sequence, 1, 2) != "GG",]), nrow(index))
      index2 <- tryCatch({inputIndex2()}, error = function(e) NULL)
      if (is.null(index2)){
        nr2 <- 1
      } else{
        nr2 <- ifelse(input$chemistry == "2", nrow(index2[substr(index2$sequence, 1, 2) != "GG",]), nrow(index2))
      }
      choices <- choices[choices <= nr*nr2]
      selectInput("multiplexingRate", label="Multiplexing rate (i.e. # of samples per pool)", 
                  choices=choices, selected=choices[length(choices)])
      # selectInput("multiplexingRate", label="Multiplexing rate (i.e. # of samples per pool)", 
      #             choices=choices, selected=ifelse(input$chemistry == "2", choices[length(choices)], choices[2]))
    }
  })
  
  # generate list(s) of indices
  generateList <- reactive({
    index <- inputIndex()
    # print(paste0('column names: ', names(index)))
    # print(paste0('index[, -5]: ', index[, -5]))
    
    return(generateListOfIndexesCombinations(index = index[, -5],
                                             nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                             completeLane = input$completeLane,
                                             selectCompIndexes = input$selectCompIndexes,
                                             chemistry = input$chemistry))
  })
  generateList2 <- reactive({
    index2 <- inputIndex2()
    if (is.null(index2)){
      return(NULL)
    } else{
      return(generateListOfIndexesCombinations(index = index2[, -5],
                                               nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                               completeLane = input$completeLane,
                                               selectCompIndexes = input$selectCompIndexes,
                                               chemistry = input$chemistry))
    }
  })
  generateListPairedIndexes <- reactive({
    index <- inputIndex()
    index2 <- inputIndex2()
    if (is.null(index) & is.null(index2)){
      return(NULL)
    } else{
      if (input$chemistry == "2"){
        mask <- which(substr(index$sequence, 1, 2) == "GG" | substr(index2$sequence, 1, 2) == "GG")
        if (length(mask) > 0){
          index <- index[-mask,]
          index2 <- index2[-mask,]
        }
      }
      names(index2) <- paste0(names(index2), "2")
      index <- cbind(index[, -5], index2[, -5])
      return(generateListOfIndexesCombinations(index = index,
                                               nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                               completeLane = input$completeLane,
                                               selectCompIndexes = input$selectCompIndexes,
                                               chemistry = input$chemistry))
    }
  })
  
  # text describing the solution
  output$textDescribingSolution <- renderText({
    if (!is.null(tryCatch({displaySolution()$solution}, error = function(e) NULL))){
      paste("Below is a solution for", as.numeric(input$nbSamples)/as.numeric(input$multiplexingRate),
            "pool(s) of", input$multiplexingRate, "samples using the parameters specified. The table contains
             one row per sample to be sequenced and several columns: pool/lane labels, index ids, index sequences,
             the corresponding colors according to the chosen Illumina chemistry and a score equal to the minimum
             number of mismatches with the other indices of the pool/lane.")
    } else{
      "Please load indices and then press the \"Search for a solution\" button."
    }
  })
  
  # display the solution
  displaySolution <- eventReactive(input$go, {
    shinyjs::hide("proposedSolution")
    shinyjs::hide("visualization")
    shinyjs::hide("colorBalancing")
    # print("Start of displaySolution")
    if (is.null(input$multiplexingRate) | (is.null(fileInput()) & is.null(fileInput2()) & testData()=="none")){
      return(NULL)
    } else{
      withProgress({
        index <- inputIndex()[, -5]
        index2 <- inputIndex2()[, -5]
        if (is.null(index2)){
          indexesList <- generateList()
          indexesList2 <- NULL
        } else{
          if (input$i7i5pairing){
            # come back to a kind of single-indexing
            if (input$chemistry == "2"){
              mask <- which(substr(index$sequence, 1, 2) == "GG" | substr(index2$sequence, 1, 2) == "GG")
              if (length(mask) > 0){
                index <- index[-mask,]
                index2 <- index2[-mask,]
              }
            }
            names(index2) <- paste0(names(index2), "2")
            index <- cbind(index, index2)
            indexesList <- generateListPairedIndexes()
            index2 <- NULL
            indexesList2 <- NULL
          } else{
            indexesList <- generateList()
            indexesList2 <- generateList2()
          }
        }
        # print(paste0("Index: ", index))
        # print(paste0("Index2: ", index2))
        solution <- findSolution(indexesList = indexesList,
                                 index = index,
                                 indexesList2 = indexesList2,
                                 index2 = index2,
                                 nbSamples = as.numeric(input$nbSamples),
                                 multiplexingRate = as.numeric(input$multiplexingRate),
                                 unicityConstraint = ifelse(is.null(inputIndex2()), input$unicityConstraint, "none"),
                                 nbMaxTrials = as.numeric(input$nbMaxTrials),
                                 completeLane = input$completeLane,
                                 selectCompIndexes = input$selectCompIndexes,
                                 chemistry = input$chemistry,
                                 i7i5pairing = input$i7i5pairing)
        
        # Calculate the percentage of each character in each position of color1 and color2
        solution_color_percentages_with_weights <- calculate_color_percentages_with_weights(solution)
        colorPercentagesFormattedOutputWide <- convert_to_formatted_output_wide(solution_color_percentages_with_weights)
        # areIndicesCompatible <- !sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry)
        areIndicesCompatible = all(solution$areIndicesCompatible)
        # print(paste0("split(solution, solution$pool): ", split(solution, solution$pool)))
        # print(paste0("solution: ", solution))
        # print(paste0("areIndicesCompatible: ", areIndicesCompatible))
        
      }, message="R is looking for a solution...", max=0)
      shinyjs::show("proposedSolution")
      shinyjs::show("visualization")
      shinyjs::show("colorBalancing")
      # return(solution)
      # print(paste0("Solution: ", solution))
      return(list(solution = solution, 
                  colorPercentagesFormattedOutputWide = colorPercentagesFormattedOutputWide,
                  areIndicesCompatible = areIndicesCompatible))
    }
  })
  
  output$solution <- DT::renderDT({
    tryCatch({displaySolution()$solution}, error = function(e) NULL)
  }, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  
  # text describing the color balancing
   output$textDescribingColorBalancing <- renderText({
    if (!is.null(tryCatch({displaySolution()$colorPercentagesFormattedOutputWide}, error = function(e) NULL))){
      paste("Below are the color percentages for", as.numeric(input$nbSamples)/as.numeric(input$multiplexingRate),
            "pool(s) of", input$multiplexingRate, "samples using the parameters specified. The table contains 
             one row for each position (1 - 8) in each color column (color1 and color2), and shows the respective 
             percentages for each of the red, green, and blue colors in that position (the '-' shows the percentage 
             where no color is present).")
    } else{
      "Please load indices and then press the \"Search for a solution\" button."
    }
  })
  
  # Render color balancing data table
  # output$colorPercentagesFormattedOutputWide <- DT::renderDT({
  #  tryCatch({displaySolution()$colorPercentagesFormattedOutputWide}, error = function(e) NULL)
  # }, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  output$colorPercentagesFormattedOutputWide <- DT::renderDT({
     tryCatch({displaySolution()$colorPercentagesFormattedOutputWide}, error = function(e) NULL)
   }, options = list(
     paging = FALSE, 
     searching = FALSE, 
     info = FALSE,
     columnDefs = list(
       list(className = 'dt-right', targets = c(4,5,6,7))  # Right align columns 4-7 (R,G,B,-)
     )
  ))
   
  # text describing the heatmap
  output$textDescribingHeatmap <- renderText({
    if (!is.null(tryCatch({displaySolution()$solution}, error = function(e) NULL))){
      paste0("The plot below allows to vizualize the proposed solution. Samples (in rows) are grouped by pool/lane
             and each nucleotide of each index is displayed with a color according to the chosen Illumina chemistry. 
             One can thus quickly check whether each color is used at each position. Note that sample ids (from 1 to ",
             as.numeric(input$nbSamples), ") are printed on the left while index ids and sample weights are printed ",
             "on the right, followed by sample weights in parentheses.")
    } else{
      "Please load indices and then press the \"Search for a solution\" button."
    }
  })
  
  output$compatibilityWarning <- renderText({
    if (!is.null(tryCatch({displaySolution()$areIndicesCompatible}, error = function(e) NULL)) && !displaySolution()$areIndicesCompatible){
      return("****  Warning: Colors in indices are not compatible in at least one pool/lane  ****")
    }
    return("")
  })
  
  output$compatibilityWarning2 <- renderText({
    if (!is.null(tryCatch({displaySolution()$areIndicesCompatible}, error = function(e) NULL)) && !displaySolution()$areIndicesCompatible){
      return("****  Warning: Colors in indices are not compatible in at least one pool/lane  ****")
    }
    return("")
  })
  
  output$compatibilityWarning3 <- renderText({
    if (!is.null(tryCatch({displaySolution()$areIndicesCompatible}, error = function(e) NULL)) && !displaySolution()$areIndicesCompatible){
      return("****  Warning: Colors in indices are not compatible in at least one pool/lane  ****")
    }
    return("")
  })
  
  # plot of the solution
  output$heatmapindex <- renderPlot({heatmapindex(displaySolution()$solution)}, res=90)
  output$heatmapindex2 <- renderUI({
    if (!is.null(tryCatch({displaySolution()$solution}, error = function(e) NULL))){
      plotOutput("heatmapindex", width=900, height=220+20*as.numeric(input$nbSamples))
    }
  })
  
  # download the solution
  output$downloadData <- downloadHandler(
    filename = "chosenIndexes.txt",
    content = function(file) write.table(displaySolution()$solution, file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  )
  output$downloadButton <- renderUI({
    if (!is.null(tryCatch({displaySolution()$solution}, error = function(e) NULL))) downloadButton("downloadData", "Download")
  })
  
  # download the color balancing solution
  output$downloadColorBalanceData <- downloadHandler(
    filename = "colorbalancepercentages.txt",
    content = function(file) {
      df <- displaySolution()$colorPercentagesFormattedOutputWide
      # Replace "NA" with spaces
      df[is.na(df) | df == "NA"] <- "  "
      write.table(df, file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    }
    # content = function(file) write.table(displaySolution()$colorPercentagesFormattedOutputWide, file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  )
  output$downloadColorBalanceButton <- renderUI({
    if (!is.null(tryCatch({displaySolution()$colorPercentagesFormattedOutputWide}, error = function(e) NULL))) downloadButton("downloadColorBalanceData", "Download")
  })
  
})
