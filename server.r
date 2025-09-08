library(shiny)

options(shiny.sanitize.errors = FALSE,   # to display informative error messages
        shiny.maxRequestSize = 5*1024^2) # limit size for the input file to upload (5Mo here)

shinyServer(function(input, output, session) {
  
  # Auto kill app when browser window closed (https://github.com/daattali/advanced-shiny/tree/master/auto-kill-app)
  session$onSessionEnded(stopApp)
  
  # reset input parameters when pressing the reset button
  observeEvent(input$reset, {
    shinyjs::reset("allParameters")
    # previous inputFile disapears but is still sent to the server
    # known issue: https://github.com/daattali/shinyjs/issues/104
    rv$inputFile <- NULL
    rv$inputFile2 <- NULL
    rv$testdata <- "none"
    
    initialSelected(data.frame(id = character(0), stringsAsFactors = FALSE))   # drop all pre-/user-selected rows
    DT::selectRows(DT::dataTableProxy("inputIndex"), NULL)                     # clear the blue
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
  
  ## a mutable holder for the pre-selected rows
  initialSelected <- reactiveVal(
    data.frame(id = character(0), stringsAsFactors = FALSE)
  )
  
  # list of input index 1
  inputIndexes <- reactive({
    # print('Input index 1 entered')
    if (testData() %in% c("simple", "dual")){
      file <- ifelse(testData() == "simple", "www/inputIndexesExample.txt", "www/testCheckMyIndex-i7-i5.txt")
      # file <- ifelse(testData() == "simple", "www/inputIndexesExample.txt", "index96_UDI-i7.txt")
      # file <- ifelse(testData() == "simple", "www/inputIndexesExample.txt", "www/index24-i7.txt")
    } else{
      if (!is.null(fileInput())) file <- fileInput()$datapath else return(NULL)
    }
    # print(paste0('file: ', file))
    
    result <- tryCatch(
      readIndexesFileWithWeights(file),
      
      error = function(e) {
        
        msg <- conditionMessage(e)      # the message from readIndexesFileWithWeights()
        
        if (grepl("cannot open the connection", msg, fixed = TRUE)) {
          ## file path wrong / permission denied
          stop("Index file could not be opened (check the path and permissions):\n", msg)
          
        } else {
          ## fallback: rethrow with extra context
          stop("Unexpected error while loading index file:\n", msg)
        }
      }
    )
    
    ## Get rows where selected == 1 ----------------------------
    if ("selected" %in% names(result$index)) {
      sel <- result$index[result$index$selected == 1, "id", drop = FALSE]   # 1-column df
      if (nrow(sel)) initialSelected(sel)
      # initialSelected(result$index[result$index$selected == 1, , drop = FALSE])
    }
    
    result$index  <- result$index [ , !(names(result$index)  == "selected"), drop = FALSE]
    if (!is.null(result$index2) && "selected" %in% names(result$index2)) {
      result$index2 <- result$index2[ , !(names(result$index2) == "selected"), drop = FALSE]
    }
    
    # result <- tryCatch({readIndexesFileWithWeights(file)}, 
    #                   error = function(e) stop("An error occured when loading index 1 file, please check its structure."))
    
    # index <- tryCatch({readIndexesFileWithWeights(file)$index}, 
    #                   error = function(e) stop("An error occured when loading index 1 file, please check its structure."))
    # index <- addColors(index, input$chemistry)
    # print(paste0('index$color: ', index$color))
    # index$score <- scores(index$sequence)
    # return(index)
    
    if (!is.null(result$index2)) {
      # Case with four columns
      index <- addColors(result$index, input$chemistry)
      index$score <- scores(index$sequence)
      index2 <- addColors(result$index2, input$chemistry)
      index2$score <- scores(index2$sequence)
    } else {
      # Case with two columns
      index <- addColors(result$index, input$chemistry)
      index$score <- scores(index$sequence)
      index2 <- NULL

      # Check inputFile2 as normal
      if (!is.null(fileInput2())) {
        file2 <- fileInput2()$datapath

        if (file.exists(file2)) {
          result <- tryCatch(
            readIndexesFileWithWeights(file2),
            error = function(e) {
              msg <- conditionMessage(e)      # the message from readIndexesFileWithWeights()
              if (grepl("cannot open the connection", msg, fixed = TRUE)) {
                ## file path wrong / permission denied
                stop("Index file could not be opened (check the path and permissions):\n", msg)
              } else {
                ## fallback: rethrow with extra context
                stop("Unexpected error while loading index file:\n", msg)
              }
            }
          )
          index2 <- result$index
          index2 <- addColors(index2, input$chemistry)
          index2$score <- scores(index2$sequence)
        }
      }
    }
    
    return(list(index = index, index2 = index2))
  })
  
  inputIndex <- reactive({
    return(inputIndexes()$index)
  })
  
  inputIndex2 <- reactive({
    return(inputIndexes()$index2)
  })
  
  output$indexUploaded <- reactive({!is.null(inputIndex())})
  outputOptions(output, "indexUploaded", suspendWhenHidden=FALSE)
  
  # --------------------------------------------------------------------
  # 1) Show the table with multi-row selection enabled
  # 2) Make the selected rows easy to use elsewhere
  # --------------------------------------------------------------------
  output$inputIndex <- DT::renderDT({
    
    ## ── 1.  get the data you want to show --------------------------------
    df <- req(inputIndexes()$index)          # already has NO “selected” column
    
    ## ── 2.  find which rows were pre-selected in the file -----------------
    pre_rows <- which(df$id %in% isolate(initialSelected()$id))  # 1-based
    
    ## ── 3.  build the widget with those rows pre-highlighted -------------
    DT::datatable(
      df,
      rownames  = FALSE,
      options   = list(
        paging    = FALSE,
        searching = FALSE,
        info      = FALSE
      ),
      selection = list(
        mode     = "multiple",
        target   = "row",
        selected = pre_rows                # <- <-- HERE: blue rows on load
      )
    )
  })
  
  # Old code
  # output$inputIndex <- DT::renderDT({inputIndex()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  
  # Any downstream code can use selectedRows() to see what was chosen
  # selectedRows <- reactive({
  #   req(input$inputIndex_rows_selected)             # wait until at least one row is clicked
  #    inputIndex()[ input$inputIndex_rows_selected , ]
  # })
  
  # selectedRows <- reactive({
  #   rows <- input$inputIndex_rows_selected        # integer vector or integer(0)
  #   
  #   if (length(rows) == 0) {
  #     # keep the same columns but zero rows
  #     inputIndex()[0, , drop = FALSE]
  #   } else {
  #     inputIndex()[rows, , drop = FALSE]
  #   }
  # })

  observeEvent(input$inputIndex_rows_selected, 
               ignoreNULL = FALSE,          # ← fire when selection becomes empty
  {
    
    df  <- req(inputIndexes()$index)                 # full table
    ids <- df$id[input$inputIndex_rows_selected]     # may be character(0)
    
    # overwrite the reactiveVal with the *current* set of IDs
    initialSelected(data.frame(id = ids, stringsAsFactors = FALSE))
  })
  
  selectedRows <- reactive({
    initialSelected()          # already a 1-column data-frame (id)
  })
  
  # Example: print the chosen rows to the console
  observeEvent(selectedRows(), {
    cat("--- rows the user clicked ---\n")
    print(selectedRows())
  })
  
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
  # inputIndex2 <- reactive({
  #   # print('Input index 2 entered')
  #   inputIndex()
  #   if (testData() == "simple") return(NULL)
  #   if (testData() == "dual"){
  #      # file2 <- "www/index24-i5.txt"
  #      file2 <- "index96_UDI-i5.txt"  # //---
  #   } else{
  #      if (!is.null(fileInput2())){
  #        if (is.null(inputIndex())) stop("Please load index 1 (i7) first.")
  #        file2 <- fileInput2()$datapath
  #      } else {
  #        if (!is.null(inputIndex())){
  #          file2 <- fileInput()$datapath
  #        } else {
  #          return(NULL)
  #        }
  #      }
  #   }
  #   print(paste0('file2: ', file2))
    
  #   # index2 <- tryCatch({readIndexesFileWithWeights(file2)$index}, 
  #   #                    error = function(e) stop("An error occured when loading index 2 file, please check its structure."))
  #   if (testData() == "dual" | !is.null(fileInput2())){
  #      index2 <- tryCatch({readIndexesFileWithWeights(file2)$index}, 
  #                     error = function(e) stop("An error occured when loading index 2 file, please check its structure."))
  #   } else {
  #      result <- tryCatch({readIndexesFileWithWeights(file2)}, 
  #                       error = function(e) stop("An error occured when loading index 2 file, please check its structure."))
  #      if (!is.null(result$index2)) {
  #          index2 <- result$index2
  #      } else {
  #          return(NULL)
  #      }
  #   }
  #   index2 <- addColors(index2, input$chemistry)
  #   # print(paste0('index2$color: ', index2$color))
  #   index2$score <- scores(index2$sequence)
  #   return(index2)
  # })
  
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
  
  # generate list(s) of indices (keep existing weights)
  generateListKeepWeights <- reactive({
    index <- inputIndex()
    # print(paste0('column names: ', names(index)))
    # print(paste0('index[, -5]: ', index[, -5]))
    # print(paste0('input$multiplexingRate: ', input$multiplexingRate))
    
    return(generateListOfIndexesCombinations(index = index[, -6],
                                             nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                             completeLane = input$completeLane,
                                             selectCompIndexes = input$selectCompIndexes,
                                             chemistry = input$chemistry,
                                             selectedRows = selectedRows()))
  })

  # generate list(s) of indices (set weights all weights to 1)
  generateListRemoveWeights <- reactive({
    index <- inputIndex()
    index$weight <- 1
    # print(paste0('column names: ', names(index)))
    # print(paste0('index[, -5]: ', index[, -5]))
    # print(paste0('input$multiplexingRate: ', input$multiplexingRate))
    
    return(generateListOfIndexesCombinations(index = index[, -6],
                                             nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                             completeLane = input$completeLane,
                                             selectCompIndexes = input$selectCompIndexes,
                                             chemistry = input$chemistry,
                                             selectedRows = selectedRows()))
  })

  generateList2 <- reactive({
    index2 <- inputIndex2()
    if (is.null(index2)){
      return(NULL)
    } else{
      index2$weight <- 1
      return(generateListOfIndexesCombinations(index = index2[, -6],
                                               nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                               completeLane = input$completeLane,
                                               selectCompIndexes = input$selectCompIndexes,
                                               chemistry = input$chemistry,
                                               selectedRows = selectedRows()))
    }
  })
  generateListPairedIndexes <- reactive({
    index <- inputIndex()
    index2 <- inputIndex2()
    if (is.null(index) & is.null(index2)){
      return(NULL)
    } else{
      # if (!input$i7i5pairing) {
      #   index$weight <- 1
      #   index2$weight <- 1
      # }
      if (input$chemistry == "2"){
        mask <- which(substr(index$sequence, 1, 2) == "GG" | substr(index2$sequence, 1, 2) == "GG")
        if (length(mask) > 0){
          index <- index[-mask,]
          index2 <- index2[-mask,]
        }
      }
      names(index2) <- paste0(names(index2), "2")
      index <- cbind(index[, -6], index2[, -6])
      return(generateListOfIndexesCombinations(index = index,
                                               nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                               completeLane = input$completeLane,
                                               selectCompIndexes = input$selectCompIndexes,
                                               chemistry = input$chemistry,
                                               selectedRows = selectedRows()))
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
  
  # Test function to print G / C counts
  test_GC_counts <- function(colors) {
    
    # Convert the strings to a character matrix (rows = strings, columns = positions)
    char_matrix <- do.call(rbind, strsplit(colors, split = ""))
    
    # Count the number of G or C characters at each position
    gc_counts <- colSums(char_matrix == "G" | char_matrix == "C")
    
    # Output the result
    print(paste0('gc_counts: ', gc_counts))
    
  }
  
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
        index <- inputIndex()[, -6]
        index2 <- inputIndex2()[, -6]
        if (is.null(index2)){
          indexesList <- generateListKeepWeights()
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
            index$weight <- 1
            index2$weight <- 1
            indexesList <- generateListRemoveWeights()
            indexesList2 <- generateList2()
          }
        }
        
        # print(paste0("indexesList: ", indexesList))  # //--- comment me out
        # testWeight <- index$weight  # //--- comment me out
        # print("Index$weight: ")  # //--- comment me out
        # print(testWeight)  # //--- comment me out
        
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
        results <- calculate_color_percentages_with_weights(solution)
        solution_color_percentages_with_weights <- results$solution_color_percentages
        colorPercentagesFormattedOutputWide <- convert_to_formatted_output_wide(solution_color_percentages_with_weights)
        # areIndicesCompatible <- !sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry)
        # areIndicesCompatible = all(solution$areIndicesCompatible)
        
        # print('displaySolution')
        # split_result <- split(solution, solution$pool)
        # print(paste0("split_result[[1]]: ", split_result[[1]]))
        # print(paste0("names(split_result[[1]]): ", names(split_result[[1]])))
        
        # result <- areIndexesCompatible(split_result[[1]], input$chemistry, column="color1")
        # print(paste0("result: ", result))
        # print(paste0("split(solution, solution$pool): ", split(solution, solution$pool)))
        # print(paste0("solution: ", solution))
        # print(paste0("solution$areIndicesCompatible: ", solution$areIndicesCompatible))
        # print(paste0("sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, color1): ",sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, "color1")))

        # print(paste0("all(sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, color1)): ",all(sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, "color1"))))
        # print(paste0("all(sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, color2)): ",all(sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, "color2"))))
        
        # test_GC_counts(solution$color1)
        # print(paste0("solution_color_percentages_with_weights: ", solution_color_percentages_with_weights))
        
        # if (all(sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, "color1")) & all(sapply(split(solution, solution$pool), areIndexesCompatible, input$chemistry, "color2"))) {
        #   areIndicesCompatible <- TRUE
        # } else {
        #   areIndicesCompatible <- FALSE
        # }
        
        # Check which columns exist
        if ("color" %in% names(solution)) {
          # Single color case
          areIndicesCompatible <- all(sapply(split(solution, solution$pool), 
                                             areIndexesCompatible, 
                                             input$chemistry, 
                                             "color"))
        } else if (all(c("color1", "color2") %in% names(solution))) {
          # Dual color case
          areIndicesCompatible <- all(sapply(split(solution, solution$pool), 
                                             areIndexesCompatible, 
                                             input$chemistry, 
                                             "color1")) & 
                                  all(sapply(split(solution, solution$pool), 
                                             areIndexesCompatible, 
                                             input$chemistry, 
                                            "color2"))
        }
        
        hammingDistanceScoreCutoff <- 3
        hammingDistanceScoreGreaterThanCuttoff <- TRUE
        
        # TRUE if every index in this subset has min Hamming distance >= cutoff
        hamming_ok <- function(df, cutoff, column = "score") {
          seqs <- df[[column]]
          seqs <- seqs[!is.na(seqs)]
          if (length(seqs) <= 1L) return(TRUE)
          # "scores" returns, for each seq, min distance to the others
          min(scores(seqs)) >= cutoff   # equivalent to all(scores(seqs) >= cutoff)
        }
        
        pool_list <- split(solution, solution$pool)
        
        if ("sequence" %in% names(solution)) {
          # Single-index case
          hammingDistanceScoreGreaterThanCuttoff <- all(vapply(
            pool_list,
            function(df) hamming_ok(df, hammingDistanceScoreCutoff, "sequence"),
            logical(1)
          ))
          
        } else if (all(c("sequence1","sequence2") %in% names(solution))) {
          # Dual-index case: require BOTH i7 and i5 to meet the cutoff
          hammingDistanceScoreGreaterThanCuttoff <- all(vapply(
            pool_list,
            function(df) {
              ok1 <- hamming_ok(df, hammingDistanceScoreCutoff, "sequence1")
              ok2 <- hamming_ok(df, hammingDistanceScoreCutoff, "sequence2")
              ok1 && ok2
            },
            logical(1)
          ))
        }
        
        # print(paste0("areIndicesCompatible: ", areIndicesCompatible))
        # test_GC_counts(c(solution$color1, solution$color2))
        
        # print(paste0("solution_color_percentages_with_weights: ", solution_color_percentages_with_weights))
        
      }, message="App is looking for a solution...", max=0)
      shinyjs::show("proposedSolution")
      shinyjs::show("visualization")
      shinyjs::show("colorBalancing")
      # return(solution)
      # print(paste0("Solution: ", solution))
      return(list(solution = solution, 
                  colorPercentagesFormattedOutputWide = colorPercentagesFormattedOutputWide,
                  areIndicesCompatible = areIndicesCompatible,
                  hammingDistanceScoreGreaterThanCuttoff = hammingDistanceScoreGreaterThanCuttoff))
    }
  })
  
  # Modify original renderDT to highlight selected rows (all rows whose id is in selectedRows() now show bold, deep-pink text)
  # output$solution <- DT::renderDT({
  #   tryCatch({displaySolution()$solution}, error = function(e) NULL)
  # }, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  output$solution <- DT::renderDT({
    
    df <- tryCatch(displaySolution()$solution, error = function(e) NULL)
    req(df)

    ## Which column holds the barcode ID?  ("id" OR "id1")
    id_col <- if ("id"  %in% names(df))      "id"
    else if ("id1" %in% names(df)) "id1"
    else                           NULL           # no ID column found

    ## Vector of IDs that must be highlighted ---------------------------
    hilite <- if (!is.null(id_col) && !is.null(selectedRows()))
      hilite <- selectedRows()$id
    else character(0)
    # print(paste0('hilite: ', hilite))
    
    ## **Index of the column to hide (0‑based for JS)**
    hide_idx <- which(names(df) == "areIndicesCompatible") - 1
    
    ## Build the table ------------------------------------------------
    DT::datatable(
      df,
      rownames = FALSE,
      options  = list(
        paging    = FALSE,
        searching = FALSE,
        info      = FALSE,
        columnDefs  = if (length(hide_idx))
          list(list(visible = FALSE, targets = hide_idx))
        else NULL        # do nothing if column not present
      )
    ) %>%
      ## Colour EVERY column in the matching rows --------------------
    {                                   # apply styling only if id_col exists
      if (!is.null(id_col) && length(hilite)) {
        DT::formatStyle(
          .,
          columns        = names(df),                 # colour the whole row
          valueColumns   = id_col,                    # "id" or "id1"
          # backgroundColor = DT::styleEqual(hilite, "#f0f0f0"), # grey fill
          color           = DT::styleEqual(hilite, "deeppink"), # dark-pink text
          fontWeight      = DT::styleEqual(hilite, "bold"),    # bold faces
          target          = "row"
        )
      } else .
    }
  })
  
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
    if (!is.null(tryCatch({displaySolution()$hammingDistanceScoreGreaterThanCuttoff}, error = function(e) NULL)) && !displaySolution()$hammingDistanceScoreGreaterThanCuttoff){
      return("****  Warning: Score (distance between indices) is less that required minimum  ****")
    }
    else if (!is.null(tryCatch({displaySolution()$areIndicesCompatible}, error = function(e) NULL)) && !displaySolution()$areIndicesCompatible){
      return("****  Warning: Colors in indices are not compatible in at least one pool/lane  ****")
    } else {
      return("")
    }
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
  output$heatmapindex <- renderPlot({
    req(selectedRows())                    # wait until at least one row selected
    req(displaySolution())                 # make sure the solution is ready
    
    heatmapindex(
      displaySolution()$solution,
      selectedRows = selectedRows()        # ← evaluated *here*, inside reactivity
    )
  }, res = 90)
  
  # output$heatmapindex <- renderPlot({heatmapindex(displaySolution()$solution)}, res=90)
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
