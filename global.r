# this file contains the functions used in both the "checkMyIndex" script and shiny application
library(parallel)
library(dplyr)
library(tidyr)
library(stringr)
# library(plotrix)            # boxedtext()


checkMyIndexVersion <- "1.0.2"
cornellCheckMyIndexVersion <- "1.4.3"

readIndexesFile <- function(file){
  index <- tryCatch({read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE, col.names=c("id","sequence"))},
                    error = function(e) stop("An error occured when loading the input file, please check its structure."))
  checkInputIndexes(index)
  return(index)
}

# The readIndexesFileWithWeights function reads the first row to detect the number of columns.
# If there are two columns, it adds a "weight" column with a default value of 1.
# If there are four columns, it creates two data frames: index with columns 1, 2, and 4, and index2 with columns 1, 3, and 4.
# Returns both data frames if applicable and checks both using checkInputIndexes.
readIndexesFileWithWeights <- function(file) {

  ## ───────────────────────────────────────────────────────────────────
  ## helper — stop if one or more columns are empty
  ## ───────────────────────────────────────────────────────────────────
  check_empty_columns <- function(df, context = "input file") {
    empty_cols <- which(vapply(df, function(col)
      length(col) == 0 || all(is.na(col)) || all(trimws(col) == ""),
      logical(1)
    ))
    if (length(empty_cols)) {
      stop(
        sprintf(
          "%s error: column(s) %s are empty. Each column must contain at least one non-empty value.",
          context,
          paste(empty_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  
  out <- tryCatch({
    # Read the first row to determine the number of columns
    sample_data <- read.table(file, header=FALSE, sep="\t", nrows=1, stringsAsFactors=FALSE)
    num_cols <- ncol(sample_data)
    # print(paste0('sample_data: ', sample_data))
    
    ## helper: TRUE only when every entry is a non-empty string of A/C/G/T
    is_valid_barcode <- function(x) {
      all(grepl("^[ACGT]+$", x))    # at least one char + only those four letters
    }
    
    if (num_cols == 2) {
      
      index <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                          col.names = c("id", "sequence"))
      check_empty_columns(index, "Two-column index file")
      
      if (!is_valid_barcode(index$sequence)) {
        stop("Two-column index file error: column 2 must contain non-empty barcode ",
             "sequences composed only of the letters A, C, G and T.")
      }
      
      index$weight <- 1
      index$selected <- 0
      colnames(index) <- c("id", "sequence", "weight", "selected")
      list(index = index, index2 = NULL)
      
    } else if (num_cols == 3) {
      
      full_data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      check_empty_columns(full_data, "Three-column index file")
      
      ## columns 2 and 3 are the barcode sequences
      if (is_valid_barcode(full_data[[2]]) && is_valid_barcode(full_data[[3]])) {
        
        index  <- full_data[, c(1, 2)]
        
        index2 <- full_data[, c(1, 3)]
        
        index$weight <- 1
        index2$weight <- 1
        index$selected <- 0
        index2$selected <- 0
        colnames(index) <- c("id", "sequence", "weight", "selected")
        colnames(index2) <- c("id", "sequence", "weight", "selected")
        
        list(index = index, index2 = index2)
        
      } else if (is_valid_barcode(full_data[[2]])) {
        index  <- full_data[, c(1, 2, 3)]
        index$selected <- 0
        colnames(index) <- c("id", "sequence", "weight", "selected")
        list(index = index, index2 = NULL)
        
      } else {
        stop("Three-column index file error: columns 2 *and* 3 must each contain ",
             "non-empty barcode sequences composed only of A, C, G and T, or column 2 ",
             "must contain a non-empty barcode sequence and column 3 must contain a weight.")
      }
      
    } else if (num_cols == 4) {
      
      full_data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      check_empty_columns(full_data, "Four-column index file")
      
      ## columns 2 and 3 are the barcode sequences
      if (is_valid_barcode(full_data[[2]]) && is_valid_barcode(full_data[[3]])) {
        
        index  <- full_data[, c(1, 2, 4)]
        index2 <- full_data[, c(1, 3, 4)]
        
        index$selected <- 0
        index2$selected <- 0
        colnames(index) <- c("id", "sequence", "weight", "selected")
        colnames(index2) <- c("id", "sequence", "weight", "selected")
        
        list(index = index, index2 = index2)
        
      } else if (is_valid_barcode(full_data[[2]])) {
        index  <- full_data[, c(1, 2, 3, 4)]
        colnames(index) <- c("id", "sequence", "weight", "selected")
        
        list(index = index, index2 = NULL)
        
      } else {
        stop("Four-column index file error: columns 2 *and* 3 must each contain ",
             "non-empty barcode sequences composed only of A, C, G and T, or column ",
             "must contain a non-empty barcode sequence and columns 3 and 4 must ",
             "contain a weight and a selected indicator.")
      }
      
    } else if (num_cols == 5) {
    
      full_data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      check_empty_columns(full_data, "Five-column index file")
      
      ## columns 2 and 3 are the barcode sequences
      if (!is_valid_barcode(full_data[[2]]) || !is_valid_barcode(full_data[[3]])) {
        stop("Five-column index file error: columns 2 *and* 3 must each contain ",
             "non-empty barcode sequences composed only of A, C, G and T.")
      }
      
      index  <- full_data[, c(1, 2, 4, 5)]
      colnames(index) <- c("id", "sequence", "weight", "selected")
      
      index2 <- full_data[, c(1, 3, 4, 5)]
      colnames(index2) <- c("id", "sequence", "weight", "selected")
      
      list(index = index, index2 = index2)
    
    } else {
      
      stop("The input file must contain either two, three, four, or five tab-separated columns.")
    }
  },
  error = function(e) {
    stop(conditionMessage(e))
  })
  
  # Assuming checkInputIndexes is applicable to both indices
  if (!is.null(out$index)) {
    checkInputIndexes(out$index)
  }
  if (!is.null(out$index2)) {
    checkInputIndexes(out$index2)
  }
  
  return(out)
}

addColors <- function(index, chemistry){
  # convert bases into colors
  if (chemistry == "X"){
    # G has no color, A is blue, C is Blue+Green (Cyan) and T is green
    index$color <- gsub("T", "G", gsub("C", "C", gsub("A", "B", gsub("G", "-", index$sequence))))
  }
  if (chemistry == "4"){
    # A/C are red and G/T are green
    index$color <- gsub("G|T", "G", gsub("A|C", "R", index$sequence))
  }
  if (chemistry == "2"){
    # G has no color, A is orange, C is red and T is green
    index$color <- gsub("T", "G", gsub("C", "R", gsub("A", "O", gsub("G", "-", index$sequence))))
  } 
  if (chemistry == "1"){
    # G has no color, A is green, C is red and T is orange
    index$color <- gsub("A", "G", gsub("C", "R", gsub("T", "O", gsub("G", "-", index$sequence))))
  } 
  # print('addColors')
  # print(paste0('chemistry: ', chemistry))
  # print(paste0('index$color: ', index$color))
  return(index)
}

checkInputIndexes <- function(index){
  if (any(duplicated(index$id))){
    stop("Input index ids are not unique, please check your input file.")
  }
  if (any(duplicated(index$sequence))){
    stop("Input indices are not unique, please check your input file.")
  }
  if (!(all(unlist(strsplit(index$sequence,"")) %in% c("A","C","T","G")))){
    stop("Input indices contain other characters than A, T, C and G, please check your input file.")
  }
  if (length(table(sapply(index$sequence, nchar))) > 1){
    stop("Input indices do not have all the same length, please check your input file.")
  }
}

# function to compute the Hamming distance between two indices
dist2indexes <- function(seq1, seq2){
  seq1 <- unlist(strsplit(seq1, ""))
  seq2 <- unlist(strsplit(seq2, ""))
  if (length(seq1)!=length(seq2)) stop("Impossible to calculate the distance as the two sequences do not have the same length.")
  return(sum(seq1!=seq2))
}

# function to compute the minimum Hamming distance between N indices/sequences
distNindexes <- function(sequences){
  dists <- mapply(dist2indexes, 
                  rep(sequences, length(sequences)), 
                  rep(sequences, each=length(sequences)))
  dists <- matrix(dists, nrow=length(sequences))
  diag(dists) <- Inf
  return(dists)
}

# function to compute the scores
# the score of each index is the minimum number of mismatches with the others
scores <- function(sequences) apply(distNindexes(sequences), 2, min)

# Test function to print G / C counts
test_GC_counts <- function(colors) {

  # Convert the strings to a character matrix (rows = strings, columns = positions)
  char_matrix <- do.call(rbind, strsplit(colors, split = ""))
  
  # Count the number of G or C characters at each position
  gc_counts <- colSums(char_matrix == "G" | char_matrix == "C")
  
  # Output the result
  print(paste0('gc_counts: ', gc_counts))
  
}

areIndexesCompatible <- function(index, chemistry, column="color"){
  # return TRUE if the input indices are compatible (i.e. can be used within the same pool/lane)
  # print('areIndexesCompatible inside')  # //--- comment me out
  # print(paste0('column: ', column))  # //--- comment me out
  # print(paste0('index: ', index))  # //--- comment me out
  # print(paste0('nrow(index): ', nrow(index)))  # //--- comment me out
  if (nrow(index)==1) return(TRUE)
  # print(paste0('column: ', column))
  # print(paste0('index$color: ', index$color))
  # print(paste0('strsplit(index[, column], "")', strsplit(index[, column], "")))
  matColors <- do.call("rbind", strsplit(index[, column], ""))
  # print(paste0('matcolors: ', matColors))  # //--- comment me out
  if (chemistry == "4"){
    sumRed <- apply(matColors, 2, function(x) sum(x=="R"))
    sumGreen <- nrow(matColors) - sumRed
    # print(paste0('all(sumRed >= 1 & sumGreen >= 1):', (all(sumRed >= 1 & sumGreen >= 1))))  # //--- comment me out
    return(all(sumRed >= 1 & sumGreen >= 1))
  }
  
  if (chemistry == "X"){
    # Count colors for each position
    sumBlue <- apply(matColors, 2, function(x) sum(x=="B" | x=="C"))  # Count both B and C as blue
    sumGreen <- apply(matColors, 2, function(x) sum(x=="G" | x=="C")) # Count both G and C as green
    sumNoColor <- apply(matColors, 2, function(x) sum(x=="-"))
    
    # For each position, check if:
    # - there's at least one green OR
    # - there's a mixture of blue and no color
    # AND ensure it's not all blue
    # positionCheck <- mapply(function(b, g, nc) {
    #  (g >= 1) || (b >= 1 && nc >= 1)
    # }, sumBlue, sumGreen, sumNoColor)
    
    # Return true only if all positions pass and at least one position has green
    # return(all(positionCheck) && any(sumGreen >= 1))
    
    # print(paste0('index: ', index))  # //--- comment me out
    # print(paste0('sumGreen: ', sumGreen))  # //--- comment me out
    # print(paste0('all(sumGreen >= 1): ', all(sumGreen >= 1)))  # //--- comment me out
    
    # print('areIndexesCompatible:')
    # test_GC_counts(index$color)
    # stop("Halting script.")
    
    # Check if each position has at least one green
    return(all(sumGreen >= 1))
  }

  if (chemistry %in% c("1","2")){
    sumNoColor <- apply(matColors, 2, function(x) sum(x=="-"))
    return(all(sumNoColor < nrow(index)))
  }
}

generateListOfIndexesCombinations <- function(index, nbSamplesPerLane, completeLane, selectCompIndexes, chemistry, selectedRows){
  # remove indices starting with GG if two-channel chemistry
  if (chemistry == "2") index <- index[!sapply(index$sequence, substr, 1, 2) == "GG",]

  ## Save the *original* lane size
  initialNbSamplesPerLane <- nbSamplesPerLane
  
  ## Cap lane size to what is actually available
  nbSamplesPerLane <- min(nrow(index), nbSamplesPerLane)

  # print('generateListOfIndexesCombinations')
  # print(paste0('nrow(index): ', nrow(index)))
  # print(paste0('nbSamplesPerLane: ', nbSamplesPerLane))
  # print(paste0('completeLane: ', completeLane))
  # print(paste0('choose(nrow(index), nbSamplesPerLane): ', choose(nrow(index), nbSamplesPerLane)))
  
  ## Optionally shrink lane size to keep the search space sane
  while (!completeLane &&
         choose(nrow(index), nbSamplesPerLane) > 2e5 &&
         nbSamplesPerLane > 2) {
    nbSamplesPerLane <- nbSamplesPerLane - 1
  }
  
  # optimize nbSamplesPerLane to generate a reduced number of combinations of indices
  # while (!completeLane & choose(nrow(index), nbSamplesPerLane) > 2e5 & nbSamplesPerLane > 2) nbSamplesPerLane <- nbSamplesPerLane - 1
  
  # print(paste0('nbSamplesPerLane: ', nbSamplesPerLane))
  # stop if too many combinations to test
  if (choose(nrow(index), nbSamplesPerLane) > 1e9) stop("Too many candidate combinations of indices to easily find a solution, please use different parameters.")
  
  ## -------------------------------------------------------------------
  ## Gather the REQUIRED IDs coming from the Shiny table
  ## (handles NULL, vector, or data-frame consistently)
  ## -------------------------------------------------------------------
   tryCatch({
    
    # print('generateListOfIndexesCombinations')
    # isolate({
    #    message("---- selectedRows (raw) ----")
    #    message("class:  ", paste(class(selectedRows), collapse = ", "))
    #    message("typeof: ", typeof(selectedRows))
    #    print(selectedRows)                    # prints the function prototype
    # })
    sr <- selectedRows
    # print(paste0('sr:', sr))
    # -------------------------------------------------------------------
    if (is.null(sr)                                             ||        # nothing
        (is.data.frame(sr) && nrow(sr) == 0)                    ||
        (length(sr) == 0)) {
      
      required_ids <- character(0)
      
    } else if (is.data.frame(sr)) {                                      # data-frame
      
      if (!"id" %in% names(sr))
        stop("selectedRows is a data-frame but has no 'id' column")
      required_ids <- as.character(sr$id)
      
    } else {                                                             # vector
      required_ids <- as.character(sr)
    }
    
  },  # <---------------------------- end of try block ------------------
  
  # Catch *any* ordinary error
  error = function(e) {
    message("❌  ERROR inside selectedRows processing: ", e$message)
    print(e)                       # full condition (call, classes, etc.)
    NULL
  },
  
  # Catch req()/validate() aborts as well
  shiny.silent.error = function(e) {
    message("⚠️  req()/validate() aborted: ", conditionMessage(e))
    NULL
  })
  
  nRequired <- length(required_ids)
  # print(paste0('required_ids: ', required_ids))
  # print(paste0('nRequired: ', nRequired))
  
  ## Helper that turns a vector of IDs into the matching data.frame rows
  df_from_ids  <- function(ids) index[index$id %in% ids, , drop = FALSE]
  
  ## -------------------------------------------------------------------
  ##  CASE LOGIC
  ## -------------------------------------------------------------------
  ## CASE 1  ─ fewer required IDs than (reduced) lane size - generate all combinations and return those that contain required ids
  if (nRequired < nbSamplesPerLane) {
    
    Clust <- parallel::makeCluster(max(c(parallel::detectCores(logical = FALSE) - 1, 1)))
    
    # All size-nbSamplesPerLane combos drawn from EVERY index
    all_combos <- combn(index$id, nbSamplesPerLane, simplify = FALSE)
    
    if (nRequired == 0) {
      ## No “must-have” indices → keep every combination
      possibleCombinations <- all_combos
    } else {
      ## At least one required index → keep only combos that contain them all
      possibleCombinations <- Filter(function(x) all(required_ids %in% x), all_combos)
    }
    
    # Turn each vector of IDs into the full rows
    indexesCombinations <-
      parallel::parLapply(Clust, possibleCombinations, df_from_ids)
    
    # Optional compatibility filter
    if (selectCompIndexes)
      indexesCombinations <-
      indexesCombinations[
        parallel::parSapply(Clust, indexesCombinations,
                            areIndexesCompatible, chemistry = chemistry)
      ]
    
    parallel::stopCluster(Clust)
    
    ## CASE 2  ─ many required IDs, and they are ≥ initial lane size
  } else if (nRequired >= initialNbSamplesPerLane) {
    
    # All combos of size *initial* lane size, but using only required IDs
    possibleCombinations <-
      combn(required_ids, initialNbSamplesPerLane, simplify = FALSE)
    
    indexesCombinations <- lapply(possibleCombinations, df_from_ids)
    
    if (selectCompIndexes)
      indexesCombinations <-
      indexesCombinations[
        vapply(indexesCombinations,
               areIndexesCompatible,
               logical(1),
               chemistry = chemistry)
      ]
    
    ## CASE 3  ─ many required IDs, but fewer than the original lane size
  } else {  # nRequired >= nbSamplesPerLane  &&  nRequired < initialNbSamplesPerLane
    
    # Return a single list element that contains *all* required IDs
    indexesCombinations <- list(df_from_ids(required_ids))
    
    if (selectCompIndexes &&
        !areIndexesCompatible(indexesCombinations[[1]], chemistry = chemistry))
      warning("The only possible set of required indices is incompatible.")
  }
  
  return(indexesCombinations)
}

searchOneSolution <- function(indexesList, index, indexesList2=NULL, index2=NULL,
                              nbLanes, multiplexingRate, unicityConstraint, chemistry,
                              i7i5pairing){
  # goal: look for a solution (i.e. a combination of combination of indices) such that
  #  - each index is used only once if required (unicityConstraint = index)
  #  - each combination of indices is used only once if required (unicityConstraint = lane)
  # two steps:
  #  1) fill the lanes with as many samples per lane as in indexesList (may be lower than multiplexingRate for optimization purposes)
  #  2) if this number is lower than the desired multiplexing rate, complete the solution returned at step 1 adding indices
  # this function can return NULL if no solution is found (need to re-run in that case)
  inputNbSamplesPerLane <- nrow(indexesList[[1]])
  compatibleCombinations <- vector(mode="list", length=nbLanes)
  areIndicesCompatible <- TRUE
  
  # print('Got to searchOneSolution')  # //--- comment me out
  
  # single-indexing or paired dual-indexing
  if (is.null(index2) | is.null(indexesList2)){
    k <- 1
    # print(paste0('nbLanes: ', nbLanes))
    while (k <= nbLanes){
      
      #  print(paste0('compatibleCombinations: ', compatibleCombinations))
      if (length(indexesList) == 0){
        # no available combination of indices anymore, try again!
        return(NULL)
      } else{
        i <- sample(1:length(indexesList), 1, FALSE)
        if (!i7i5pairing){
          # if (areIndexesCompatible(indexesList[[i]], chemistry, "color")){
          # areIndicesCompatible <- areIndexesCompatible(indexesList[[i]], chemistry, "color")
          # print(paste0('areIndicesCompatible before: ', areIndexesCompatible(indexesList[[i]], chemistry, "color")))
          # print('Came from searchOneSolution1')
          if (is.null(areIndicesCompatible) || areIndicesCompatible) {
            areIndicesCompatible <- areIndexesCompatible(indexesList[[i]], chemistry, "color")
          }
          # print(paste0('areIndicesCompatible: ', areIndicesCompatible))
          compatibleCombinations[[k]] <- indexesList[[i]]
          # remove either the combination used or all the combinations for which an index has already been selected
          if (unicityConstraint=="index") indexesList <- indexesList[!sapply(indexesList, function(tab) any(tab$id %in% indexesList[[i]]$id))]
          if (unicityConstraint=="lane") indexesList <- indexesList[-i]
          k <- k+1
          # } else{
          #   indexesList <- indexesList[-i]
          # }
        } else{
          # if (areIndexesCompatible(indexesList[[i]], chemistry, "color") & areIndexesCompatible(indexesList[[i]], chemistry, "color2")){
          # areIndicesCompatible <- (areIndexesCompatible(indexesList[[i]], chemistry, "color") & areIndexesCompatible(indexesList[[i]], chemistry, "color2"))
          # print(paste0('areIndicesCompatible before: ', (areIndexesCompatible(indexesList[[i]], chemistry, "color") & areIndexesCompatible(indexesList[[i]], chemistry, "color2"))))
          # print('Came from searchOneSolution2')
          # print(paste0('i: ', i))
          # print(paste0('indexesList[[i]]: ', indexesList[[i]]))
          if (is.null(areIndicesCompatible) || areIndicesCompatible) {
            areIndicesCompatible <- (areIndexesCompatible(indexesList[[i]], chemistry, "color") & areIndexesCompatible(indexesList[[i]], chemistry, "color2"))
          }
          # print(paste0('areIndicesCompatible: ', areIndicesCompatible))
          compatibleCombinations[[k]] <- indexesList[[i]]
          k <- k+1
          # } else{
          #   indexesList <- indexesList[-i]
          # }
        }
      }
      
    }
    
    solution <- data.frame(sample=1:(nbLanes*inputNbSamplesPerLane), 
                           pool=rep(1:nbLanes, each=inputNbSamplesPerLane), 
                           do.call("rbind", compatibleCombinations))
    if (multiplexingRate > inputNbSamplesPerLane){
      # only a partial solution has been found, need to complete it
      if (chemistry == "2"){
        index <- index[!(sapply(index$sequence, substr, 1, 2) == "GG"),]
        if (i7i5pairing) index <- index[!(sapply(index$sequence2, substr, 1, 2) == "GG"),]
      }
      solution <- completeSolution(partialSolution = solution,
                                   index = index,
                                   multiplexingRate = multiplexingRate,
                                   unicityConstraint = unicityConstraint)
    }
    if (!is.null(solution)){
      if (i7i5pairing){
        names(solution)[3:6] <- paste0(names(solution)[3:6], "1")
        solution$score1 <- unlist(tapply(solution$sequence1, solution$pool, scores))
        solution$score2 <- unlist(tapply(solution$sequence2, solution$pool, scores))
        # re-order columns
        # print(paste0('Solution names: ', names(solution)))
        # print(paste0('Solution before reorder columns: ', solution$weight1))
        solution <- solution[, c("sample", "pool", 
                                 "id1", "sequence1", "color1", "score1", "weight1", 
                                 "id2", "sequence2", "color2", "score2", "weight2")]
      } else{
        solution$score <- unlist(tapply(solution$sequence, solution$pool, scores))
      }
    }
    
    solution$areIndicesCompatible <- areIndicesCompatible
    
    # if (is.null(solution$areIndicesCompatible) || all(solution$areIndicesCompatible)) {
    #   solution$areIndicesCompatible <- ave(rep(TRUE, nrow(solution)), solution$pool, 
    #                                        FUN=function(x) areIndicesCompatible)
    # }
    
    # solution$areIndicesCompatible <- unlist(tapply(solution$sequence1, solution$pool, 
    #                                                function(x) areIndicesCompatible))
    
    # if (is.null(solution$areIndicesCompatible) || solution$areIndicesCompatible) {
    #   solution$areIndicesCompatible <- areIndicesCompatible
    # }
    return(solution)
    # dual-indexing without pairing
  } else {
    # reminder: unicityConstraint has been set to "none" to simplify the algorithm
    # and also because it is not so important with dual-indexing
    inputNbSamplesPerLane2 <- nrow(indexesList2[[1]])
    compatibleCombinations2 <- vector(mode="list", length=nbLanes)
    k <- 1
    # print('Second while loop: ')  # //--- comment me out
    # print(paste0('nbLanes: ', nbLanes))  # //--- comment me out
    
    while (k <= nbLanes){
      # print(paste0("Start of while loop, k: ", k))
      # print(paste0("length(indexesList) in global.r: ", length(indexesList)))
      # print(paste0("length(indexesList2) in global.r: ", length(indexesList2)))
      if (length(indexesList) == 0 | length(indexesList2) == 0){
        # no available combination of indices anymore, try again!
        return(NULL)
      } else{
        i <- sample(1:length(indexesList), 1, FALSE)
        i2 <- ifelse(i7i5pairing, i, sample(1:length(indexesList2), 1, FALSE))
        # if (areIndexesCompatible(indexesList[[i]], chemistry)){
        # areIndexesCompatible <- ((areIndexesCompatible(indexesList[[i]], chemistry)) & (areIndexesCompatible(indexesList2[[i2]], chemistry)))
        # print('Before is.null(areIndicesCompatible')  # //--- comment me out
        if (is.null(areIndicesCompatible) || areIndicesCompatible) {
          areIndexesCompatible <- ((areIndexesCompatible(indexesList[[i]], chemistry)) & (areIndexesCompatible(indexesList2[[i2]], chemistry)))
        }
        # print('After is.null(areIndicesCompatible')  # //--- comment me out
        compatibleCombinations[[k]] <- indexesList[[i]]
        # if (areIndexesCompatible(indexesList2[[i2]], chemistry)){
        compatibleCombinations2[[k]] <- indexesList2[[i2]]
        k <- k+1
        # } else{
        #   indexesList2 <- indexesList2[-i2]
        # }
        # } else{
          # indexesList <- indexesList[-i]
        # }
      }
      # print(paste0("End of while loop, k: ", k))
    }
    # print('Got past second while loop')  # //--- comment me out
    solution <- data.frame(sample=1:(nbLanes*inputNbSamplesPerLane), 
                           pool=rep(1:nbLanes, each=inputNbSamplesPerLane), 
                           do.call("rbind", compatibleCombinations))
    solution2 <- data.frame(sample=1:(nbLanes*inputNbSamplesPerLane2), 
                            pool=rep(1:nbLanes, each=inputNbSamplesPerLane2), 
                            do.call("rbind", compatibleCombinations2))
    
    # print(paste0('Solution: ', solution))  # //--- comment me out
    
    # remove indices starting with GG before completing the solutions if four-channel chemistry
    if (chemistry == "2"){
      index <- index[!sapply(index$sequence, substr, 1, 2) == "GG",]
      index2 <- index2[!sapply(index2$sequence, substr, 1, 2) == "GG",]
    }
    # only a partial solution has been found, need to complete it
    if (min(multiplexingRate, nrow(index)) > inputNbSamplesPerLane){
      solution <- completeSolution(partialSolution = solution,
                                   index = index,
                                   multiplexingRate = min(multiplexingRate, nrow(index)),
                                   unicityConstraint = unicityConstraint)
    }
    if (min(multiplexingRate, nrow(index2)) > inputNbSamplesPerLane2){
      solution2 <- completeSolution(partialSolution = solution2,
                                    index = index2,
                                    multiplexingRate = min(multiplexingRate, nrow(index2)),
                                    unicityConstraint = unicityConstraint)
    }
    names(solution) <- paste0(names(solution), "1")
    names(solution2) <- paste0(names(solution2), "2")
    # generate all the possible couple of indices
    solution.merged <- merge(solution, solution2, by.x="pool1", by.y="pool2")
    solution <- list()
    for (pool in unique(solution.merged$pool1)){
      # look for the most diverse couple of indices can be slightly difficult
      # we thus encapsulate it into a while() loop
      solution.pool.OK <- FALSE
      while (!solution.pool.OK){
        solution.merged.pool <- solution.merged[which(solution.merged$pool1 == pool),]
        counts.id1 <- numeric(length(unique(solution.merged.pool$id1)))
        names(counts.id1) <- unique(solution.merged.pool$id1)
        counts.id2 <- numeric(length(unique(solution.merged.pool$id2)))
        names(counts.id2) <- unique(solution.merged.pool$id2)
        solution.tmp.pool <- NULL
        for (i in 1:multiplexingRate){
          for (k in sample(1:nrow(solution.merged.pool), nrow(solution.merged.pool), FALSE)){
            counts.id1.k <- counts.id1 + I(names(counts.id1) == solution.merged.pool[k, "id1"])
            counts.id2.k <- counts.id2 + I(names(counts.id2) == solution.merged.pool[k, "id2"])
            maxdiff <- function(x) max(x) - min(x)
            # make sure no index is used n+2 times and another n times
            # i.e. indices must be used a homogeneous number of times
            if (maxdiff(counts.id1.k) <= 1 & maxdiff(counts.id2.k) <= 1){
              counts.id1 <- counts.id1.k
              counts.id2 <- counts.id2.k
              solution.tmp.pool <- rbind(solution.tmp.pool, solution.merged.pool[k,])
              solution.merged.pool <- solution.merged.pool[-k,]
              break
            }
          }
        }
        if (nrow(solution.tmp.pool) == multiplexingRate) solution.pool.OK <- TRUE
      }
      solution[[pool]] <- solution.tmp.pool
    }
    solution <- do.call("rbind", solution)
    solution <- solution[order(solution$pool1, solution$id1, solution$id2),]
    # add scores
    solution$score1 <- unlist(tapply(solution$sequence1, solution$pool, scores))
    solution$score2 <- unlist(tapply(solution$sequence2, solution$pool, scores))
    solution <- data.frame(sample=1:(nbLanes*multiplexingRate),
                           pool=solution$pool1,
                           id1=solution$id1,
                           sequence1=solution$sequence1,
                           color1=solution$color1,
                           score1=solution$score1,
                           weight1=solution$weight1,
                           id2=solution$id2,
                           sequence2=solution$sequence2,
                           color2=solution$color2,
                           score2=solution$score2,
                           weight2=solution$weight2,
                           stringsAsFactors = FALSE)
    
    solution$areIndicesCompatible <- areIndicesCompatible

    return(solution)
  }
}

completeSolution <- function(partialSolution, index, multiplexingRate, unicityConstraint){
  nbSamplesToAdd <- multiplexingRate - nrow(partialSolution)/max(partialSolution$pool) # to each lane
  partialSolution$sample <- NULL
  for (l in unique(partialSolution$pool)){
    if (unicityConstraint == "index"){
      # remove all the indices already used
      index.remaining <- index[!(index$id %in% partialSolution$id),]
    } else{
      # remove the indices already used in the current lane
      index.remaining <- index[!(index$id %in% partialSolution$id[partialSolution$pool==l]),]
    }
    if (nbSamplesToAdd > nrow(index.remaining)) return(NULL) # not enough remaining indices to complete the solution
    
    indexesToAdd <- index.remaining[sample(1:nrow(index.remaining), nbSamplesToAdd, FALSE),]
    indexesToAdd$pool <- l
    partialSolution <- rbind.data.frame(partialSolution, 
                                        indexesToAdd[, c(ncol(indexesToAdd), 1:(ncol(indexesToAdd)-1))])
  }
  finalSolution <- data.frame(sample=1:nrow(partialSolution), 
                              partialSolution[order(partialSolution$pool, partialSolution$id),])
  return(finalSolution)
}

calculateSolutionScore <- function(solution, chemistry, percentage_threshold, distance_threshold) {
  
  solution_color_percentages_list <- calculate_color_percentages_with_weights(solution)

  if (chemistry == "X") {
    percentages_vector <- solution_color_percentages_list$vector_gb
  } else {
    percentages_vector <- solution_color_percentages_list$vector_rg
  }

  # Check condition
  if ("score" %in% names(solution)) {
    # Single score case
    distance_metric <- if (any(solution$score < distance_threshold)) 0 else 1
  } else if (all(c("score1", "score2") %in% names(solution))) {
    # Dual score case
    distance_metric <- if (any(solution$score1 < distance_threshold) || 
                           any(solution$score2 < distance_threshold)) 0 else 1
  } 
  
  # Check condition
  # if (any(solution$score1 < distance_threshold) || any(solution$score2 < distance_threshold)) {
  #   distance_metric <- 0
  # } else {
  #   distance_metric <- 1
  # }
  
  # cat(solution$score1, solution$score2, "\n")
  # print(paste0('distance_metric: ', distance_metric))
  
  count_of_values_below_threshold <- sum(percentages_vector < percentage_threshold)
  rank_score = (mean(percentages_vector) + min(percentages_vector)) / (1 + count_of_values_below_threshold)
  
  # cat(percentages_vector, "\n")
  # print(paste0('count_of_values_below_threshold: ', count_of_values_below_threshold))
  # print(paste0('rank_score: ', rank_score))
  
  return(list(rank_score = rank_score, distance_metric = distance_metric))

  # For now, just return a random integer between 1 and 100
  # return(sample(1:100, 1))
}

findSolution <- function(indexesList, index, indexesList2=NULL, index2=NULL,
                         nbSamples, multiplexingRate, unicityConstraint, nbMaxTrials, 
                         completeLane, selectCompIndexes, chemistry, i7i5pairing){

  percentage_threshold <- 25
  distance_threshold <- 3
  
  # this function run searchOneSolution() nbMaxTrials times until finding a solution based on the parameters defined
  if (!unicityConstraint %in% c("none", "index", "lane")) stop("unicityConstraint parameter must be equal to 'none', 'lane' or 'index'.")
  if (unicityConstraint=="index" & nbSamples > nrow(index)) stop("More samples than available indices: cannot use each index only once.")
  if (nbSamples %% multiplexingRate != 0) stop("Number of samples must be a multiple of the multiplexing rate.")
  if (I(!is.null(index2) | !is.null(indexesList2)) & unicityConstraint != "none") stop("unicityConstraint parameter must be equal to 'none' with dual-indexing.")
  nbLanes <- nbSamples/multiplexingRate
  if (unicityConstraint!="none" & completeLane & selectCompIndexes & length(indexesList) < nbLanes){
    stop("There are only ", length(indexesList), " combinations of compatible indices to fill ", nbLanes, " lanes.")
  }
  
  # Initialize an empty list to store solutions and their scores
  solutions_list <- list()
  last_non_null_solution <- NULL  # Track the last non-null solution
  
  nbTrials <- 1
  
  # Need to compare multiplexingRate to total number of indices (not number of samples)
  # print(paste0('nbSamples: ', nbSamples))
  # print(paste0('multiplexingRate: ', multiplexingRate))
  # print(paste0('nrow(index): ', nrow(index)))
  # print(paste0('length(indexesList): ', length(indexesList)))
  # print(paste0('nbLanes: ', nbLanes))
  
  if (nrow(index) == multiplexingRate) {
    nbMaxTrials <- 1  # Only need a single trial if there is only one pool
  }
  while (nbTrials <= nbMaxTrials){
    solution <- searchOneSolution(indexesList = indexesList,
                                  index = index,
                                  indexesList2 = indexesList2,
                                  index2 = index2,
                                  nbLanes = nbLanes,
                                  multiplexingRate = multiplexingRate,
                                  unicityConstraint = unicityConstraint,
                                  chemistry = chemistry,
                                  i7i5pairing = i7i5pairing)
    
    if (!is.null(solution)){
      last_non_null_solution <- solution  # Update the last non-null solution
      # check solution only for single-indexing as dual-indexing has no constraint
      if (is.null(indexesList2) & is.null(index2) & !i7i5pairing){
        if (checkProposedSolution(solution = solution, unicityConstraint = unicityConstraint, chemistry = chemistry)) {
          # Store valid solution with its score
          calculateSolutionScore_result = calculateSolutionScore(solution, chemistry, percentage_threshold, distance_threshold)
          solutions_list[[nbTrials]] <- list(
            solution = solution,
            score = calculateSolutionScore_result$rank_score,
            distance_metric = calculateSolutionScore_result$distance_metric
          )
        }
      } else{
        if (checkProposedSolution2(solution = solution, chemistry = chemistry)) {
          # Store valid solution with its score
          calculateSolutionScore_result = calculateSolutionScore(solution, chemistry, percentage_threshold, distance_threshold)
          solutions_list[[nbTrials]] <- list(
            solution = solution,
            score = calculateSolutionScore_result$rank_score,
            distance_metric = calculateSolutionScore_result$distance_metric
          )
        }
      }
    }
    nbTrials <- nbTrials + 1
  }

  # First, check if we have any solutions in our list
  valid_solutions <- solutions_list[!sapply(solutions_list, is.null)]
  
  # print(paste0("solutions_list: ", solutions_list))
  # print(paste0("valid_solutions: ", valid_solutions))
  # print(paste0("length(valid_solutions): ", length(valid_solutions)))
  # print(paste0("last_non_null_solution: ", last_non_null_solution))
  
  # If solutions were found, return the one with highest score
  # Otherwise, return the last solution and show a warning message
  if (length(valid_solutions) > 0) {
    # Only create solutions_df if we have valid solutions
    solutions_df <- do.call(rbind, lapply(valid_solutions, function(x) {
      # print(str(x))
      # print('\nfindSolution1:')
      # test_GC_counts(c(x$solution$color1, x$solution$color2))
      data.frame(
        solution = I(list(x$solution)),
        score = as.numeric(x$score),
        distance_metric = as.numeric(x$distance_metric)
      )
    }))
    
    # solutions_df <- tryCatch({
    #   do.call(rbind, lapply(valid_solutions, function(x) {
    #     print(str(x))
    #     data.frame(
    #       solution = I(list(if (!is.null(x$solution)) x$solution else NA)),
    #       score = if (!is.null(x$score) && length(x$score) > 0) x$score else NA,
    #       distance_metric = if (!is.null(x$distance_metric) && length(x$distance_metric) > 0) x$distance_metric else NA
    #     )
    #   }))
    # }, error = function(e) {
    #   print(paste("Error in solutions_df:", e$message))  # Console
    #   showNotification(paste("Error:", e$message), type = "error")  # UI
    #   return(NULL)
    # })
    
    # best_solution <- solutions_df[which.max(solutions_df$score), "solution"][[1]]
    # First, try to find the best among distance_metric == 1
    preferred_solutions <- subset(solutions_df, distance_metric == 1)
    
    # print(paste0("solutions_df$score: ", solutions_df$score))
    # print(paste0("solutions_df$distance_metric: ", solutions_df$distance_metric))
    # print(paste0("preferred_solutions$score: ", preferred_solutions$score))
    # print(paste0("preferred_solutions$distance_metric: ", preferred_solutions$distance_metric))
    
    if (nrow(preferred_solutions) > 0) {
      best_solution <- preferred_solutions[which.max(preferred_solutions$score), "solution"][[1]]
      # print(paste0("best_solution: ", which.max(preferred_solutions$score)))
      # print(paste0("best_solution: ", best_solution))
      # cat(best_solution$color1, best_solution$color2, "\n")
      
      # print('findSolution2:')
      # test_GC_counts(c(best_solution$color1, best_solution$color2))

    } else {
      # Fallback: use distance_metric == 0
      fallback_solutions <- subset(solutions_df, distance_metric == 0)
      best_solution <- fallback_solutions[which.max(fallback_solutions$score), "solution"][[1]]
    }
    
    return(best_solution)
  } else {
    # print(paste0("last_non_null_solution: ", last_non_null_solution))
    return(last_non_null_solution)
  }
  
  print("Didn't find a solution in findSolution in global.r")
  stop(paste("No solution found after", nbMaxTrials, "trials using these parameters, you can modify them or increase the number of trials."))
}

# single-indexing solution checking
checkProposedSolution <- function(solution, unicityConstraint, chemistry){
  # indices not unique
  if (unicityConstraint=="index" & any(duplicated(solution$id))){
    stop("The solution proposed uses some indices several times. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # indices not unique within a pool/lane
  if (any(sapply(split(solution$id, solution$pool), function(x) any(duplicated(x))))){
    stop("The solution proposed uses some indices several times within a lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # several pools/lanes with the same combination of indices
  if (unicityConstraint=="lane" & any(duplicated(split(solution$id, solution$pool)))){
    stop("The solution proposed uses some combinations of indices several times. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # different number of samples on the pools/lanes
  if (length(table(table(solution$pool))) > 1){
    stop("The solution proposed uses different numbers of samples per lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # indices not compatible in at least one pool/lane
  # print('Came from checkProposedSolution')
  if (any(!sapply(split(solution, solution$pool), areIndexesCompatible, chemistry))){
    # stop("The solution proposed uses incompatible indices Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
    return(FALSE)
  }
  # two-channel chemistry and indices starting with GG
  if (chemistry == "2" && any(sapply(solution$sequence, substr, 1, 2) == "GG")){
    stop("Indices starting with GG are not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # one-channel chemistry and all indices starting with GG
  if (chemistry == "1" && any(sapply(lapply(split(solution$sequence, solution$pool), substr, 1, 2), function(x) all(x=="GG")))){
    stop("Having all the indices of a pool starting with GG is not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  return(TRUE)
}

# dual-indexing solution checking
checkProposedSolution2 <- function(solution, chemistry){
  # indices not unique within a pool/lane
  if (any(sapply(split(paste(solution$id1, solution$id2, sep="-"), solution$pool), function(x) any(duplicated(x))))){
    stop("The solution proposed uses some dual-index combinations several times within a lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # different number of samples on the pools/lanes
  if (length(table(table(solution$pool))) > 1){
    stop("The solution proposed uses different numbers of samples per lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # indices not compatible in at least one pool/lane
  # print('Came from checkProposedSolution2')
  sol1 <- data.frame(id=solution$id1, sequence=solution$sequence1, color=solution$color1, stringsAsFactors = FALSE)
  if (any(!sapply(split(sol1, solution$pool), areIndexesCompatible, chemistry))){
    # stop("The solution proposed uses incompatible indices. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
    return(FALSE)
  }
  sol2 <- data.frame(id=solution$id2, sequence=solution$sequence2, color=solution$color2, stringsAsFactors = FALSE)
  if (any(!sapply(split(sol2, solution$pool), areIndexesCompatible, chemistry))){
    # stop("The solution proposed uses incompatible indices. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
    return(FALSE)
  }
  # two-channel chemistry and indices starting with GG
  if (chemistry == "2" && any(sapply(solution$sequence1, substr, 1, 2) == "GG")){
    stop("Indices starting with GG are not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  if (chemistry == "2" && any(sapply(solution$sequence2, substr, 1, 2) == "GG")){
    stop("Indices starting with GG are not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # one-channel chemistry and all indices starting with GG
  if (chemistry == "1" && any(sapply(lapply(split(solution$sequence1, solution$pool), substr, 1, 2), function(x) all(x=="GG")))){
    stop("Having all the indices of a pool starting with GG is not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  if (chemistry == "1" && any(sapply(lapply(split(solution$sequence2, solution$pool), substr, 1, 2), function(x) all(x=="GG")))){
    stop("Having all the indices of a pool starting with GG is not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  return(TRUE)
}

heatmapindex <- function(solution, selectedRows){
  
  percentage_threshold <- 25
  
  # print('heatmapindex')
  # print(paste0('selectedRows:', selectedRows$id))
  
  lengthSequence <- nchar(solution$sequence[1])
  lengthSequence1 <- nchar(solution$sequence1[1])
  lengthSequence2 <- nchar(solution$sequence2[1])
  
  # print(paste0('solution: ', solution))
  solution_color_percentages_list <- calculate_color_percentages_with_weights(solution)
  solution_color_percentages <- solution_color_percentages_list$solution_color_percentages
  # print(paste0('solution_color_percentages: ', colnames(solution_color_percentages)))
  
  if ("id1" %in% names(solution)){
    # dual indexing, need to format the solution data.frame as for single-indexing
    solution <- data.frame(sample=solution$sample, pool=solution$pool,
                           id=ifelse(solution$id1 == solution$id2, 
                                     paste(solution$id1, paste("(",round(solution$weight1),")", sep="")),
                                     paste(solution$id1, solution$id2, sep="-")),
                           sequence=paste(solution$sequence1, solution$sequence2),
                           color=paste(solution$color1, solution$color2),
                           stringsAsFactors = FALSE)
  }
  
  # Get unique pools
  pools <- sort(as.numeric(unique(str_extract(colnames(solution_color_percentages), "^pool(\\d+)") %>% 
                                    str_replace("pool", ""))))
  
  num_pools <- length(pools)
  numberPositions <- 1 + 
    sum(as.numeric(lengthSequence %||% 0), 
        as.numeric(lengthSequence1 %||% 0), 
        as.numeric(lengthSequence2 %||% 0), 
        na.rm = TRUE)
  
  lengthSequence_first_col <- sum(as.numeric(lengthSequence %||% 0), as.numeric(lengthSequence1 %||% 0), na.rm = TRUE)
  # print(paste0('lengthSequence_first_col: ', lengthSequence_first_col))
  # print(paste0('numberPositions: ', numberPositions))
  
  color_percentages.pool <- matrix(" ", nrow = 4 * num_pools, ncol = numberPositions)
  
  for(p in seq_along(pools)) {
    pool <- pools[p]
    row_indices <- (4*(p-1) + 1):(4*p)
    
    # Process Color1 (positions 1-lengthSequence1)
    for(pos in 1:lengthSequence_first_col) {
      # Get columns for this position and pool
      r_col <- grep(paste0("pool", pool, "_pos", pos, ".*_position",pos,"_R$"), colnames(solution_color_percentages), value = TRUE)
      g_col <- grep(paste0("pool", pool, "_pos", pos, ".*_position",pos,"_G$"), colnames(solution_color_percentages), value = TRUE)
      b_col <- grep(paste0("pool", pool, "_pos", pos, ".*_position",pos,"_B$"), colnames(solution_color_percentages), value = TRUE)
      dash_col <- grep(paste0("pool", pool, "_pos", pos, ".*_position",pos,"_-$"), colnames(solution_color_percentages), value = TRUE)
      
      # print(paste0('pos: ', pos))
      # print(paste0('r_col: ', r_col))
      
      # Fill values for Color1 (first set of columns)
      if(length(r_col) > 0) {
        val <- solution_color_percentages[1, r_col[1]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[1], pos] <- paste0(round(val), "%")
          # numeric_values[row_indices[1], pos] <- val
        }
      }
      if(length(g_col) > 0) {
        val <- solution_color_percentages[1, g_col[1]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[2], pos] <- paste0(round(val), "%")
          # numeric_values[row_indices[2], pos] <- val
        }
      }
      if(length(b_col) > 0) {
        val <- solution_color_percentages[1, b_col[1]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[3], pos] <- paste0(round(val), "%")
          # numeric_values[row_indices[3], pos] <- val
        }
      }
      if(length(dash_col) > 0) {
        val <- solution_color_percentages[1, dash_col[1]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[4], pos] <- paste0(round(val), "%")
          # numeric_values[row_indices[4], pos] <- val
        }
      }
      
      # Fill values for Color2 (second set of columns)
      if(length(r_col) > 1) {
        val <- solution_color_percentages[1, r_col[2]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[1], pos + lengthSequence_first_col + 1] <- paste0(round(val), "%")
          # numeric_values[row_indices[1], pos + 8] <- val
        }
      }
      if(length(g_col) > 1) {
        val <- solution_color_percentages[1, g_col[2]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[2], pos + lengthSequence_first_col + 1] <- paste0(round(val), "%")
          # numeric_values[row_indices[2], pos + 8] <- val
        }
      }
      if(length(b_col) > 1) {
        val <- solution_color_percentages[1, b_col[2]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[3], pos + lengthSequence_first_col + 1] <- paste0(round(val), "%")
          # numeric_values[row_indices[3], pos + 8] <- val
        }
      }
      if(length(dash_col) > 1) {
        val <- solution_color_percentages[1, dash_col[2]]
        if(!is.na(val) && is.numeric(val)) {
          color_percentages.pool[row_indices[4], pos + lengthSequence_first_col + 1] <- paste0(round(val), "%")
          # numeric_values[row_indices[4], pos + 8] <- val
        }
      }
    }
  }
  
  # Set row and column names
  row_names <- c()
  for(pool in pools) {
    row_names <- c(row_names, paste0("Pool", pool, c("_R", "_G", "_B", "_-")))
  }
  rownames(color_percentages.pool) <- row_names
  colnames(color_percentages.pool) <- paste0("Pos", 1:numberPositions)
  
  # Find rows that have all zeros
  # non_zero_rows <- which(rowSums(numeric_values) > 0)
  
  # Keep only non-zero rows
  # color_percentages.pool <- color_percentages.pool[non_zero_rows, ]
  
  # Update row names to reflect only the kept rows
  # row_names <- rownames(color_percentages.pool)
  # row_names <- row_names[non_zero_rows]
  # rownames(color_percentages.pool) <- row_names
  
  splitsol <- split(solution, solution$pool)
  # build a matrix containing the index bases and NAs to separate the pools/lanes
  seqmat <- lapply(splitsol, function(splitsol.pool) do.call("rbind", strsplit(splitsol.pool$sequence, "")))
  # seqmat <- do.call("rbind", lapply(seqmat, function(seqmat.pool) rbind(seqmat.pool, rep(NA, ncol(seqmat.pool)))))
  
  seqmat <- do.call("rbind", lapply(1:length(seqmat), function(i) {
   row_indices <- (4*(i-1) + 1):(4*i)
   rbind(
     color_percentages.pool[row_indices[1],],  # Add corresponding row of percentages
     color_percentages.pool[row_indices[2],],  # Add corresponding row of percentages
     color_percentages.pool[row_indices[3],],  # Add corresponding row of percentages
     color_percentages.pool[row_indices[4],],  # Add corresponding row of percentages
     seqmat[[i]],  # Original sequence matrix for this pool
     rep(NA, ncol(seqmat[[i]]))  # Separator row
   )
  }))
  
  seqmat <- seqmat[-nrow(seqmat),]
  # build a matrix containing the colors and NAs to separate the pools/lanes
  seqcol <- lapply(splitsol, function(splitsol.pool) do.call("rbind", strsplit(splitsol.pool$color, "")))
  # seqcol <- do.call("rbind", lapply(seqcol, function(seqcol.pool) rbind(seqcol.pool, rep(NA, ncol(seqcol.pool)))))
  
  seqcol <- do.call("rbind", lapply(1:length(seqcol), function(i) {
    rbind(
      rep(NA, ncol(seqcol[[i]])),  # Add row for percentages
      rep(NA, ncol(seqcol[[i]])),  # Add row for percentages
      rep(NA, ncol(seqcol[[i]])),  # Add row for percentages
      rep(NA, ncol(seqcol[[i]])),  # Add row for percentages
      seqcol[[i]],  # Original sequence matrix for this pool
      rep(NA, ncol(seqcol[[i]]))  # Separator row
    )
  }))
  seqcol <- seqcol[-nrow(seqcol),]
  
  # plot
  par(mar=c(2, 6, 1, 6) + 0.1, xpd=TRUE, xaxs="i", yaxs="i") # ,
  #    pin=c(8, 23.5))  # Width=6 inches, Height=calculated value
  plot.new()
  # plot.window(xlim=c(-1.5, ncol(seqmat)+1.5), ylim=c(-2, nrow(seqmat)))
  plot.window(xlim=c(-1.5, ncol(seqmat)+1.5), ylim=c(-2, nrow(seqmat)))
  
  for (i in 1:nrow(seqmat)){
    # Check if entire row contains only zeros
    row_values <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", seqmat[i,])))
    # print(paste0('row_values: ', row_values))
    is_all_zeros <- all(row_values[!is.na(row_values)] == 0)
    
    for (j in 1:ncol(seqmat)){
      rect(xleft=j-1, ybottom=nrow(seqmat)-(i-1), xright=j, ytop=(nrow(seqmat)-(i)),
           col=switch(seqcol[i,j], "R"="orangered", "G"="darkseagreen4", "-"="white", "O"="orange", "NA"="white", " "="white", "B"="blue", "C"="cyan"),
           border=switch(seqcol[i,j], "R"="orangered", "G"="darkseagreen4", "-"="white", "O"="orange", "NA"="white", " "="white", "B"="blue", "C"="cyan"))
           # border=switch(seqcol[i,j], "R"="orangered", "orangered"))
      # Extract number from percentage string (if it exists)
      value <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", seqmat[i,j])))
      
      # Determine text color based on value
      text_color <- "black"  # default color
      if (!is.na(value)) {  # check if it's actually a number
        if (value == 0 && is_all_zeros) {
          text_color <- "lightgray"
        } else if (value < percentage_threshold) {  # adjust threshold as needed
          text_color <- "red"
        }
      }
      text(x=j-0.5, y=nrow(seqmat)-(i-0.5), labels=seqmat[i,j], col=text_color)
    }
  }
  
  # extract and print sample ids
  tmp <- unlist(lapply(splitsol, function(sol) c("R", "G", "B", "-", sol$sample, NA)))
  text(x=-0.25, y=(nrow(seqmat):1 - 0.5), labels=tmp[-length(tmp)], pos=2)
  text(x=-1.5, y=nrow(seqmat)/2, labels="Sample", srt=90)
  
  # extract and print index ids
  tmp <- unlist(lapply(splitsol, function(sol) c(NA, NA, NA, NA, sol$id, NA)))
  # text(x=ncol(seqmat)+0.25, y=(nrow(seqmat):1 - 0.5), labels=tmp[-length(tmp)], pos=4)
  labs  <- tmp[-length(tmp)]
  
  ## bare ID = the part after “: ” and before any space / parentheses
  ids   <- sub(" .*", "", sub("^.*: ", "", labs))
  
  ## character vector of colours: red for selected IDs, black otherwise
  cols  <- ifelse(ids %in% selectedRows$id, "deeppink", "black")
  
  ## 1 ── build a font vector: 2 = bold, 1 = plain
  font_vec <- ifelse(ids %in% selectedRows$id, 2, 1)   # ids is bare “P4_A01” etc.
  
  ## keep colours vector aligned with labels (NA labels get NA colours)
  cols[is.na(labs)] <- NA
  
  ## draw the labels
  text(
    x     = ncol(seqmat) + 0.25,
    y     = (nrow(seqmat):1) - 0.5,
    labels = labs,
    pos    = 4,          # left-align at (x, y)
    col    = cols,        # <- red for selected, black otherwise
    font   = font_vec      # bold if selected, plain otherwise
  )
  
  # ## --------------------------------------------------------------------
  # ## PREP
  # ## --------------------------------------------------------------------
  # tmp <- unlist(lapply(splitsol, function(sol) c(NA, NA, NA, NA, sol$id, NA)))
  # labels <- tmp[-length(tmp)]                           # one row per cell
  # sel_id <- selectedRows$id                            # data-frame column
  # 
  # ## helper: pull the “P4_A01” part out of “tmo: P4_A01 (206)”
  # get_id <- function(z) sub(" .*", "", sub("^.*: ", "", z))
  # 
  # ## ----------------------------------------------------------------------
  # ## NEW ─ filter out the placeholder rows  -------------------------------
  # ## ----------------------------------------------------------------------
  # keep      <- !is.na(labels)                       # TRUE only for real labels
  # labels    <- labels[keep]                         # drop the NAs
  # is_sel    <- get_id(labels) %in% sel_id
  # 
  # fg_col    <- ifelse(is_sel, "white", "black")
  # bg_col    <- ifelse(is_sel, "grey25", "white")
  # 
  # ## Coordinates that match the *kept* rows -------------------------------
  # y_all     <- nrow(seqmat):1 - 0.5                 # original y-seq
  # x_all     <- rep(ncol(seqmat) + 2.0, length(y_all))
  # 
  # x <- x_all[keep]                                  # keep same rows only
  # y <- y_all[keep]
  # 
  # print(paste0('tmp: ', tmp))
  # print(paste0('keep: ', keep))
  # print(paste0('labels: ', labels))
  # print(paste0('sel_id: ', sel_id))
  # print(paste0('is_sel: ', is_sel))
  # print(paste0('x: ', x))
  # print(paste0('y: ', y))
  # 
  # ## ----------------------------------------------------------------------
  # ## Draw the boxed labels: all rows now line up perfectly
  # ## ----------------------------------------------------------------------
  # boxed.labels(
  #   x, y, labels,
  #   col    = fg_col,
  #   bg     = bg_col,
  #   adj    = c(0, 0.5),       # left-align at x
  #   xpad   = 1.2,               # padding (char widths)
  #   ypad   = 1.8
  # )  
}

# border = NA,                # no rectangle outline

# Function to calculate the percentage of each character in each position of color1 and color2
calculate_color_percentages <- function(solution) {
  possible_characters <- c("R", "G", "B", "-", "C", "O")
  
  # Helper function to split characters and pad with NA if lengths are unequal
  split_and_pad <- function(color_str, max_length) {
    chars <- unlist(strsplit(color_str, ""))
    if (length(chars) < max_length) {
      chars <- c(chars, rep(NA, max_length - length(chars)))
    }
    return(chars)
  }
  
  # Split color1 and color2 into separate columns for each position
  # print(paste0('nchar(solution$color1): ', nchar(solution$color1)))
  # print(paste0('nchar(solution$color2): ', nchar(solution$color2)))
  max_length <- max(nchar(solution$color), nchar(solution$color1), nchar(solution$color2))
  
  expanded_colors <- solution %>%
    rowwise() %>%
    mutate(
      color1_chars = list(split_and_pad(color1, max_length)),
      color2_chars = list(split_and_pad(color2, max_length))
    ) %>%
    ungroup() %>%
    select(color1_chars, color2_chars) %>%
    unnest_wider(color1_chars, names_sep = "_") %>%
    unnest_wider(color2_chars, names_sep = "_")
  
  # Function to calculate adjusted percentages for each position
  calculate_adjusted_percentages <- function(df, col_prefix, max_length) {
    result_list <- list()
    for (i in 1:max_length) {
      col_name <- paste0(col_prefix, "_", i)
      
      position_percentages <- df %>%
        count(!!sym(col_name)) %>%
        filter(!is.na(!!sym(col_name))) %>%
        mutate(percentage = n / sum(n) * 100) %>%
        select(-n) %>%
        pivot_wider(names_from = !!sym(col_name), values_from = percentage, values_fill = list(percentage = 0))
      
      # Ensure all possible characters are present
      for (character in possible_characters) {
        if (!(character %in% colnames(position_percentages))) {
          position_percentages[[character]] <- 0
        }
      }
      
      # Adjust based on presence of 'C' and 'O'
      if ("C" %in% colnames(position_percentages)) {
        position_percentages <- position_percentages %>%
          mutate(
            B = B + C,
            G = G + C,
            C = 0
          )
      }
      if ("O" %in% colnames(position_percentages)) {
        position_percentages <- position_percentages %>%
          mutate(
            R = R + O,
            G = G + O,
            O = 0
          )
      }
      
      # Add position indicator to column names
      position_colnames <- paste0(col_prefix, "_position", i, "_", colnames(position_percentages))
      colnames(position_percentages) <- position_colnames
      result_list[[i]] <- position_percentages
    }
    return(do.call(cbind, result_list))
  }
  
  # Calculate percentages for color1 and color2
  color1_percentages <- calculate_adjusted_percentages(expanded_colors, "color1_chars", max_length)
  color2_percentages <- calculate_adjusted_percentages(expanded_colors, "color2_chars", max_length)
  
  # Combine all percentages into a single dataframe
  solution_color_percentages <- bind_cols(color1_percentages, color2_percentages)
  
  return(solution_color_percentages)
}

# Function to calculate the percentage of each character in each position of color1 and color2
calculate_color_percentages_with_weights <- function(solution) {
  possible_characters <- c("R", "G", "B", "-", "C", "O")
  
  # Helper function to split characters and pad with NA if lengths are unequal
  split_and_pad <- function(color_str, max_length) {
    chars <- unlist(strsplit(color_str, ""))
    if (length(chars) < max_length) {
      chars <- c(chars, rep(NA, max_length - length(chars)))
    }
    return(chars)
  }
  
  # Max length of color sequences
  # print(paste0('names(solution): ', names(solution)))
  # print(paste0('solution$color: ', solution$color))
  max_length <- max(nchar(solution$color), nchar(solution$color1), nchar(solution$color2))
  # print(paste0('nchar(solution$color): ', max_length))
  
  # Expand color1 and color2 into separate columns for each position
  expanded_colors <- solution %>%
    rowwise() %>%
    {
      if ("color" %in% names(.)) {
        # Case for single color and weight
        mutate(., 
               color_chars = list(split_and_pad(color, max_length))) %>%
          ungroup() %>%
          select(pool, color_chars, weight) %>%
          unnest_wider(color_chars, names_sep = "_")
      } else {
        # Case for color1, color2, weight1, weight2
        mutate(.,
               color1_chars = list(split_and_pad(color1, max_length)),
               color2_chars = list(split_and_pad(color2, max_length))) %>%
          ungroup() %>%
          select(pool, color1_chars, color2_chars, weight1, weight2) %>%
          unnest_wider(color1_chars, names_sep = "_") %>%
          unnest_wider(color2_chars, names_sep = "_")
      }
    }
  
  # Expand color1 and color2 into separate columns for each position
  # expanded_colors <- solution %>%
  #   rowwise() %>%
  #   mutate(
  #     color1_chars = list(split_and_pad(color1, max_length)),
  #     color2_chars = list(split_and_pad(color2, max_length))
  #   ) %>%
  #   ungroup() %>%
  #   select(pool, color1_chars, color2_chars, weight1, weight2) %>%  # Added pool to select
  #   unnest_wider(color1_chars, names_sep = "_") %>%
  #   unnest_wider(color2_chars, names_sep = "_")
  
  # Function to calculate weighted adjusted percentages for each position
  calculate_adjusted_percentages <- function(df, col_prefix, weight_col, max_length) {
    result_list <- list()
    result_list_R_G <- list()
    result_list_G_B <- list()
    
    # Process each pool separately
    pools <- unique(df$pool)
    
    for (current_pool in pools) {
      pool_df <- df %>% filter(pool == current_pool)
      
      for (i in 1:max_length) {
        col_name <- paste0(col_prefix, "_", i)
        
        # Calculate weighted count for each character within this pool
        position_percentages <- pool_df %>%
          count(!!sym(col_name), wt = !!sym(weight_col)) %>%
          filter(!is.na(!!sym(col_name))) %>%
          mutate(percentage = n / sum(n) * 100) %>%
          select(-n) %>%
          pivot_wider(names_from = !!sym(col_name), values_from = percentage, values_fill = list(percentage = 0))
        
        # Ensure all possible characters are present
        for (character in possible_characters) {
          if (!(character %in% colnames(position_percentages))) {
            position_percentages[[character]] <- 0
          }
        }
        
        # Adjust based on presence of 'C' and 'O'
        if ("C" %in% colnames(position_percentages)) {
          position_percentages <- position_percentages %>%
            mutate(
              B = B + C,
              G = G + C
            )
        }
        if ("O" %in% colnames(position_percentages)) {
          position_percentages <- position_percentages %>%
            mutate(
              R = R + O,
              G = G + O
            )
        }
        
        # Remove C and O from the results
        position_percentages <- position_percentages %>%
          select(-C, -O)
        
        # Set up position percentages for R and G only
        position_percentages_R_G <- position_percentages %>%
          select(R, G)
        
        # Set up position percentages for G and B only
        position_percentages_G_B <- position_percentages %>%
          select(G, B)
        
        # Add pool and position indicator to column names
        position_colnames <- paste0(col_prefix, "_pool", current_pool, "_position", i, "_", colnames(position_percentages))
        colnames(position_percentages) <- position_colnames
        
        # print(paste0('position_colnames: ', position_colnames))
        # print(paste0('position_percentages: ', position_percentages))
        
        # Store results
        result_key <- paste0("pool", current_pool, "_pos", i)
        result_list[[result_key]] <- position_percentages
        
        # Add pool and position indicator to column names for R/G chemistrys
        position_colnames_R_G <- paste0(col_prefix, "_pool", current_pool, "_position", i, "_", colnames(position_percentages_R_G))
        colnames(position_percentages_R_G) <- position_colnames_R_G
        
        # print(paste0('position_colnames_R_G: ', position_colnames_R_G))
        # print(paste0('position_percentages_R_G: ', position_percentages_R_G))
        
        # Store results
        result_key_R_G <- paste0("pool", current_pool, "_pos", i)
        result_list_R_G[[result_key_R_G]] <- position_percentages_R_G
        
        # Add pool and position indicator to column names for G/B chemistrys
        position_colnames_G_B <- paste0(col_prefix, "_pool", current_pool, "_position", i, "_", colnames(position_percentages_G_B))
        colnames(position_percentages_G_B) <- position_colnames_G_B
        
        # print(paste0('position_colnames_G_B: ', position_colnames_G_B))
        # print(paste0('position_percentages_G_B: ', position_percentages_G_B))
        
        # Store results
        result_key_G_B <- paste0("pool", current_pool, "_pos", i)
        result_list_G_B[[result_key_G_B]] <- position_percentages_G_B
      }
    }

    result_all <- do.call(cbind, result_list)
    # result_rg <- do.call(cbind, result_list_R_G)
    # result_gb <- do.call(cbind, result_list_G_B)
    result_vector_rg <- unlist(result_list_R_G)
    result_vector_gb <- unlist(result_list_G_B)
    
    # print(paste0('result_list_G_B: ', result_list_G_B))  
    # cat(result_vector_gb, "\n")
    
    return(list(result_all = result_all, result_vector_rg = result_vector_rg, result_vector_gb = result_vector_gb))
    
    # return(do.call(cbind, result_list))
  }
  
  # Calculate percentages for color1 and color2 using dynamic weight columns
  # print(paste0('names(solution): ', names(solution)))
  if ("color" %in% names(solution)) {
    adjusted_percentages <- calculate_adjusted_percentages(expanded_colors, "color_chars", "weight", max_length)
    color_percentages <- adjusted_percentages$result_all
    color_percentages_vector_rg <- adjusted_percentages$result_vector_rg
    color_percentages_vector_gb <- adjusted_percentages$result_vector_gb

    # Combine all percentages into a single dataframe
    solution_color_percentages <- color_percentages
    solution_color_percentages_vector_rg <- color_percentages_vector_rg
    solution_color_percentages_vector_gb <- color_percentages_vector_gb
    
    # print_color_percentages(solution_color_percentages)
  } else {
    adjusted_percentages <- calculate_adjusted_percentages(expanded_colors, "color1_chars", "weight1", max_length)
    color1_percentages <- adjusted_percentages$result_all
    color1_percentages_vector_rg <- adjusted_percentages$result_vector_rg
    color1_percentages_vector_gb <- adjusted_percentages$result_vector_gb
    adjusted_percentages <- calculate_adjusted_percentages(expanded_colors, "color2_chars", "weight2", max_length)
    color2_percentages <- adjusted_percentages$result_all
    color2_percentages_vector_rg <- adjusted_percentages$result_vector_rg
    color2_percentages_vector_gb <- adjusted_percentages$result_vector_gb

    # print(paste0('color1_percentages: ', color1_percentages))

    # Combine all percentages into a single dataframe
    solution_color_percentages <- bind_cols(color1_percentages, color2_percentages)
    solution_color_percentages_vector_rg <- c(color1_percentages_vector_rg, color2_percentages_vector_rg)
    solution_color_percentages_vector_gb <- c(color1_percentages_vector_gb, color2_percentages_vector_gb)
  
    # print_color_percentages(solution_color_percentages)
    # cat(solution_color_percentages_vector_gb, "\n")
  }
  
  # print(paste0('solution_color_percentages: ', solution_color_percentages))
  return(list(solution_color_percentages = solution_color_percentages, vector_rg = solution_color_percentages_vector_rg, vector_gb = solution_color_percentages_vector_gb))
}

# Format and print the results
print_color_percentages <- function(solution_color_percentages) {
  # Get all unique pools from column names
  pool_numbers <- unique(gsub(".*_pool(\\d+)_position.*", "\\1", 
                              colnames(solution_color_percentages)))
  
  # Get the number of positions
  positions <- unique(gsub(".*_position(\\d+)_.*", "\\1", 
                           colnames(solution_color_percentages)))
  
  # For each pool and color type (color1 and color2)
  for (pool in pool_numbers) {
    cat("\n=== Pool", pool, "===\n")
    
    # Process color1
    cat("\nColor 1:\n")
    cat("Position  |   R    |   G    |   B    |   -    |\n")
    cat("----------+--------+--------+--------+--------+\n")
    
    for (pos in positions) {
      # Extract relevant columns for this position and pool
      cols <- grep(paste0("color1_chars_pool", pool, "_position", pos, "_[RGB-]$"), 
                   colnames(solution_color_percentages), value = TRUE)
      if (length(cols) > 0) {
        values <- solution_color_percentages[1, cols]
        cat(sprintf("%8s  | %5.1f%% | %5.1f%% | %5.1f%% | %5.1f%% |\n",
                    pos,
                    values[[grep("_R$", cols)]],
                    values[[grep("_G$", cols)]],
                    values[[grep("_B$", cols)]],
                    values[[grep("_-$", cols)]]
        ))
      }
    }
    
    # Process color2
    cat("\nColor 2:\n")
    cat("Position  |   R    |   G    |   B    |   -    |\n")
    cat("----------+--------+--------+--------+--------+\n")
    
    for (pos in positions) {
      # Extract relevant columns for this position and pool
      cols <- grep(paste0("color2_chars_pool", pool, "_position", pos, "_[RGB-]$"), 
                   colnames(solution_color_percentages), value = TRUE)
      if (length(cols) > 0) {
        values <- solution_color_percentages[1, cols]
        cat(sprintf("%8s  | %5.1f%% | %5.1f%% | %5.1f%% | %5.1f%% |\n",
                    pos,
                    values[[grep("_R$", cols)]],
                    values[[grep("_G$", cols)]],
                    values[[grep("_B$", cols)]],
                    values[[grep("_-$", cols)]]
        ))
      }
    }
  }
}

# Function to convert percentages into the specified format for output
convert_to_formatted_output <- function(solution_color_percentages) {
  # Reshape and format the dataframe for output
  long_df <- solution_color_percentages %>%
    pivot_longer(cols = everything(), names_to = "column", values_to = "percentage") %>%
    separate(column, into = c("Color", "Position_Char"), sep = "_position") %>%
    mutate(
      Position = str_extract(Position_Char, "\\d+"),
      Character = str_remove(Position_Char, "\\d+_"),
      percentage = paste0(round(percentage, 1), "%")
    ) %>%
    select(Color, Position, Character, percentage) %>%
    rename(Column = Color, Percentage = percentage) %>%
    mutate(Column = str_remove(Column, "_chars")) %>%
    arrange(Column, Position, Character)
  
  return(long_df)
}

# Function to convert percentages into the specified wide format for output
convert_to_formatted_output_wide <- function(solution_color_percentages) {
  possible_characters <- c("R", "G", "B", "-")
  
  # Extract pool numbers from column names
  pools <- sort(as.numeric(unique(str_extract(colnames(solution_color_percentages), "^pool(\\d+)") %>% 
                                    str_replace("pool", ""))))
  
  # Create an empty list to store dataframes for each pool
  pool_dfs <- list()
  
  for(pool in pools) {
    # Filter columns for current pool
    pool_cols <- grep(paste0("^pool", pool, "_"), colnames(solution_color_percentages), value = TRUE)
    pool_data <- solution_color_percentages[, pool_cols]
    
    # Create long format for current pool
    long_df <- pool_data %>%
      pivot_longer(cols = everything(), names_to = "column", values_to = "percentage") %>%
      mutate(
        Pool = pool,  # Add pool number
        Color = str_extract(column, "color\\d"),
        Position = as.numeric(str_extract(column, "position(\\d+)") %>% str_replace("position", "")),
        Character = str_extract(column, "[RGB-]$"),
        Percentage = paste0(round(percentage, 1), "%")
      ) %>%
      select(Pool, Color, Position, Character, Percentage)
    
    # Convert to wide format
    wide_df <- long_df %>%
      filter(Character %in% possible_characters) %>%
      pivot_wider(names_from = Character, values_from = Percentage, values_fill = list(Percentage = "0%")) %>%
      arrange(Color, Position) %>%
      select(Pool, Color, Position, all_of(possible_characters))
    
    pool_dfs[[as.character(pool)]] <- wide_df
  }
  
  # Combine all pools
  final_df <- do.call(rbind, pool_dfs)
  
  # Rename Color column to Index for consistency
  names(final_df)[names(final_df) == "Color"] <- "Index"
  
  # Remove X column if it exists
  if("X." %in% colnames(final_df)) {
    final_df <- final_df %>% select(-X.)
  }
  
  # Add row.names = FALSE to prevent automatic row numbering
  rownames(final_df) <- NULL
  
  return(final_df)
}

# Function to convert percentages into the specified wide format for output
convert_to_formatted_output_wide_before_adding_pool_as_column <- function(solution_color_percentages) {
  possible_characters <- c("R", "G", "B", "-")
  
  # Extract pool numbers from column names
  pools <- sort(as.numeric(unique(str_extract(colnames(solution_color_percentages), "^pool(\\d+)") %>% 
                                    str_replace("pool", ""))))
  
  # Create an empty list to store dataframes for each pool
  pool_dfs <- list()
  
  for(pool in pools) {
    # Filter columns for current pool
    pool_cols <- grep(paste0("^pool", pool, "_"), colnames(solution_color_percentages), value = TRUE)
    pool_data <- solution_color_percentages[, pool_cols]
    
    # Create long format for current pool
    long_df <- pool_data %>%
      pivot_longer(cols = everything(), names_to = "column", values_to = "percentage") %>%
      mutate(
        Color = str_extract(column, "color\\d"),
        Position = as.numeric(str_extract(column, "position(\\d+)") %>% str_replace("position", "")),
        Character = str_extract(column, "[RGB-]$"),
        Percentage = paste0(round(percentage, 1), "%")
      ) %>%
      select(Color, Position, Character, Percentage)
    
    # Convert to wide format
    wide_df <- long_df %>%
      filter(Character %in% possible_characters) %>%
      pivot_wider(names_from = Character, values_from = Percentage, values_fill = list(Percentage = "0%")) %>%
      arrange(Color, Position) %>%
      select(Color, Position, all_of(possible_characters))
    
    # Add pool header
    pool_header <- data.frame(
      Color = paste("=== Pool", pool, "==="),
      Position = NA,
      R = NA, G = NA, B = NA, `-` = NA
    )
    
    # Combine header with data
    pool_dfs[[as.character(pool)]] <- bind_rows(pool_header, wide_df)
    
    # Add blank row if this isn't the last pool
    if(pool != max(pools)) {
      blank_row <- data.frame(
        Color = "",
        Position = NA,
        R = NA, G = NA, B = NA, `-` = NA
      )
      pool_dfs[[as.character(pool)]] <- bind_rows(pool_dfs[[as.character(pool)]], blank_row)
    }
  }
  
  # Combine all pools
  final_df <- do.call(rbind, pool_dfs)
  
  # Rename Color column to Column for consistency
  names(final_df)[names(final_df) == "Color"] <- "Column"
  
  # Remove X column if it exists
  if("X." %in% colnames(final_df)) {
    final_df <- final_df %>% select(-X.)
  }
  
  # Add row.names = FALSE to prevent automatic row numbering
  rownames(final_df) <- NULL
  
  return(final_df)
}

# Function to convert percentages into the specified wide format for output
convert_to_formatted_output_wide_before_split_by_pool <- function(solution_color_percentages) {
  possible_characters <- c("R", "G", "B", "-")
  
  # Reshape and format the dataframe for output
  long_df <- solution_color_percentages %>%
    pivot_longer(cols = everything(), names_to = "column", values_to = "percentage") %>%
    separate(column, into = c("Color", "Position_Char"), sep = "_position") %>%
    mutate(
      Position = as.numeric(str_extract(Position_Char, "\\d+")),
      Character = str_remove(Position_Char, "\\d+_"),
      Percentage = paste0(round(percentage, 1), "%")
    ) %>%
    select(Color, Position, Character, Percentage)
  
  # Remove the "_chars" suffix from the Column values
  long_df <- long_df %>%
    mutate(Column = str_remove(Color, "_chars")) %>%
    select(-Color)
  
  # Convert to wide format with characters as column headings
  wide_df <- long_df %>%
    filter(Character %in% possible_characters) %>%
    pivot_wider(names_from = Character, values_from = Percentage, values_fill = list(Percentage = "0%")) %>%
    arrange(Column, Position) %>%
    select(Column, Position, all_of(possible_characters))
  
  return(wide_df)
}

# Function to write both dataframes to a file
write_to_file <- function(solution, formatted_color_percentages, outputFile) {
  if (!is.null(outputFile)) {
    # Write the solution dataframe
    write.table(solution, outputFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    
    # Format and write the solution_color_percentages dataframe
    outputFile_percentages <- paste0(unlist(strsplit(outputFile, "\\."))[1], "_percentages.txt")
    # header_text <- "\nColumn\tPosition\tR\tG\tB\t-\tPercentage\tColors"
    # write(header_text, file = outputFile_percentages)
    write.table(formatted_color_percentages, outputFile_percentages, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
    
    cat("\nThe proposed sequencing design has been exported into", outputFile)
    cat("\nThe the respective percentages for color balancing have been exported into", outputFile_percentages)
  }
}
