library(shiny)
library(shinyjs)

shinyUI(fluidPage(theme = "bootstrap.min.css", shinyjs::useShinyjs(),

                  # Add the CSS styling
                  tags$head(
                    tags$style(
                      "#compatibilityWarning,
                       #compatibilityWarning2,
                       #compatibilityWarning3 {
                          color: red;
                          text-align: center;
                          font-weight: bold;
                          margin-top: 20px;  /* This adds space above the warning */
                      }"
                    )
                  ),
                  
                  titlePanel(title=div(img(src="dna.png", width=50), 
                            strong("Search for a set of compatible indices for your sequencing experiment"), 
                            span(style = "font-size: 80%;", paste0("(version ", cornellCheckMyIndexVersion, ")"))), 
                            windowTitle="checkMyIndex"),
                  
                  sidebarLayout(
                    
                    # parameters
                    sidebarPanel(div(id="allParameters",
                      
                      # input indices
                      # //--- additional text for combined file
                      conditionalPanel(condition="!output.testdataProvided", {
                        fileInput("inputFile", label="Select your tab-delimited file containing the index 1 (i7) ids and sequences, or combined indices (i7 and i5)", accept="text")
                      }),
                      conditionalPanel(condition="!output.testdataProvided", {
                        fileInput("inputFile2", label="Optional index 2 (i5) file for dual-indexing", accept="text")
                      }),
                      conditionalPanel(condition="!output.inputFileProvided & !output.inputFile2Provided", {
                        selectizeInput("testdata", label="Load test indices",
                                       choices=c("None"="none",
                                                 "24 indices 1 (i7)"="simple",
                                                 "27 indices 1 (i7) and 27 indices 2 (i5)"="dual"))
                      }),
                      
                      # parameters
                      selectizeInput("chemistry", label="Illumina chemistry", 
                                     choices=c("Four-channels (HiSeq & MiSeq)" = 4, 
                                               "Two-channels (original SBS)" = 2,
                                               "Two-channels (XLEAP-SBS)" = "X",
                                               "One-channel (iSeq 100)" = 1)),
                      conditionalPanel(condition="output.i7i5sameLength", {uiOutput("i7i5pairing")}),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("nbSamples")}),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("multiplexingRate")}),
                      conditionalPanel(condition="!output.indexUploaded2", {selectInput("unicityConstraint", label="Constraint on the indices (single-indexing only)", 
                                                                                      choices=c("None" = "none", "Use each combination only once" = "lane", "Use each index only once" = "index"),
                                                                                      selected="none")}),
                      conditionalPanel(condition="output.indexUploaded & !output.indexUploaded2", {checkboxInput("completeLane", "Directly look for a solution with the desired multiplexing rate", value=FALSE)}),
                      conditionalPanel(condition="output.indexUploaded & !output.indexUploaded2", {checkboxInput("selectCompIndexes", "Select compatible indices before looking for a solution", value=FALSE)}),
                      selectInput("nbMaxTrials", label="Maximum number of trials to find a solution", choices=10^(1:4)),
                      
                      # go and reset buttons
                      actionButton("go", label="Search for a solution"),
                      actionButton("reset", "Reset parameters")
                      
                    )),
                    
                    # output
                    mainPanel(
                      
                      tabsetPanel(id="mainPanel",
                        # 1st panel: input
                        tabPanel("Input indexes",
                                 value="inputIndexes",
                                 p(textOutput("textIndex")),
                                 DT::DTOutput("inputIndex"),
                                 p(textOutput("textIndex2")),
                                 DT::DTOutput("inputIndex2"),
                                 p(textOutput("textPairingTable")),
                                 DT::DTOutput("pairingTable")),
                        
                        # 2nd panel: results
                        tabPanel("Proposed flowcell design",
                                 value="proposedSolution",
                                 shinyjs::hidden(div(
                                   id="proposedSolution",
                                   p(textOutput("textDescribingSolution"),
                                     textOutput("compatibilityWarning")),
                                   DT::DTOutput("solution"),
                                   br(""),
                                   uiOutput("downloadButton")
                                 ))
                        ),
                        
                        # 3rd panel: color balancing
                        tabPanel("Color balancing",
                                 value="colorBalancing",
                                 shinyjs::hidden(div(
                                   id="colorBalancing",
                                   p(textOutput("textDescribingColorBalancing"),
                                     textOutput("compatibilityWarning2")),
                                   DT::DTOutput("colorPercentagesFormattedOutputWide"),
                                   br(""),
                                   uiOutput("downloadColorBalanceButton")
                                 ))
                        ),
                        
                        # 4th panel: plot results
                        tabPanel("Visualization of the design",
                                 shinyjs::hidden(div(
                                   id="visualization",
                                   p(textOutput("textDescribingHeatmap"),
                                     textOutput("compatibilityWarning3")),
                                   uiOutput("heatmapindex2")
                                 ))
                                 ),
                        
                        # 5th panel: help
                        # //--- add text for combined file and XLEAP chemistry
                        tabPanel("Help",
                                 
                                 h3("What checkMyIndex does"),
                                 p("Searches for a set of compatible indices for your sequencing experiment according to the number of samples and the desired multiplexing rate (i.e. number 
                                   of samples per pool/lane). The app returns the best, color balanced solution it can find (assuming one exists) from the number of trials you select (note: 
                                   better solutions may be found by running for more trials, but this takes longer to run). The app will also allow the selection of required indices (within 
                                   a larger set of indices) and find best the solution that includes those required indices."),

                                 h3("Input index file(s)"),
                                 p("The list of available indices are supplied in one, or two, tab delimited text files without headers. This app requires each row in the file(s) to consist of 
                                   an ID, a single sequence (in the case of a single index file), or two sequences (in the case of a dual index file), a weight, and a flag indicating whether 
                                   the index is required to be in the final solution (1 = required; 0 = not required). If the weight is not supplied it will default to '1'; if the required 
                                   indicator is not supplied it will default to '0'. Possible file formats and column contents are given below:"),
                                 tags$ul(
                                   tags$li("Two column, single index file, where each row has the format: <ID> <Sequence>"),
                                   tags$li("Three column, single index file, where each row has the format: <ID> <Sequence> <Weight>"),
                                   tags$li("Four column, single index file, where each row has the format: <ID> <Sequence> <Weight> <Required indicator>"),
                                   tags$li("Three column, dual index file, where each row has the format: <ID> <Sequence> <Sequence>"),
                                   tags$li("Four column, dual index file, where each row has the format: <ID> <Sequence> <Sequence> <Weight>"),
                                   tags$li("Five column, dual index file, where each row has the format: <ID> <Sequence> <Sequence> <Weight> <Required indicator>")
                                 ),
                                 p("Any other formats or column orders are not valid and will cause the app to error out or generate incorrect results. 
                                   Examples of both a four-column and a two-column file are available in the GitHub repository ", 
                                   a("here", href="https://github.com/Cornell-Genomics-Facility/checkMyIndex/blob/master/www/testCheckMyIndex-i7-i5.txt", target="blank", download="testCheckMyIndex-i7-i5.txt"), 
                                   "and ", a("here", href="https://github.com/Cornell-Genomics-Facility/checkMyIndex/blob/master/www/inputIndexesExample.txt", target="blank", download="inputIndexesExample.txt"), 
                                   "to test the application."),
                                 
                                 h3("Parameters"),
                                 p(strong("Total number of samples"), "in your experiment (can be greater than the number of available indices)."),
                                 p(strong("Multiplexing rate"), "i.e. number of samples per pool/lane (only divisors of the total number of samples are proposed)."),
                                 p(strong("i7 and i5 pairing"), "(only for dual-indexing) is proposed if there are as many i5 as i7 indices to deal with Illumina Unique Dual-Indices (UDI).
                                   Note that the pairing is done using the order of the indices in the input files."),
                                 p(strong("Constraint on the indices"), "(only for single-indexing) to avoid having two samples or two pools/lanes with the same index(es)."),
                                 p(strong("Directly look for a solution with the desired multiplexing rate"), "(only for single-indexing) instead of looking for a partial solution 
                                           with a few samples per pool/lane and then add some of the remaining indices to reach the desired multiplexing rate."),
                                 p(strong("Select compatible indices"), "(only for single-indexing) before looking for a (partial) solution can take some time but then speed up the algorithm."),
                                 p(strong("Maximum number of trials"), "can be increased if a solution is difficult to find with the parameters chosen."),
                                 
                                 h3("How the algorithm works"),
                                 p("There can be many combinations of indices to check according to the number of input indices and the multiplexing rate. Thus, testing for 
                                    the compatibility of all the combinations may be long or even impossible. The trick is to find a partial solution with the desired number 
                                    of pools/lanes but with fewer samples than asked and then to complete each pool/lane with some of the remaining indices to reach the desired 
                                    multiplexing rate. Indeed, adding indices to a combination of compatible indices will give a compatible combination. Briefly, a lower 
                                    number of samples per pool/lane generates a lower number of combinations to test and thus makes the research of a partial solution very fast. 
                                    Adding some indices to complete each pool/lane is fast too and gives the final solution."),
                                 p("Unfortunately, the research of a final solution might become impossible as the astuteness reduces the number of combinations of indices.
                                    In such a case, one can look for a solution using directly the desired multiplexing rate (see parameters), the only risk is to increase 
                                    the computational time."),
                                 
                                 h3("Background on Illumina chemistry and color balancing"),
                                 p(strong("Illumina chemistry"), "can be either four-channels (HiSeq & MiSeq), two-channels (original SBS and XLEAP-SBS) or one-channel (iSeq 100).
                                   With the four-channel chemistry, a red laser detects A/C bases and a green laser detects G/T bases and the indices are compatible if there is at 
                                   least one red light and one green light at each position. With the two-channel chemistry (original SBS), G bases have no color, A bases are 
                                   orange, C bases are red and T bases are green and indices are compatible if there is at least one color at each position. For two-channel XLEAP-SBS chemistry, 
                                   G bases have no color, A bases are blue, C bases are Blue+Green (Cyan) and T bases are green and indices are compatible if there is at least one color at each 
                                   position. Note that indices starting with GG are not compatible with the two-channel chemistry. With the one-channel chemistry, compatibility cannot be 
                                   defined with colors and indices are compatible if there is at least one A or C or T base at each position. Please refer to the Illumina documentation for more 
                                   detailed information on the different chemistries."),

                                 h3("About"),
                                 p("The original application was developed at the Biomics pole of the Institut Pasteur by Hugo Varet and an ", a("Application Note", href="https://doi.org/10.1093/bioinformatics/bty706"), 
                                   "describing it has been published in 2018 in Bioinformatics. Send an e-mail to", a("hugo.varet@pasteur.fr"), "for any suggestion or bug report."),
                                 p("Source code and instructions to run the original application locally are available at the ", a("PF2 - Institut Pasteur GitHub", href="https://github.com/PF2-pasteur-fr/checkMyIndex"), "repository. "),

                                 p("Modifications to include Illumina XLEAP-SBS chemistry have been made by the Genomics Facility at Cornell. These modifcations are available at the ", a("Cornell Genomics Facility GitHub", href="https://github.com/Cornell-Genomics-Facility/checkMyIndex"), "repository. "),

                                 p("Please note that checkMyIndex is provided without any guarantees as to its accuracy."),
                                 
                                 
                                 h3("Version"),
                                 p(paste0("Software based on the Institut Pasteur checkMyIndex version ", checkMyIndexVersion, ".")),
                                 p(paste0("Genomics Facility at Cornell checkMyIndex version ", cornellCheckMyIndexVersion, ", with modifications to include XLEAP-SBS chemistry.")),
                                 div(img(src="logo_c3bi_citech.jpg", width=300), style="text-align: center;"))
                        
                      )
                    )
                  )
))
