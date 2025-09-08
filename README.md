# checkMyIndex

Search for a set of compatible indices for your sequencing experiment according to:

* the number of samples
* the desired multiplexing rate (i.e. number of samples per pool/lane)
* the constraint on the indices (none, use each one or each combination only once)
* the inclusion of any indices that are req to appear in the solution

## Modifications in this fork (made by Cornell Genomics Facility)

* Added XLEAP-SBS chemistry (used by the Illumina NovaSeq X Series and NextSeq 1000/2000)
* Added color balancing tab to show the respective percentages for each of the red, green, and blue colors at each position
* Added ability to rank solutions and return the best one from a given number of trials
* Modifications to enable selection of required indices (within a larger set of indices) and find the best solution that includes those required indices

## Input index file(s)

The list of available indices are supplied in one, or two, tab delimited text files without headers. This app requires each row in the file(s) to consist of an ID, a single sequence (in the case of a single index file), or two sequences (in the case of a dual index file), a weight, and a flag indicating whether the index is required to be in the final solution (1 = required; 0 = not required). If the weight is not supplied it will default to "1"; if the required indicator is not supplied it will default to "0". Possible file formats and column contents are given below:
* Two column, single index file, where each row has the format: <ID> <Sequence>
* Three column, single index file, where each row has the format: <ID> <Sequence> <Weight>
* Four column, single index file, where each row has the format: <ID> <Sequence> <Weight> <Required indicator>
* Three column, dual index file, where each row has the format: <ID> <Sequence> <Sequence>
* Four column, dual index file, where each row has the format: <ID> <Sequence> <Sequence> <Weight>
* Five column, dual index file, where each row has the format: <ID> <Sequence> <Sequence> <Weight> <Required indicator>

Any other formats or column orders are not valid and will cause the app to error out or generate incorrect results. Examples of both a four column (for use a Unique Dual-Indexing (UDI) context) and a two-column, single index, file are available in this repository (testCheckMyIndex-i7-i5.txt and inputIndexesExample.txt respectively) to test the application. 

## Shiny application

### Private Cornell version (you must be behind Cornell's firewall to access this version)

Click [Cornell version](http://cbsugg02.biohpc.cornell.edu:8025/CMI2/) to use the Cornell shiny interface of *checkMyIndex*.

### Public website for the original Pasteur Institute version

Click [original version](https://checkmyindex.pasteur.fr/) to use the original shiny interface of *checkMyIndex*.

### Locally

If both *shiny* and *shinyjs* R packages are already installed, one can use the modified application locally running the two following lines in R:

`library(shiny)`

`runGitHub("Cornell-Genomics-Facility/checkMyIndex", launch.browser=TRUE)`

## Rscript command line

*checkMyIndex* can be executed calling `checkMyIndex.r` with `Rscript`. As `checkMyIndex.r` sources `global.r`, both files must be placed in the same directory. Here are 4 examples using `testCheckMyIndex-i7-i5.txt`:

* List of 27 indices for 27 samples distributed on 3 lanes with the two-channel XLEAP-SBS Illumina chemistry:

`Rscript checkMyIndex.r --inputFile7=testCheckMyIndex-i7-i5 -C X -n 27 -m 9`

* List of 27 indices for 27 samples distributed on 3 lanes with the four-channel Illumina chemistry using each lane combination only once:

`Rscript checkMyIndex.r --inputFile7=testCheckMyIndex-i7-i5 -C 4 -n 27 -m 9 -u lane`

* List of 27 indices for 27 samples distributed on 3 lanes with the two-channel original SBS Illumina chemistry using each index only once:

`Rscript checkMyIndex.r --inputFile7=testCheckMyIndex-i7-i5 -C 2 -n 27 -m 9 -u index`

* List of 27 dual-indices for 27 samples distributed on 3 lanes with the two-channel original SBS Illumina chemistry:

`Rscript checkMyIndex.r --inputFile7=testCheckMyIndex-i7-i5 -C 2 -n 27 -m 9`

The help page of the script can be displayed running the following command: 

`Rscript checkMyIndex.r --help`

## Requirements

Here is the list of the R packages needed to run *checkMyIndex*:

* *shiny* and *shinyjs* to run the shiny application locally

* *optparse* to interpret the input parameters when using the Rscript command line

* *parallel* to speed up the calculations

One can install each of these packages running `install.packages(packageName)` in R.

## Illumina chemistry

Illumina has developed four types of chemistry: the four-channel for the HiSeq and MiSeq sequencing devices, the two-channel original SBS and two-channel XLEAP-SBS for the NovaSeq, NextSeq and MiniSeq devices and the one-channel for the iSeq 100 device. With the four-channel chemistry, a red laser detects A/C bases and a green laser detects G/T bases and the indices are compatible if there is at least one red light and one green light at each position. With the two-channel original SBS chemistry, G bases have no color, A bases are orange, C bases are red and T bases are green and indices are compatible if there is at least one color at each position. With the two-channel XLEAP-SBS chemistry, G bases have no color, A bases are blue, C bases are Blue+Green (Cyan) and T bases are green and indices are compatible if there is at least one color at each position. Note that indices starting with GG are not compatible with the two-channel chemistry. With the one-channel chemistry, compatibility cannot be defined with colors and indices are compatible if there is at least one A or C or T base at each position. Please refer to the Illumina documentation for more detailed information on the different chemistries.

## Updates

* Institut Pasteur version 1.0.1: Update of the user interface with *shinyjs*
* Institut Pasteur version 1.0.2: Deal with the new Illumina unique dual-indices (UDI)
* Genomics Facility at Cornell version 1.2.0: Modifications to include XLEAP-SBS chemistry
* Genomics Facility at Cornell version 1.2.1: Bug fix in check for index compatibility for XLEAP-SBS chemistry
* Genomics Facility at Cornell version 1.2.2: Bug fix in check for index compatibility when using multiple pools
* Genomics Facility at Cornell version 1.2.3: Bug fix in server code to use renderDT
* Genomics Facility at Cornell version 1.3.0: Modifications to find the best solution from a given number of trials
* Genomics Facility at Cornell version 1.4.0: Modifications to enable selection of required indices and find best solution that includes those indices
* Genomics Facility at Cornell version 1.4.1: Bug fix to ignore weights when using separate i7 and i5 index files unless in UMI mode
* Genomics Facility at Cornell version 1.4.2: Modifications to add version number to app title text
* Genomics Facility at Cornell version 1.4.3: Automatically stop processing when browser window is closed
* Genomics Facility at Cornell version 1.4.4: Add a warning message to proposed flowcell design tab if any score in a solution is less than 3

## About checkMyIndex

The original tool was developed at the Biomics pole of the Institut Pasteur by Hugo Varet (<hugo.varet@pasteur.fr>) and an [Application Note](https://doi.org/10.1093/bioinformatics/bty706) describing it has been published in 2018 in Bioinformatics. Modifications to include Illumina XLEAP-SBS chemistry have been made by the Genomics Facility at Cornell. Please note that *checkMyIndex* is provided without any guarantees as to its accuracy.
