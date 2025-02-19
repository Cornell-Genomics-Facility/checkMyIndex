# checkMyIndex

Search for a set of compatible indexes for your sequencing experiment according to:

* the number of samples
* the desired multiplexing rate (i.e. number of samples per pool/lane)
* the constraint on the indexes (none, use each one or each combination only once)

## Modifications in this fork (made by Cornell Genomics Facility)

* Added XLEAP-SBS chemistry (used by the Illumina NovaSeq X Series and NextSeq 1000/2000)
* Added color balancing tab to show the respective percentages for each of the red, green, and blue colors at each position

## Input indexes file

The list of the available indexes must be stored in a text file containing two tab-separated columns (without header): index ids are in the first column and the corresponding sequences in the second. `inputIndexesExample.txt` is an example of such a file and can be used to test *checkMyIndex*. For dual-indexing sequencing experiments the user needs to provide two files: the first one for the indexes 1 (i7) and the second for the indexes 2 (i5). `index24-i7.txt` and `index24-i5.txt` are available to test the research of compatible dual-indexes. Moreover, two files `index96_UDI-i5.txt` and `index96_UDI-i7.txt` are available to test the application in a Unique Dual-Indexing (UDI) context.

## Shiny application

### Public website for the original Pasteur Institute version

Click [here](https://checkmyindex.pasteur.fr/) to use the shiny interface of *checkMyIndex*.

### Locally

If both *shiny* and *shinyjs* R packages are already installed, one can use the modified application locally running the two following lines in R:

`library(shiny)`

`runGitHub("Cornell-Genomics-Facility/checkMyIndex", launch.browser=TRUE)`

## Rscript command line

*checkMyIndex* can be executed calling `checkMyIndex.r` with `Rscript`. As `checkMyIndex.r` sources `global.r`, both files must be placed in the same directory. Here are 4 examples using either `inputIndexesExample.txt` or `index24-i7.txt` and `index24-i5.txt`:

* List of 9 indexes for 9 samples distributed on 3 lanes with the two-channel XLEAP-SBS Illumina chemistry:

`Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C X -n 9 -m 3`

* List of 12 indexes for 12 samples distributed on 4 lanes with the four-channel Illumina chemistry using each lane combination only once:

`Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C 4 -n 12 -m 3 -u lane`

* List of 12 indexes for 12 samples distributed on 4 lanes with the two-channel Illumina chemistry using each index only once:

`Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C 2 -n 12 -m 3 -u index`

* List of 24 dual-indexes for 24 samples distributed on 2 lanes with the two-channel Illumina chemistry:

`Rscript checkMyIndex.r --inputFile7=index24-i7.txt --inputFile5=index24-i5.txt -C 2 -n 24 -m 12`

The help page of the script can be displayed running the following command: 

`Rscript checkMyIndex.r --help`

## Requirements

Here is the list of the R packages needed to run *checkMyIndex*:

* *shiny* and *shinyjs* to run the shiny application locally

* *optparse* to interpret the input parameters when using the Rscript command line

* *parallel* to speed up the calculations

One can install each of these packages running `install.packages(packageName)` in R.

## Illumina chemistry

Illumina has developed four types of chemistry: the four-channel for the HiSeq and MiSeq sequencing devices, the two-channel original SBS and two-channel XLEAP-SBS for the NovaSeq, NextSeq and MiniSeq devices and the one-channel for the iSeq 100 device. With the four-channel chemistry, a red laser detects A/C bases and a green laser detects G/T bases and the indexes are compatible if there is at least one red light and one green light at each position. With the two-channel original SBS chemistry, G bases have no color, A bases are orange, C bases are red and T bases are green and indexes are compatible if there is at least one color at each position. With the two-channel XLEAP-SBS chemistry, G bases have no color, A bases are blue, C bases are Blue+Green (Cyan) and T bases are green and indices are compatible if there is at least one color at each position. Note that indexes starting with GG are not compatible with the two-channel chemistry. With the one-channel chemistry, compatibility cannot be defined with colors and indexes are compatible if there is at least one A or C or T base at each position. Please refer to the Illumina documentation for more detailed information on the different chemistries.

## News

* Version 1.0.1: update of the user interface with *shinyjs*
* Version 1.0.2: deal with the new Illumina unique dual-indexes (UDI)
* Version 1.0.2: forked and modified to add the XLEAP-SBS chemistry

## About checkMyIndex

The original tool was developed at the Biomics pole of the Institut Pasteur by Hugo Varet (<hugo.varet@pasteur.fr>) and an [Application Note](https://doi.org/10.1093/bioinformatics/bty706) describing it has been published in 2018 in Bioinformatics. Modifications to include Illumina XLEAP-SBS chemistry have been made by the Genomics Facility at Cornell. Please note that *checkMyIndex* is provided without any guarantees as to its accuracy.
