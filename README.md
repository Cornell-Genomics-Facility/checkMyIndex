# checkMyIndex

Search for a set of compatible indexes for your sequencing experiment according to:

* the number of samples
* the desired multiplexing rate (i.e. number of samples per pool/lane)
* the constraint on the indexes (none, use each one or each combination only once)

## Input indexes file

The list of the available indexes must be stored in a text file containing two tab-separated columns (without header): index ids are in the first column and the corresponding sequences in the second. `inputIndexesExample.txt` is an example of such a file and can be used to test *checkMyIndex*. For dual-indexing sequencing experiments the user needs to provide two files: the first one for the indexes 1 (i7) and the second for the indexes 2 (i5). `index24-i7.txt` and `index24-i5.txt` are available to test the research of compatible dual-indexes.

## Shiny application

### Public website

Click [here](https://checkmyindex.pasteur.fr/) to use the shiny interface of *checkMyIndex*.

### Locally

One can use the application locally running the two following lines in R:

`library(shiny)`

`runGitHub("PF2-pasteur-fr/checkMyIndex", launch.browser=TRUE)`

## Rscript command line

*checkMyIndex* can be executed calling `checkMyIndex.r` with `Rscript`. As `checkMyIndex.r` sources `global.r`, both files must be placed in the same directory. Here are 4 examples using either `inputIndexesExample.txt` or `index24-i7.txt` and `index24-i5.txt`:

* List of 9 indexes for 9 samples distributed on 3 lanes with the four-channel Illumina chemistry:

`Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C 4 -n 9 -m 3`

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

* *shiny* to run the shiny application locally

* *optparse* to interpret the input parameters when using the Rscript command line

* *parallel* to speed up the calculations

One can install each of these packages running `install.packages(packageName)` in R.

## Illumina chemistry

Illumina has developed three types of chemistry: the four-channels for the HiSeq and MiSeq sequencing devices, the two-channels for the NovaSeq, NextSeq and MiniSeq devices and the one-channel for the iSeq 100 device. With the four-channel chemistry A/C are red and G/T are green and indexes are compatible if there are at least one red light and one green light at each position. With the wo-channel chemistry G has no color, A is orange, C is red and T is green and indexes are compatible if there is at least one color at each position. Note that indexes starting with GG are not compatible with the two-channel chemistry. With the one-channel chemistry compatibility cannot be defined with colors and indexes are compatible is there is at least one A or C or T at each position. Please refer to the Illumina documentation for more details.

## About checkMyIndex

This tool has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet (<hugo.varet@pasteur.fr>). Please note that *checkMyIndex* is provided without any guarantees as to its accuracy.
