MIST
===

___Mass spectrometry interaction STatistics (MiST)___

[Check out the paper](http://www.ncbi.nlm.nih.gov/pubmed/25754993) for details about how to use MIST.

## Necessary Resources

__Hardware__

Workstation running any current OS, Unix environment recommended

__Software__

- [R](http://www.r-project.org)
- R packages: getopt, optparse, reshape2, pheatmap, RcolorBrewer, ggplot2, MESS, yaml
- [MiST source code](https://github.com/kroganlab/private.mist)
- Git (optional, but highly recommended)

## Installation

- Download the MiST package as a .zip archive from the public GitHub repository by clicking on the “Download ZIP” button on the bottom right, unzip the files and
- move the directory to a permanent location.
    + Alternatively, you can check out the MiST package through Git as follows:
    + `git clone https://github.com/kroganlab/private.mist MiST`
- The MiST pipeline is designed to run from a terminal using R. This requires the user to have executable permissions. To set these permissions in a Unix environment, navigate in the terminal to the MiST directory, hereafter referred to as the `$INSTALL_DIR`, then type: `sudo chmod -R 775 *`

## Config file

Check the sample file in `test/sample.yaml`

## Run MIST

`~/github/kroganlab/private.mist/main.R -c sample.yaml`


