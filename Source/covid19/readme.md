# Readme

The code here can process sequenced genomes in FASTA format, retrieve meta data about each sample from the NCBI data base, perform ABC analysis using `SubClonalSelection.jl` on the data, and display a number of plots for further analysis. The order in which they notebooks should be run is described below in *Build steps*. `covid_utilities.py` only contains some variables and constants which are shared by several notebooks.

## Build instructions


### Requirements

* Python 3.9.7
* Python packages listed in `requirements.txt`
* WSL
* Julia programming language
* `SubClonalSelection.jl` Julia package
* `bcftools`
* `samtools`
* `mafft`
* `hisat2`
* `ncbi_cmdtools`
* Tested on Windows 10

### Build steps

1. Install all the above software.
2. Run `pip install -r requirements.txt` to install all Python dependencies.
3. To download the `SubClonalSelection.jl` package, run `] add https://github.com/marcjwilliams1/SubClonalSelection.jl` in a Julia session


### Test steps

If all above steps work, then everything should be fine.
