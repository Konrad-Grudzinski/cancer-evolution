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
2. Download the sequence data in FASTA format using ncbi_cmdtools (at the time of download in early March 2022, >20GB)
3. Set the paths at the top of `covid_utilities.py` to wherever the source code lies and the sequences are downloaded to
4. Split the downloaded sequences into individual files by running `split_files_by_sample.ipynb`
5. Download the meta data with `download_metadata.ipynb`
6. Get EUClusters.json from Covariants.org
7. Group samples by country and collection date using `sort_data.ipynb`
8. Build an index on the reference sequence NC_045512.2 using hisat2-build (included with hisat2). This can be done with `hisat2-build NC_045512.2.fa NC_045512.2.fa`
9. To process all the samples, run `createVCF.ipynb`
10. Create a Julia sysimage with SubClonalSelection.jl included, for example using the `PackageCompiler.jl` package. This step is technically *optional*, but the ABC analysis will have to be run sequentially (i.e. with a single thread instead of 6). For more information see [here](https://docs.julialang.org/en/v1/devdocs/sysimg/)
11.  The inference process on all grouped samples is performed by running `run_abc_analysis.ipynb`. Reduce the number of threads if required. This may take 30 hours or more 
12. To view plots of frequencies and probabilities, run `plot_probabilities_and_frequencies.ipynb`. To view predictions, run `predictive_simulation.ipynb`.


### Test steps

If all above steps work, then everything works.
