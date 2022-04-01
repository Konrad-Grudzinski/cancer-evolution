# User manual 

1. Install all the software as listed in the readme file.
2. Download the sequence data in FASTA format using ncbi_cmdtools (at the time of download in early March 2022, >20GB)
3. Set the paths at the top of `covid_utilities.py` to wherever the source code lies and the sequences are downloaded to
4. Split the downloaded sequences into individual files by running `split_files_by_sample.ipynb`
5. Download the meta data with `download_metadata.ipynb`
6. Get EUClusters.json from [Covariants.org](https://github.com/hodcroftlab/covariants/blob/master/cluster_tables/EUClusters_data.json)
7. Group samples by country and collection date using `sort_data.ipynb`
8. Build an index on the reference sequence NC_045512.2 using hisat2-build (included with hisat2). This can be done with `hisat2-build NC_045512.2.fa NC_045512.2.fa`
9. To process all the samples, run `createVCF.ipynb`
10. Create a Julia sysimage with `SubClonalSelection.jl` included, for example using the `PackageCompiler.jl` package. This step is technically *optional*, but the ABC analysis will have to be run sequentially (i.e. with a single thread instead of 6). For more information see [here](https://docs.julialang.org/en/v1/devdocs/sysimg/)
11.  The inference process on all grouped samples is performed by running `run_abc_analysis.ipynb`. Reduce the number of threads if required. This may take 30 hours or more 
12. To view plots of frequencies and probabilities, run `plot_probabilities_and_frequencies.ipynb`. To view predictions, run `predictive_simulation.ipynb`.
