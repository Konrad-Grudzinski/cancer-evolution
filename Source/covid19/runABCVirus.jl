using SubClonalSelection
using Random

input = "m_Netherlands_2020-3_671-filtered0.01"

Random.seed!(123)
out = fitABCmodels(
  "D:/ncbi_dataset/ncbi_dataset/data/vafs/$input.vaf",
  input * "Nmax500",
  read_depth = 600,        # that's basically in the filename, can read the VCF again
  minreads = 5,           # in the VCF, but can also set new
  minvaf = 0.01,           # should do something here, like at least n/500, to filter VAFs of 1/n. OVERRIDES minreads and fmin!!
  fmin = 0.01,            # Can use this instead of minvaf
  fmax = 1.0,            # Maybe set to 1? Could be that mutation becomes dominant after some time
  maxiterations = 10^5, # probably 1000 is enough
  maxclones = 2,          # leave this, want to check all possibilities
  nparticles = 500,       # we are dealing with smaller samples in general
  Nmax = 500,            # Maximum population size used to fit data, adapt to read depth and Nmaxinf 
  resultsdirectory = "D:/ncbi_dataset/ncbi_dataset/data/ABC_results",
  progress = true,        # launch terminal to show progress?
  verbose = true,
  save = true,
  Nmaxinf = 10^6,         # Maybe adapt to actual samples size, together with max iterations? Probably better to leave it, realistic for most smaller countries
  ploidy = 1,             # Virus is monoploid
  d = 0.0,                # unknown, default value
  b = log(3),             # According to https://pubmed.ncbi.nlm.nih.gov/32498136/, R0 is 1.9-6.7, median of 3.38 (=> basic reproduction number). 
  savepopulations = false,
  convergence = 0.005,
  mincellularity = 1.0
  );

println(out)