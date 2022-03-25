using SubClonalSelection
using Random

Random.seed!(123)
out = fitABCmodels("C:/Users/Konrad Grudzinski/OneDrive - University of Glasgow/Computing/4th Year/Individual Project/Source/covid19/SubClonalSelection/example/oneclone.txt",
	"oneclone",
    read_depth = 300,
    nparticles = 100,
    maxiterations = 10^3,
    maxclones = 1, #only identify 0 or 1 subclones
    save = false,
	Nmax = 100,
    adaptpriors = true,
    verbose = false,
    Nmaxinf = 10^6,
    fmin = 0.01)