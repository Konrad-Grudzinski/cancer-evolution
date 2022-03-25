import datetime
import numpy as np
from scipy.stats import poisson
from Classes import cancercell, RawOutput, SimResult, SampledData

def newmutations(cancercell, mu, mutID):
    cancercell.mutations.append(mutID)
    mutID += 1
    return cancercell, mutID


def newmutationsinit(cancercell, mu, mutID):
    cancercell.mutations = []
    return cancercell, mutID


def initializesim(clonalmutations):
    #initialize time to zero
    t = 0.0
    tvec = [t]

    #population starts with one cell
    N = 1
    Nvec = [N]

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population, lowest fitness
    cells = [cancercell([], 0)]

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutID = 1
    cells[0], mutID = newmutationsinit(cells[0], clonalmutations, mutID)

    return t, tvec, N, Nvec, cells, mutID


# time step functions
exptime = lambda: -np.log(np.random.random())
meantime = lambda: 1


def tumourgrow_birthdeath(b, d, Nmax, mu, numclones = 1, clonalmutations = None, s = [0], tevent=[0], maxclonefreq = 200, timefunction = meantime):

    if clonalmutations is None:
        clonalmutations = mu

    #set array of birthrates
    birthrates = [b]
    deathrates = [d]
    times = np.append(tevent, 0)
    #time is defined in terms of population doublings
    timesN = np.round( np.append( np.exp( np.log(2) * times[:-1] ), 0.0 ) )

    #depending on number of clones add birthrates to model, fitness is randomly distributed between death and birth rates
    for i in range(numclones):
        deathrates.append( np.random.random() * deathrates[0] )
        birthrates.append((1 + s[i]) * (birthrates[0] - deathrates[0]) + deathrates[i + 1])

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units
    Rmax = np.array(b) + np.array(d)

    #initialize arrays and parameters
    t, tvec, N, Nvec, cells, mutID = initializesim(clonalmutations)
    muts = [mutID]

    #we only want to introduce mutant once so have variable that keeps track of how many mutants have been introduced, keep track of which type of which clone aquires new clone
    fitmutant = 0
    clonetype = []
    clonetime = []
    subclonemutations = []
    cloneN = []
    Ndivisions = []
    aveDivisions = []

    clonefreq = np.zeros(numclones + 1)
    clonefreq[0] = 1

    executed = False
    changemutrate = np.array([True for i in range(numclones+1)])
    a = datetime.datetime.now()
    while N < Nmax:

        #pick a random cell
        randcell = np.random.randint(0, N)
        r = np.random.uniform(0, Rmax)
        Nt = N

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < birthrates[cells[randcell].fitness]:

            #population increases by one
            N += 1
            #copy cell and mutations for cell that reproduces
            cells.append(cells[randcell].copy())
            #add new mutations to both new cells
            if mu > 0.0:
              cells[randcell], mutID = newmutations(cells[randcell],mu,mutID)
              cells[-1], mutID = newmutations(cells[-1],mu,mutID)

            muts.append(mutID)
            clonefreq[cells[randcell].fitness] += 1
            Nvec.append(N)
            dt = 1/(Rmax * Nt) * timefunction()
            t += dt
            tvec.append(t)

            #if population time is tevent, cell is mutated into fitter cell
            if N >= timesN[fitmutant]:
                if fitmutant != numclones:
                    #one mutant turns into another "type" so decreases in frequency

                    clonefreq[cells[randcell].fitness] -= 1
                    #keep track of how many clones
                    fitmutant += 1
                    clonetype.append(cells[randcell].fitness)
                    #change one mutant to fitter mutant
                    cells[randcell].fitness = fitmutant
                    #new type increases in frequency
                    clonefreq[cells[randcell].fitness] += 1

                    #change Rmax given that we now have a new fitter mutant
                    Rmax = np.max(birthrates[:fitmutant+1]) + np.max(deathrates[:fitmutant+1])

                    clonetime.append(t)
                    subclonemutations.append(cells[randcell].mutations.copy())
                    cloneN.append(N)
                    Ndivisions.append(len(cells[randcell].mutations))
                    divs = list(map(lambda x: len(x.mutations), cells))
                    aveDivisions.append(np.mean(divs))

        if (birthrates[cells[randcell].fitness] + deathrates[cells[randcell].fitness]) <= r:
            Nvec.append(N)
            dt =  1/(Rmax * Nt) * timefunction()
            t += dt
            tvec.append(t)
        elif birthrates[cells[randcell].fitness] <= r: #death event if b<r<b+d
            #population decreases by 1
            N -= 1
            #frequency of cell type decreases
            clonefreq[cells[randcell].fitness] -= 1
            #remove deleted cell
            cells.pop(randcell)
            Nvec.append(N)
            dt = 1/(Rmax * Nt) * timefunction()
            t += dt
            tvec.append(dt)

        #every cell dies reinitialize simulation
        if (N == 0):
            t, tvec, N, Nvec, cells, mutID = initializesim(clonalmutations)
            muts = [mutID]

        if (executed == False) and np.all((clonefreq > maxclonefreq) == changemutrate):
            #if population of all clones is sufficiently large no new mutations
            #are acquired, can use this approximation as only mutations above 1%
            # frequency can be reliably detected
            mu = 0
            executed = True
    diff = (datetime.datetime.now() - a).total_seconds()
    return RawOutput(Nvec, tvec, muts, cells, birthrates, deathrates, np.array(clonetype), clonetime, subclonemutations, cloneN, Ndivisions, aveDivisions), diff


def cellsconvert(cells):
    #convert from array of cell types to one array with mutations and one array with cell fitness
    l = len(cells)
    mutations = np.hstack([cells[i].mutations for i in range(l)])

    fitness = np.zeros(l)
    for i in range(l):
        fitness[i] = cells[i].fitness

    return mutations, fitness


def allelefreq(mutations):
    #create dictionary that maps mutation ID to allele frequency
    muts = np.array(mutations)[np.argsort(mutations)]
    idx = np.searchsorted(muts,muts)
    idx[idx==len(muts)] = 0
    mask = muts[idx]==muts
    out = np.bincount(idx[mask])
    c = np.unique(muts)
    return dict(zip(c, out))


def getresults(tevent, s, b, d, mu, Nmax, ploidy = 2, clonalmutations = 100, nc = 0, timefunction = exptime, maxclonefreq = 200):
    #Nvec,tvec,mvec,cells,br,dr,ct,clonetime
    sresult, diff = tumourgrow_birthdeath(b, d, Nmax, mu, numclones = nc, s = s, tevent = tevent, clonalmutations = 0, timefunction = timefunction, maxclonefreq = maxclonefreq)
    M, fitness = cellsconvert(sresult.cells)
    return M, fitness, sresult.tvec[-1], sresult.clonetime, sresult.subclonemutations, sresult.birthrates, sresult.deathrates, sresult.cloneN, sresult.clonetype, sresult.Ndivisions, sresult.cells, sresult.aveDivisions, diff


def findallin(a, b):
    indices = []
    barray = np.array(b)
    am = np.max(a)
    g = np.unique(sorted(b))
    g = barray[barray<=am]
    for v in g:
        inds = np.where(a==v)[0]
        if inds.shape != 0:
            indices.extend(inds)
    return np.unique(sorted(indices))


def allelefreqexpand(AFDict, mu, subclonemutations, fixedmu = False):
    #expand allele frequency given mutation rate and calculate number of mutations in the subclones
    if not fixedmu:
        cmuts = np.zeros(len(subclonemutations))
        mutids, mutfreqs = list(zip(*AFDict.items()))
        mutations = poisson.rvs(mu=mu, size=len(mutfreqs))
        AFnew = np.zeros(sum(mutations))

        for i in range(len(cmuts)):
            idx = findallin(mutids, subclonemutations[i])
            cmuts[i] = sum(mutations[idx])

        j = 0
        for f in range(len(mutfreqs)):
            AFnew[j: j + mutations[f]] = np.repeat(mutfreqs[f], mutations[f])
            j = j + mutations[f]

    else:
        mutids, mutfreqs = list(zip(*AFDict.items()))
        muint = round(mu)
        mutations = np.repeat(muint, len(mutfreqs))
        AFnew = np.zeros(sum(mutations))
        cmuts = np.zeros(len(subclonemutations))

        for i in range(len(cmuts)):
            idx = findallin(mutids, subclonemutations[i])
            cmuts[i] = sum(mutations[idx])

        j = 0
        for f in range(len(mutfreqs)):
            AFnew[j: j + mutations[f]] = np.repeat(mutfreqs[f], mutations[f])
            j = j + mutations[f]

    return AFnew, cmuts


def calculateclonefreq(pctfit, cmuts, clonetype):

  cf = pctfit.copy()
  parent = clonetype
  for i in range(len(clonetype),0,-1):
    if parent[i] == 0:
        continue
    cf[parent[i]] = cf[parent[i]] + cf[i]
    cmuts[i] = cmuts[i] - cmuts[parent[i]]

  return cf, cmuts


def run1simulation(IP, minclonesize, maxclonesize):

    M, fitness, tend, clonetime, subclonemutations, br, dr, cloneN, clonetype, Ndivisions, cells, aveDivisions, diff = getresults(IP.tevent, IP.selection, IP.b, IP.d, IP.mu, IP.Nmax, ploidy = IP.ploidy, clonalmutations = IP.clonalmutations, nc = IP.numclones, timefunction = IP.timefunction, maxclonefreq = IP.maxclonefreq)

    if len(clonetime) != IP.numclones:
        IP.numclones = len(clonetime)
        IP.tevent = IP.tevent[:IP.numclones]
        IP.selection = IP.selection[:IP.numclones]
        br = br[:IP.numclones]
        dr = dr[:IP.numclones]
        clonetype = clonetype[:IP.numclones]

    AF = allelefreq(M)#, IP.Nmax)
    AF, cmuts = allelefreqexpand(AF, IP.mu, subclonemutations, fixedmu = IP.fixedmu)
    AF = np.hstack([np.repeat(IP.Nmax, IP.clonalmutations), AF])

    pctfit = []
    for i in range(IP.numclones):
        pctfit.append(sum(fitness==(i+1))/IP.Nmax)

    #remove clones that have frequency < detectionlimit
    if len(pctfit) > 1:
      clonefreq, cmuts = calculateclonefreq(pctfit, cmuts, clonetype - 1)
      if sum(clonefreq > 1) > 0:
          raise Exception(f"There is a clone with frequency greater than 1, this should be impossible ({clonefreq}), {clonetype}, {pctfit}")
    else:
      clonefreq = pctfit


    clonefreq = np.array(clonefreq)
    detectableclones = (clonefreq > minclonesize) & (clonefreq < maxclonesize)

    clonefreq = clonefreq[detectableclones]

    if sum(detectableclones) != IP.numclones:
        IP.numclones = sum(detectableclones)
        IP.tevent = IP.tevent[detectableclones]
        IP.selection = IP.selection[detectableclones]
        clonetype = clonetype[detectableclones]
        detectableclones = np.hstack([True, detectableclones])
        detectableclones = detectableclones[:len(br)]
        br = np.array(br)[detectableclones]
        dr = np.array(dr)[detectableclones]

    return SimResult(clonefreq, pctfit, clonetime, cmuts, br, dr, tend, AF, cloneN, clonetype - 1, Ndivisions, cells, aveDivisions), IP, diff



def betabinom(p, n, rho):
    mu = p * n
    shape1 = ((mu / n) * ((1 / rho) - 1))
    shape2 = abs(n * shape1/mu - shape1)
    return np.random.binomial(n, np.random.beta(shape1, shape2))


def sampledhist(AF, cellnum, rho = None, detectionlimit = 0.1, ploidy = 2.0, read_depth = 100.0, cellularity = 1.0):

    AF = (AF / ploidy) * cellularity
    AF_array = np.array(AF)
    AF = AF_array[AF_array > detectionlimit * cellnum]
    #samp_percent = read_depth / cellnum
    #depth = np.random.binomial(cellnum,samp_percent, size=len(AF))

    depth = poisson.rvs(mu=read_depth, size=len(AF))
    depths = np.repeat(depth, len(AF))
    if rho is not None:
        samp_alleles = list( map(lambda x: betabinom(x[0], x[1], rho), zip(AF / cellnum, depths)) )
    else:
        samp_alleles = list( map(lambda x: np.random.binomial(x[0], x[1]), zip(depths, AF / cellnum)) )
    VAF = np.array(samp_alleles) / depth

    # data for histogram
    edges = np.arange(0.005, 1.006, 0.01)
    hist, _ = np.histogram(VAF, bins=edges)
    DFhist = np.vstack([edges[:-1], hist]).T

    return SampledData(DFhist, VAF, samp_alleles, depth)


def getCDF(VAF, step_size, fmin = 0.05, fmax = 0.7):
  #fast way to calculate CDF
  edges = np.linspace(fmin, fmax, (fmax-fmin)/step_size)
  hist, _ = np.histogram(VAF, bins=edges)
  out = np.cumsum(hist)
  out = out - out[0]
  return out
