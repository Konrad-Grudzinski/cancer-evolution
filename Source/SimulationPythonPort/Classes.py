from pprint import pformat as pp

############# objects #################

class cancercell:

    def __init__(self, mutations, fitness):
        self.mutations = mutations
        self.fitness = fitness

    def copy(self):
        return cancercell(self.mutations.copy(), self.fitness)


class InputParameters:
    def __init__(self, numclones, Nmax, detectionlimit, ploidy, read_depth,
                 clonalmutations, selection, mu, b, d, tevent, rho,
                 cellularity, fixedmu, timefunction, maxclonefreq):
        self.numclones = numclones
        self.Nmax = Nmax
        self.detectionlimit = detectionlimit
        self.ploidy = ploidy
        self.read_depth = read_depth
        self.clonalmutations = clonalmutations
        self.selection = selection
        self.mu = mu
        self.b = b
        self.d = d
        self.tevent = tevent
        self.rho = rho
        self.cellularity = cellularity
        self.fixedmu = fixedmu
        self.timefunction = timefunction
        self.maxclonefreq = maxclonefreq


class RawOutput:
    def __init__(self, Nvec, tvec, muts, cells, birthrates, deathrates,
                 clonetype, clonetime, subclonemutations, cloneN,
                 Ndivisions, aveDivisions):
        self.Nvec = Nvec
        self.tvec = tvec
        self.muts = muts
        self.cells = cells
        self.birthrates = birthrates
        self.deathrates = deathrates
        self.clonetype = clonetype
        self.clonetime = clonetime
        self.subclonemutations = subclonemutations
        self.cloneN = cloneN
        self.Ndivisions = Ndivisions
        self.aveDivisions = aveDivisions


class SimResult:
    def __init__(self, clonefreq, clonefreqp, clonetime, subclonemutations,
                birthrates, deathrates, tend, trueVAF, cloneN, clonetype,
                Ndivisions, cells, aveDivisions):
        self.clonefreq = clonefreq
        self.clonefreqp = clonefreqp
        self.clonetime = clonetime
        self.subclonemutations = subclonemutations
        self.birthrates = birthrates
        self.deathrates = deathrates
        self.tend = tend
        self.trueVAF = trueVAF
        self.cloneN = cloneN
        self.clonetype = clonetype
        self.Ndivisions = Ndivisions
        self.cells = cells
        self.aveDivisions = aveDivisions

    def __str__(self):
        return "\n".join([
            f"clonefreq: {pp(self.clonefreq)}", f"clonefreqp: {pp(self.clonefreqp)}", f"clonetime: {pp(self.clonetime)}",
            f"subclonemutations: {pp(self.subclonemutations)}", f"birthrates: {pp(self.birthrates)}", f"deathrates: {pp(self.deathrates)}",
            f"tend: {pp(self.tend)}", f"trueVAF: {pp(self.trueVAF)}", f"cloneN: {pp(self.cloneN)}", f"clonetype: {pp(self.clonetype)}",
            f"Ndivisions: {pp(self.Ndivisions)}", f"#cells: {len(self.cells)}",f"aveDivisions: {pp(self.aveDivisions)}"])


class SampledData:
    def __init__(self, DFhist, VAF, samp_alleles, depth):
        self.DFhist = DFhist
        self.VAF = VAF
        self.samp_alleles = samp_alleles
        self.depth = depth


class Simulation:
    def __init__(self, InputParameters, SimResult, SampledData):
        self.InputParameters = InputParameters
        self.SimResult = SimResult
        self.SampledData = SampledData