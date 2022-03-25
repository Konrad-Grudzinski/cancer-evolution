# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 15:00:22 2022

@author: Konrad Grudzinski
"""

from Simulation import run1simulation, exptime, sampledhist
from Classes import InputParameters, Simulation
import numpy as np

def simulate(nclones = 1, ploidy = 2, read_depth = 100.0, detectionlimit = None, clonalmutations = 100.0, mu = 10.0, d = 0.0, b = np.log(2), rho = 0.0, Nmax = 10**4, s = None, tevent = None, cellularity = 1.0, fixedmu = False, timefunction = exptime, maxclonefreq = 200):

    #get simulation data
    simresult, IP, diff = simulate2(nclones, ploidy, read_depth, detectionlimit, clonalmutations, mu, d, b, rho, Nmax, s, tevent, cellularity, fixedmu, timefunction, maxclonefreq)

    #get sampled VAFs
    if IP.rho > 0.0:
        sampleddata = sampledhist(simresult.trueVAF, IP.Nmax, IP.rho, detectionlimit = IP.detectionlimit, ploidy = IP.ploidy, read_depth = IP.read_depth, cellularity = IP.cellularity)
    else:
        sampleddata = sampledhist(simresult.trueVAF, IP.Nmax, None, detectionlimit = IP.detectionlimit, ploidy = IP.ploidy, read_depth = IP.read_depth, cellularity = IP.cellularity)

    return Simulation(IP, simresult, sampleddata), diff



def simulate2(nclones = 1, ploidy = 2, read_depth = 100.0, detectionlimit = None, clonalmutations = 100.0, mu = 10.0, d = 0.0, b = np.log(2), rho = 0.0, Nmax = 10**4, s = None, tevent = None, cellularity = 1.0, fixedmu = False, timefunction = exptime, maxclonefreq = 200):
    if detectionlimit is None:
        detectionlimit = 5 / read_depth
    if s is None:
        s = np.repeat(1.0, nclones)
    if tevent is None:
        tevent = np.arange(1, 100, 0.5)[:nclones]

    IP = InputParameters(
        nclones,
        Nmax,
        detectionlimit,
        ploidy,
        read_depth,
        clonalmutations,
        s,
        (mu/2) * ploidy,
        b,
        d,
        tevent,
        rho,
        cellularity,
        fixedmu,
        timefunction,
        maxclonefreq
    )

    #get simulation data
    return run1simulation(IP, 0.0, 1.0)