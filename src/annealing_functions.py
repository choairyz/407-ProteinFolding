# Protein folding final project for PHY407, by Victor Choairy - 1011464149
from random import random
import numpy as np
import matplotlib.pyplot as plt
from protein import Protein
from random import randrange, choice, random, seed


seed_num = 1
seed(seed_num)


def acceptance(E_now, E_test, T):
    """Accept or reject the change from E_test to E_now using the Marvok Chain Monte Carlo algorithm, at temperature T"""
    if E_now <= E_test: # always accept energy lowerings
        return True
    else:
        P = np.exp(-(E_now - E_test)/(T))
        
    num = random()
    if num <= P:
        return True
    else:
        return False
    

def cooledT(T0, tau, t):
    """Returns temperature T at time t given initial temperature T0 and decay constant tau, for simulated annealing."""
    return T0 * np.exp(-t/tau)


def discreteAnnealing(initialProtein: Protein, tau, Tstart, Tfinal):
    """TODO"""
    stateList = [] # all guesses (states) during annealing
    stateList = [Protein(initialProtein.coords, initialProtein.width, initialProtein.height)]

    # track temperature and time for plotting cooling schedule
    tempList = []
    tList = []
    # track energy for plotting as well
    EList = []

    T = Tstart # temperature
    t = 0.

    tempList.append(Tstart) 
    tList.append(0.)
    
    E = initialProtein.calc_energy() # energy
    EList.append(E)
    
    while T > Tfinal:
        t += 1.0
        tList.append(t)
        T = cooledT(Tstart, tau, t) # update T
        tempList.append(T)
        while True:
            i = randrange(len(initialProtein.coords)) # random monomer
            dx, dy = choice([-1, 1]), choice([-1, 1]) # random diagonal step
            x, y = initialProtein.coords[i]
            new = (x + dx, y + dy)

            inside = 0 <= new[0] < initialProtein.width and 0 <= new[1] < initialProtein.height
            empty = new not in initialProtein.coords
            # must stay one step from the monomers it's bonded to (i-1 and i+1, if they exist)
            bonded_ok = True
            for j in (i - 1, i + 1):
                if 0 <= j < len(initialProtein.coords):
                    nx, ny = initialProtein.coords[j]
                    if abs(new[0] - nx) + abs(new[1] - ny) != 1:
                        bonded_ok = False

            if inside and empty and bonded_ok:
                break

        old = initialProtein.coords[i]
        initialProtein.coords[i] = new
        E_new = initialProtein.calc_energy()
        if acceptance(E_new, E, T):
            E = E_new 
        else:
            initialProtein.coords[i] = old 

        EList.append(E)
        stateList.append(Protein(initialProtein.coords, initialProtein.width, initialProtein.height))

    # last entry of stateList is the state corresponding to global min
    return stateList, tList, tempList, EList







