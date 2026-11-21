# Protein folding final project for PHY407, by Victor Choairy - 1011464149

from random import random
import numpy as np
import matplotlib.pyplot as plt


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
    

def draw_normal(sigma):
    """ Draw two random numbers from a normal distribution of mean 0, standard dev sigma. From Newman section 10.1.6 """
    theta = 2*np.pi*random()
    z = random()
    r = (-2*sigma**2 * np.log(1 - z))**.5
    return r*np.cos(theta), r*np.sin(theta)


def cooledT(T0, tau, t):
    """Returns temperature T at time t given initial temperature T0 and decay constant tau, for simulated annealing."""
    return T0 * np.exp(-t/tau)


def annealing(f, initialGuess, tau, Tstart, Tfinal):
    """TODO"""
    stateList = [] # all guesses (states) during annealing
    stateList.append(initialGuess)

    # track temperature and time for plotting cooling schedule
    tempList = []
    tList = []
    # track energy for plotting as well
    EList = []

    T = Tstart # temperature
    t = 0.

    tempList.append(Tstart) 
    tList.append(0.)
    
    x, y = initialGuess
    E = f(x, y) # energy
    EList.append(E)
    
    while T > Tfinal:
        t += 1.0
        tList.append(t)
        T = cooledT(Tstart, tau, t) # update T
        tempList.append(T)
        dx, dy = draw_normal(1.0) # pick a change in x and y from a normal distribution of mean 0 and std 1
        x_new, y_new = x + dx, y + dy
        E_new = f(x_new, y_new)
        # decide whether to accept or reject the new state
        if acceptance(E_new, E, T):
            x, y, E = x_new, y_new, E_new
        EList.append(E)
        # append current state to stateList
        stateList.append((x, y))

    # last entry of stateList is the state corresponding to global min
    return stateList, tList, tempList







