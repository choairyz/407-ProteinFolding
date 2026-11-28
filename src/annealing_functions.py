"""
Protein folding final project for PHY407, by Victor Choairy - 1011464149

The Monte Carlo part of the simulation: the Metropolis acceptance rule, three cooling schedules (exponential, linear
and inverse), and discreteAnnealing, which runs the simulated annealing of a lattice protein.

Requirements: NumPy, and protein.py in the same folder.
Outputs: none when run on its own; protein_folding_toolkit.py calls these functions.
"""
import numpy as np
from protein import Protein
from random import randrange, choice, random


def acceptance(E_new, E_old, T):
    """Metropolis rule at temperature T: always accept a move that does not raise the energy (E_new <= E_old),
    and accept a move that raises it with probability exp(-(E_new - E_old)/T)"""
    if E_new <= E_old: # always accept energy lowerings
        return True
    else:
        P = np.exp(-(E_new - E_old)/(T))

    num = random()
    if num <= P:
        return True
    else:
        return False
    

def exponential_cooling(Tstart, tau, t):
    """T = Tstart * exp(-t/tau). Reaches Tfinal after tau*ln(Tstart/Tfinal) steps"""
    return Tstart * np.exp(-t/tau)


def linear_cooling(Tstart, tau, t):
    """T = Tstart * (1 - t/tau), which would reach 0 at step tau. Reaches Tfinal after tau*(1 - Tfinal/Tstart) steps"""
    return Tstart * (1 - t/tau)


def inverse_cooling(Tstart, tau, t):
    """T = Tstart / (1 + t/tau). Reaches Tfinal after tau*(Tstart/Tfinal - 1) steps"""
    return Tstart / (1 + t/tau)


def discreteAnnealing(initialProtein: Protein, tau, Tstart, Tfinal, schedule=exponential_cooling, max_steps=500000):
    """Simulated annealing for a protein chain, starting from initialProtein. The temperature follows schedule(Tstart, tau, t)
    from Tstart to Tfinal. Raises a ValueError if the temperatures or schedule would give a run that never ends
    (more than max_steps steps). Returns lists of the protein, step number, temperature and energy at every step"""
    # check the temperature inputs, so a run always ends
    if Tstart <= 0 or Tfinal <= 0:
        raise ValueError(f"Error: temperatures must be positive, got Tstart = {Tstart}, Tfinal = {Tfinal}")
    if Tfinal >= Tstart:
        raise ValueError(f"Error: Tfinal must be lower than Tstart, got Tstart = {Tstart}, Tfinal = {Tfinal}")
    if tau <= 0:
        raise ValueError(f"Error: tau must be positive, got tau = {tau}")
    if schedule(Tstart, tau, max_steps) > Tfinal:
        raise ValueError(f"Error: the cooling schedule doesn't reach Tfinal = {Tfinal} within max_steps = {max_steps} steps; "
                         "use a smaller tau, a higher Tfinal or a larger max_steps")

    protein = Protein(initialProtein.coords, initialProtein.width, initialProtein.height) 
    stateList = [Protein(protein.coords, protein.width, protein.height)] # the protein at every step

    # track temperature and time for plotting cooling schedule
    tempList = []
    tList = []
    # track energy for plotting as well
    EList = []

    T = Tstart # temperature
    t = 0.

    tempList.append(Tstart)
    tList.append(0.)

    E = protein.calc_energy() # energy
    EList.append(E)

    while T > Tfinal:
        t += 1.0
        tList.append(t)
        T = max(schedule(Tstart, tau, t), Tfinal) # update T
        tempList.append(T)
        # propose a move
        i = randrange(len(protein.coords)) # random monomer
        dx, dy = choice([-1, 1]), choice([-1, 1]) # random diagonal step
        x, y = protein.coords[i]
        new = (x + dx, y + dy)

        inside = 0 <= new[0] < protein.width and 0 <= new[1] < protein.height
        empty = new not in protein.coords
        # must stay one step from the monomers it's bonded to (i-1 and i+1, if they exist)
        bonded_ok = True
        for j in (i - 1, i + 1):
            if 0 <= j < len(protein.coords):
                nx, ny = protein.coords[j]
                if abs(new[0] - nx) + abs(new[1] - ny) != 1:
                    bonded_ok = False

        # an invalid move counts as a rejected step
        if inside and empty and bonded_ok:
            old = protein.coords[i]
            protein.coords[i] = new
            E_new = protein.calc_energy()
            if acceptance(E_new, E, T):
                E = E_new
            else:
                protein.coords[i] = old

        EList.append(E)
        stateList.append(Protein(protein.coords, protein.width, protein.height))

    return stateList, tList, tempList, EList
