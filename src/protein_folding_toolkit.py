# header

from annealing_functions import *
from protein import *

Tstart = 100 # initial temperature
Tfinal = 0.1 # final temperature
tau = 1000 # decay constant

def build_protein(moves, width=25, height=25):
    """Build a Protein chain from a string of moves, one letter per bond: R = right, L = left, U = up, D = down. 
    The chain is centred in the lattice"""
    N = len(moves)
    step = {'R': (1, 0), 'L': (-1, 0), 'U': (0, 1), 'D': (0, -1)}
    sites = [(0, 0)]
    for letter in moves.upper():
        if letter not in step:
            print(f"Error: moves can only use R, L, U or D, got {letter}")
            return None
        dx, dy = step[letter]
        x, y = sites[-1]
        sites.append((x + dx, y + dy))

    # shift the chain so it sits in the middle of the lattice
    xs = [x for x, _ in sites]
    ys = [y for _, y in sites]
    shift_x = (width - (max(xs) - min(xs) + 1)) // 2 - min(xs)
    shift_y = (height - (max(ys) - min(ys) + 1)) // 2 - min(ys)
    return Protein([(x + shift_x, y + shift_y) for x, y in sites], width, height), N

def plot_states(stateList, skip=100):
    """Plots the protein structure every 100 steps during simulated annealing."""
    plt.figure()
    shown = list(range(0, len(stateList), skip))
    if shown[-1] != len(stateList) - 1:
        shown.append(len(stateList) - 1) # always end on the final state
    for step in shown:
        plt.cla() # clear the last frame
        stateList[step].Plot(f"Step {step}", show=False)
        plt.pause(0.1) # draw it and wait 0.1 s
    plt.show()

if __name__ == "__main__":
    initialProtein, N = build_protein('RUUULURRRDDRDRR')
    initialProtein.Plot(title='Protein: Initial State')
    stateList, tList, tempList, EList = discreteAnnealing(initialProtein, tau, Tstart, Tfinal)
    #plot_states(stateList)
    stateList[-1].Plot(title='Protein: Final State')
    
    plt.plot(tList, tempList, color='orange')
    plt.xlabel('Monte Carlo step')
    plt.ylabel('Temperature')
    plt.title('Cooling Schedule')
    plt.show()

    plt.plot(tList, EList, color='red', ms=1)
    plt.xlabel('Monte Carlo step')
    plt.ylabel('Protein chain energy')
    plt.title('Protein energy during simulated annealing')
    plt.show()




