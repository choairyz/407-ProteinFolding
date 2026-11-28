"""
Protein folding final project for PHY407, by Victor Choairy - 1011464149

Folds a simplified 2D lattice protein with a Markov chain Monte Carlo simulation with simulated annealing,
starting from a shape chosen by the user. This is the file to run.

Requirements: Python 3 with NumPy and Matplotlib, and protein.py and annealing_functions.py in the same folder.
Inputs: the settings below, and the starting shape passed to build_protein as a string of moves.
Outputs: figures of the starting and final structures, the cooling schedule (temperature against Monte Carlo
step) and the energy against Monte Carlo step. No files are written.
"""
from random import seed
import matplotlib.pyplot as plt
from annealing_functions import discreteAnnealing, exponential_cooling, linear_cooling, inverse_cooling
from protein import Protein

# settings
Tstart = 100 # initial temperature, in units of |epsilon|/k_B
Tfinal = 0.1 # final temperature: the run stops once T reaches this
tau = 1000 # time constant of the cooling schedule, in Monte Carlo steps
schedule_name = 'exponential' # or 'linear' or 'inverse'
seed_num = 1 # seed for the random number generator

if schedule_name == 'exponential':
    schedule = exponential_cooling
elif schedule_name == 'linear':
    schedule = linear_cooling
elif schedule_name == 'inverse':
    schedule = inverse_cooling
else:
    raise ValueError('Not a valid cooling schedule. Has to be one of: exponential, linear or inverse')



def build_protein(moves, width=25, height=25):
    """Build a Protein chain from a string of moves, one letter per bond: R = right, L = left, U = up, D = down,
    so a string of N-1 letters gives N monomers. The chain is centred in the width x height lattice.
    Returns the Protein and its number of monomers N"""
    step = {'R': (1, 0), 'L': (-1, 0), 'U': (0, 1), 'D': (0, -1)}
    sites = [(0, 0)]
    for letter in moves.upper():
        if letter not in step:
            raise ValueError(f"Error: moves can only use R, L, U or D, got {letter}")
        dx, dy = step[letter]
        x, y = sites[-1]
        sites.append((x + dx, y + dy))

    # shift the chain so it sits in the middle of the lattice
    xs = [x for x, _ in sites]
    ys = [y for _, y in sites]
    shift_x = (width - (max(xs) - min(xs) + 1)) // 2 - min(xs)
    shift_y = (height - (max(ys) - min(ys) + 1)) // 2 - min(ys)
    return Protein([(x + shift_x, y + shift_y) for x, y in sites], width, height), len(sites)


def plot_states(stateList, skip=100):
    """Animate the protein in one window, drawing every skip-th Monte Carlo step and ending on the final state"""
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
    seed(seed_num)
    initialProtein, N = build_protein('RUUULURRRDDRDRR')
    initialProtein.Plot(title='Protein: Initial State')
    stateList, tList, tempList, EList = discreteAnnealing(initialProtein, tau, Tstart, Tfinal, schedule)
    #plot_states(stateList)
    stateList[-1].Plot(title='Protein: Final State')

    plt.plot(tList, tempList, color='orange')
    plt.xlabel('Monte Carlo step')
    plt.ylabel(r'Temperature $T$ ($|\epsilon|/k_B$)')
    plt.title(f'Cooling Schedule ({schedule_name})')
    plt.show()

    plt.plot(tList, EList, color='red')
    plt.xlabel('Monte Carlo step')
    plt.ylabel(r'Energy $E$ ($|\epsilon|$)')
    plt.title('Protein energy during simulated annealing')
    plt.show()