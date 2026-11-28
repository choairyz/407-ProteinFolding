"""
Protein folding final project for PHY407, by Victor Choairy - 1011464149

The Protein class: a chain of identical monomers on a 2D rectangular lattice, stored as a list of (x, y) sites in
chain order. It checks that a chain is allowed, computes its energy, and plots it.

Requirements: Matplotlib.
Outputs: none when run on its own; protein_folding_toolkit.py and annealing_functions.py use this class.
"""
import matplotlib.pyplot as plt

# lattice steps between monomers for each side
steps = {'left': (-1, 0), 'right': (1, 0), 'top': (0, 1), 'bot': (0, -1)}
epsilon = -1


class Protein:
    """A single chain of monomers on a width x height rectangular lattice, stored as a list of (x, y) sites
    in chain order: monomer i is bonded to monomers i-1 and i+1. Sites run from (0, 0) to (width-1, height-1)"""
    def __init__(self, coords: list[tuple[int, int]], width=25, height=25):
        self.coords = [tuple(site) for site in coords] # monomer coordinates
        self.width = width # lattice width
        self.height = height # lattice height
        self.is_valid() # check if specified protein arrangement is allowed
        self.energy = self.calc_energy() # protein energy: epsilon (negative) for each pair of adjacent monomers that aren't bonded

    def is_valid(self):
        """Check the chain fits the model: every monomer inside the lattice, no two on the same site, and each
        monomer one lattice step from the one before it. Raises a ValueError describing the first problem found"""
        for i, (x, y) in enumerate(self.coords):
            if not (0 <= x < self.width and 0 <= y < self.height):
                raise ValueError(f"Error: monomer {i} at {(x, y)} is outside the {self.width} x {self.height} lattice")
        if len(set(self.coords)) != len(self.coords): # sets do not count duplicates, so duplicates would be caught here
            raise ValueError("Error: at least two monomers are on the same site")
        for i in range(1, len(self.coords)):
            (x0, y0), (x1, y1) = self.coords[i - 1], self.coords[i]
            if abs(x1 - x0) + abs(y1 - y0) != 1:
                raise ValueError(f"Error: monomers {i - 1} and {i} are bonded but not one lattice step apart")
        return True

    def calc_energy(self):
        """epsilon times the number of adjacent monomer pairs that aren't bonded, counting each pair once"""
        index = {site: i for i, site in enumerate(self.coords)} # which monomer sits on each site
        contacts = 0
        for i, (x, y) in enumerate(self.coords):
            for dx, dy in steps.values():
                j = index.get((x + dx, y + dy))
                # j > i + 1 skips the bonded neighbour (j = i + 1) and counts each pair only from its lower index
                if j is not None and j > i + 1:
                    contacts += 1
        return epsilon * contacts

    def Plot(self, title="Protein", show=True):
        """Plot the chain on its lattice: bonds as lines, monomers as red dots"""
        xs = [x for x, _ in self.coords]
        ys = [y for _, y in self.coords]
        if show:
            plt.figure()
        plt.plot(xs, ys, '-o', color='red', markersize=8)
        plt.xlim(-0.5, self.width - 0.5)
        plt.ylim(-0.5, self.height - 0.5)
        ax = plt.gca()
        ax.set_xticks(range(self.width), minor=True)
        ax.set_yticks(range(self.height), minor=True)
        ax.grid(which='both', color='lightgray', linewidth=0.5)
        plt.title(rf"{title}, Energy = ${-1 * self.energy}\epsilon$")
        if show:
            plt.show()
    
