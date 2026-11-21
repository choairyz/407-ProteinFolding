# protein class, etc
from __future__ import annotations
import matplotlib.pyplot as plt

# lattice steps between monomers for each side
steps = {'left': (-1, 0), 'right': (1, 0), 'top': (0, 1), 'bot': (0, -1)}
epsilon = -1


class Protein:
    """A single chain of monomers on a width x height rectangular lattice, stored as a list of (x, y) sites
    in chain order: monomer i is bonded to monomers i-1 and i+1. Sites run from (0, 0) to (width-1, height-1)"""
    def __init__(self, coords: list[tuple[int, int]], width=25, height=25):
        self.coords = [tuple(site) for site in coords]
        self.width = width
        self.height = height
        self.is_valid() 
        self.energy = self.calc_energy()

    def is_valid(self):
        """Check the chain fits the model: every monomer inside the lattice, no two on the same site, and each
        monomer one lattice step from the one before it. Prints the first problem found"""
        for i, (x, y) in enumerate(self.coords):
            if not (0 <= x < self.width and 0 <= y < self.height):
                print(f"Error: monomer {i} at {(x, y)} is outside the {self.width} x {self.height} lattice")
                return False
        if len(set(self.coords)) != len(self.coords):
            print("Error: at least two monomers are on the same site")
            return False
        for i in range(1, len(self.coords)):
            (x0, y0), (x1, y1) = self.coords[i - 1], self.coords[i]
            if abs(x1 - x0) + abs(y1 - y0) != 1:
                print(f"Error: monomers {i - 1} and {i} are bonded but not one lattice step apart")
                return False
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
        plt.title(f"{title}, E = {self.energy}")
        if show:
            plt.show()
    
