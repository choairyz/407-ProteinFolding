# protein class, etc
from __future__ import annotations
import matplotlib.pyplot as plt

# lattice steps between monomers for each side
steps = {'left': (-1, 0), 'right': (1, 0), 'top': (0, 1), 'bot': (0, -1)}
epsilon = -1


class Monomer:
    """A monomer used to build a starting shape by hand. Each side holds (neighbour, link) or None,
    where link is True for a chain bond. Turn the finished chain into a Protein with protein_from_monomers"""
    def __init__(self) -> None:
        self.left: tuple[Monomer, bool] | None = None
        self.right: tuple[Monomer, bool] | None = None
        self.top: tuple[Monomer, bool] | None = None
        self.bot: tuple[Monomer, bool] | None = None

    def addMonomer(self, location: str, link: bool, monomer: Monomer) -> None:
        """Put monomer on the given side (location) of self. Sets the neighbour on both monomers"""
        if location == 'left':
            self.left = (monomer, link)
            monomer.right = (self, link)
        elif location == 'right':
            self.right = (monomer, link)
            monomer.left = (self, link)
        elif location == 'top':
            self.top = (monomer, link)
            monomer.bot = (self, link)
        elif location == 'bot':
            self.bot = (monomer, link)
            monomer.top = (self, link)
        else:
            print(f"Error: location must be left, right, top or bot, got {location}")

    def neighbours(self, linked_only: bool = False) -> list[tuple[str, Monomer]]:
        """(side, monomer) for each occupied side of self. With linked_only, returns chained neighbors only"""
        neighbours = []
        for side, n in (('left', self.left), ('right', self.right), ('top', self.top), ('bot', self.bot)):
            if n is not None and (n[1] or not linked_only):
                neighbours.append((side, n[0]))
        return neighbours

    def chain_list(self) -> list[Monomer] | None:
        """Walk along chain bonds starting at self until a monomer with only the previous monomer in the path as
        a neighbor is reached. Returns None if the bonds don't form a single open chain"""
        chain = [self]
        prev, curr = None, self
        while True:
            bonded = curr.neighbours(linked_only=True)
            if len(bonded) > 2:
                print("Error: a monomer has more than two chain bonds, so this is not a linear chain")
                return None
            following = [m for _, m in bonded if m is not prev]
            if not following: # end of chain
                return chain
            prev, curr = curr, following[0]
            if curr in chain:
                print("Error: a closed loop is invalid")
                return None
            chain.append(curr)

    def makeChain(self) -> list[Monomer] | None:
        """Return every monomer in the chain containing self, from one end to the other (None if invalid)"""
        chain = self.chain_list() # go to one end of the chain
        if chain is None:
            return None
        return chain[-1].chain_list() # go back to the other end


def protein_from_monomers(monomers: list[Monomer], width=25, height=25) -> Protein | None:
    """Turn Monomers linked with addMonomer into a Protein, centred in a width by height lattice"""
    chain = monomers[0].makeChain()
    if chain is None:
        return None
    if set(chain) != set(monomers):
        print("Error: at least one monomer listed isn't in the chain")
        return None

    # follow the bonds one lattice step at a time, starting from (0, 0)
    sites = [(0, 0)]
    for curr, following in zip(chain, chain[1:]):
        side = next(s for s, m in curr.neighbours(linked_only=True) if m is following)
        dx, dy = steps[side]
        x, y = sites[-1]
        sites.append((x + dx, y + dy))

    # shift the chain so it sits in the middle of the lattice
    xs = [x for x, _ in sites]
    ys = [y for _, y in sites]
    shift_x = (width - (max(xs) - min(xs) + 1)) // 2 - min(xs)
    shift_y = (height - (max(ys) - min(ys) + 1)) // 2 - min(ys)
    return Protein([(x + shift_x, y + shift_y) for x, y in sites], width, height)


class Protein:
    """A single chain of monomers on a width x height rectangular lattice, stored as a list of (x, y) sites
    in chain order: monomer i is bonded to monomers i-1 and i+1. Sites run from (0, 0) to (width-1, height-1)"""
    def __init__(self, coords: list[tuple[int, int]], width: int = 25, height: int = 25):
        self.coords = [tuple(site) for site in coords]
        self.width = width
        self.height = height
        self.is_valid() # prints the problem if the starting shape breaks the model
        self.energy = self.calc_energy()

    def is_valid(self) -> bool:
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

    def calc_energy(self) -> float:
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

    def Plot(self, title="Protein"):
        """Plot the chain on its lattice: bonds as lines, monomers as red dots"""
        xs = [x for x, _ in self.coords]
        ys = [y for _, y in self.coords]
        plt.figure()
        plt.plot(xs, ys, '-o', color='black', markerfacecolor='red', markersize=8)
        plt.xlim(-0.5, self.width - 0.5)
        plt.ylim(-0.5, self.height - 0.5)
        plt.title(f"{title}, E = {self.energy}")
        plt.show()
