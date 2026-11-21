# protein class, etc

import matplotlib.pyplot as plt

# lattice steps between monomers for each side
steps = {'left': (-1, 0), 'right': (1, 0), 'top': (0, 1), 'bot': (0, -1)}

class Monomer:
    """A monomer on a 2D rectangular lattice. Each side holds (neighbour, link) or None,
    where link is True for a chain bond and False for a non-bonded lattice neighbour"""
    def __init__(self) -> None:
        self.left: tuple[Monomer, bool] | None = None
        self.right: tuple[Monomer, bool] | None = None
        self.top: tuple[Monomer, bool] | None = None
        self.bot: tuple[Monomer, bool] | None = None
        self.coords: tuple | None = None

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
        neighbors = []
        for side, n in (('left', self.left), ('right', self.right), ('top', self.top), ('bot', self.bot)):
            if n is not None and (n[1] or not linked_only):
                neighbors.append((side, n[0]))
        return neighbors

    def chain_list(self) -> list[Monomer]:
        """Walk along chain bonds starting at self until a monomer with only the previous monomer in the path as 
        a neighbor is reached"""
        chain = [self]
        prev, curr = None, self
        while True:
            bonded = curr.neighbours(linked_only=True)
            if len(bonded) > 2:
                print("A monomer has more than two chain bonds, so this is not a linear chain")
                return None
            next = [m for _, m in bonded if m is not prev]
            if not next: # end of chain
                return chain
            prev, curr = curr, next[0]
            if curr in chain:
                print("A closed loop is invalid")
                return None
            chain.append(curr)

    def makeChain(self) -> list[Monomer]:
        """Return every monomer in the chain containing self, from one end to the other"""
        end = self.chain_list()[-1] # go to one end of the chain
        return end.chain_list() # then go back to the other end


class Protein:
    """A single chain of monomers on a 2D rectangular lattice"""
    def __init__(self, monomers: list[Monomer]):
        self.monomers = monomers

    def set_coords(self):
        """Set coordinates for each monomer in the protein chain, with the left end of the chain at (0, 0)"""
        chain = self.monomers[0].makeChain()
        if set(chain) != set(self.monomers):
            print("At least one monomer listed isn't in the chain")
            return

        # walk the chain in  first and scale by the spacings at the end
        # so that rounding error in the spacings can't hide two monomers sharing a site
        sites = [(0, 0)]
        for curr, next in chain, chain[1:]:
            side = next(s for s, m in curr.neighbours(linked_only=True) if m is next)
            dx, dy = steps[side]
            x, y = sites[-1]
            sites.append((x + dx, y + dy))

        if len(set(sites)) != len(sites):
            print("Error: at leasdt two monomers are in the same place")
            return
        
    def Plot(self):
        """Plot protein chain"""
        self.set_coords()
        for m in self.monomers:
            xpoint, ypoint = m.coords
            for side, n in m.neighbors(linked_only=True):
                if side == 'left':
                    plt.axhline(y=ypoint, xmin=n.coords[0], xmax=m.coords[0])
                if side == 'right':
                    plt.axhline(y=ypoint, xmin=m.coords[0], xmax=n.coords[0])
                if side == 'top':
                    plt.axvline(x=xpoint, ymin=m.coords[1], ymax=n.coords[1])
                if side == 'bot':
                    plt.axvline(x=xpoint, ymin=n.coords[1], ymax=m.coords[1])
            plt.plot(xpoint, ypoint, '.', color='red')
        plt.show()

