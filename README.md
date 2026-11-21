# 407-ProteinFolding

Final project for PHY407 (Fall 2026), topic 3: Protein Folding.

A Markov Chain Monte Carlo simulation with simulated annealing that folds a simplified 2D protein, starting from a shape you choose.

## The model

- A protein is a chain of N identical monomers on a rectangular `width` × `height` lattice (25 × 25 by default).
- Monomer `i` is bonded to monomers `i-1` and `i+1`. Bonded monomers are always one lattice step apart, and no two monomers share a site.
- Every pair of monomers that are next to each other on the lattice but **not** bonded adds ε = −1 to the energy, so E = ε × (number of such contacts).

Temperatures are in units of |ε|/k<sub>B</sub>.

## How the simulation works

Each Monte Carlo step:

1. Lower the temperature: T = T<sub>start</sub> · e<sup>−t/τ</sup>, where t is the step number.
2. Pick a random monomer and a random diagonal step (±1, ±1). Keep picking until the move is allowed: the new site has to be inside the lattice, empty, and one step from the monomer's bonded neighbours. For a monomer in the middle of the chain, this is a "corner flip" across a bend. For an end monomer, it's a 90° turn around its neighbour.
3. Accept or reject the move with the Metropolis rule: always accept if the energy doesn't go up, otherwise accept with probability e<sup>−ΔE/T</sup>. A rejected move is undone.
4. Record T, E and a copy of the protein.

The run stops when T falls below T<sub>final</sub>, which takes about τ · ln(T<sub>start</sub>/T<sub>final</sub>) steps.

## Files

| File | Contents |
|---|---|
| [`src/protein_folding_toolkit.py`](src/protein_folding_toolkit.py) | The settings, `build_protein(moves)` for the starting shape, `plot_states()` for the animation. Running it runs the simulation. |
| [`src/protein.py`](src/protein.py) | The `Protein` class: monomer positions, `is_valid()`, `calc_energy()` and `Plot()`. |
| [`src/annealing_functions.py`](src/annealing_functions.py) | `acceptance()` (Metropolis rule), `cooledT()` (cooling schedule) and `discreteAnnealing()` (the simulation loop). |

## Running it

You need Python 3, NumPy and Matplotlib. To install them:

```bash
pip install -r requirements.txt
```

Then run:

```bash
python src/protein_folding_toolkit.py
```

A window opens and shows the protein every 100 steps.

### Settings

The settings are at the top of `protein_folding_toolkit.py`:

```python
Tstart = 100   # initial temperature
Tfinal = 0.1   # final temperature
tau = 1000     # decay constant, in Monte Carlo steps
```

### Starting shape

Give the starting shape as a string of moves, one letter per bond: **R**ight, **L**eft, **U**p, **D**own. You start at one end of the chain and say which way each next monomer goes. A string of length N − 1 gives N monomers, centred in the lattice.

```python
initialProtein, N = build_protein('RURRRURRRRRRRURRRRRRRR')   # 23 monomers
initialProtein, N = build_protein('R' * 14)                  # straight chain of 15
```

A string can't produce a broken chain, since every letter is one lattice step. If the shape runs into itself or leaves the lattice, `Protein` prints what's wrong.
