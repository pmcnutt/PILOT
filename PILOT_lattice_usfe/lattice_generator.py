#!/usr/bin/env python3
"""
lattice_generator.py -- build a random solid-solution (RSS) atomic lattice.

This is the lattice-building stage of the PILOT USFE workflow.  It replaces the
earlier OTIS "Order Through Informed Swapping" generator: because PILOT models
*random* solid solutions only, there is no atom swapping, no bond counting and
no short-range-order (SRO) matrix.  Atoms are placed on the grid in the exact
target proportions and then shuffled, which is the definition of an ideal RSS.

Two lattices are supported, each tied to one slip system downstream:

    BCC                -> {110} slip   (orthogonal-in-index bcc grid)
    CUSTOM_NON_ORTHO   -> {112} slip   (non-orthogonal grid for the {112} interp)

Input file format (whitespace-delimited, one item per line; trailing "#..."
comments are ignored):

    BCC                 # lattice type: BCC or CUSTOM_NON_ORTHO
    32 32 32            # N1 N2 N3  (primitive cells in each direction)
    2.86                # nearest-neighbour distance (lattice parameter)
    2                   # number of components
    Ta W                # element names
    0.5 0.5             # concentrations (renormalised if they don't sum to 1)
    my_lattice          # output file name

Run:

    python3 lattice_generator.py <input_file>

The output file matches the format the interpolation scripts and the viewer
expect: an extended-XYZ-style header (atom count, then a Lattice="..." / pbc
line) followed by one "element x y z" line per site, in the site ordering
idx = N1*N2*k + N1*j + i.
"""

import sys
import numpy as np
import numpy.random as rand


# Primitive lattice vectors for each supported lattice type.  These are kept
# byte-for-byte identical to the previous generator so that existing USFE files,
# interpolation neighbour tables and the viewer's projection all stay valid.
LATTICE_VECTORS = {
    "BCC": (
        np.array([1, 1, -1]) / np.sqrt(3),
        np.array([-1, 1, 1]) / np.sqrt(3),
        np.array([1, -1, 1]) / np.sqrt(3),
    ),
    "CUSTOM_NON_ORTHO": (
        np.array([0, 1, -1]) / np.sqrt(2),
        np.array([-1, 1, 1]) / np.sqrt(3),
        np.array([1, 1, 0]) / np.sqrt(2),
    ),
}


def _clean(line):
    """Strip a trailing '#...' comment and surrounding whitespace."""
    return line.split("#", 1)[0].strip()


def read_input(path):
    """Parse a generator input file into its fields."""
    with open(path) as f:
        raw = [f.readline() for _ in range(7)]
    if any(r == "" for r in raw):
        raise ValueError("Input file ended early; expected 7 non-empty lines.")

    lattice_type = _clean(raw[0]).upper()
    if lattice_type not in LATTICE_VECTORS:
        raise ValueError("Lattice type must be BCC or CUSTOM_NON_ORTHO "
                         f"(got {lattice_type!r}).")

    cell = _clean(raw[1]).split()
    if len(cell) < 3:
        raise ValueError("Second line must give N1 N2 N3.")
    N1, N2, N3 = (int(x) for x in cell[:3])
    if min(N1, N2, N3) <= 0:
        raise ValueError("N1 N2 N3 must all be positive.")

    lattice_parameter = float(_clean(raw[2]).split()[0])

    components = int(_clean(raw[3]).split()[0])
    if components <= 0:
        raise ValueError("Number of components must be positive.")

    labels = _clean(raw[4]).split()
    if len(labels) < components:
        raise ValueError("Fewer element names than components.")
    labels = labels[:components]

    conc = np.array(_clean(raw[5]).split()[:components], dtype=float)
    if len(conc) < components:
        raise ValueError("Fewer concentrations than components.")
    if abs(conc.sum() - 1.0) > 0.1:
        raise ValueError(f"Concentrations sum to {conc.sum():.3f}; must be ~1.")
    conc = conc / conc.sum()

    output_file = _clean(raw[6]).split()[0]

    return dict(lattice_type=lattice_type, cell=(N1, N2, N3),
                lattice_parameter=lattice_parameter, components=components,
                labels=labels, conc=conc, output_file=output_file)


def build_rss(cell, conc, seed=None):
    """Return an integer array (one component index per site) for a random
    solid solution with the requested concentrations, in idx order."""
    N1, N2, N3 = cell
    ntot = N1 * N2 * N3
    components = len(conc)

    lattice = np.zeros(ntot, dtype=int)
    for c in range(components):
        a = int(np.sum(conc[:c]) * ntot)
        b = int(np.sum(conc[:c + 1]) * ntot)
        lattice[a:b] = c
    # any sites left as 0 by integer rounding are already component 0; fine.

    rng = rand.default_rng(seed) if seed is not None else rand
    rng.shuffle(lattice)
    return lattice


def site_positions(cell, a1, a2, a3):
    """Real-space (unscaled) position of every site, in idx order."""
    N1, N2, N3 = cell
    pos = np.zeros((N1 * N2 * N3, 3))
    for i in range(N1):
        for j in range(N2):
            for k in range(N3):
                idx = N1 * N2 * k + N1 * j + i
                pos[idx] = i * a1 + j * a2 + k * a3
    return pos


def write_lattice(path, cell, lattice, labels, pos, lattice_parameter,
                  a1, a2, a3):
    """Write the extended-XYZ-style lattice file the rest of PILOT reads."""
    N1, N2, N3 = cell
    lp = lattice_parameter
    with open(path, "w") as f:
        f.write("%d\n" % (N1 * N2 * N3))
        f.write('Lattice="%0.8f %0.8f %0.8f ' % tuple(a1 * N1 * lp))
        f.write('%0.8f %0.8f %0.8f ' % tuple(a2 * N2 * lp))
        f.write('%0.8f %0.8f %0.8f" pbc="T T T"\n' % tuple(a3 * N3 * lp))
        for idx in range(N1 * N2 * N3):
            f.write("%s %f %f %f\n" % (labels[lattice[idx]],
                                       pos[idx, 0] * lp,
                                       pos[idx, 1] * lp,
                                       pos[idx, 2] * lp))


def generate(cfg, seed=None):
    """Build and write one RSS lattice from a parsed-input dict."""
    a1, a2, a3 = LATTICE_VECTORS[cfg["lattice_type"]]
    N1, N2, N3 = cfg["cell"]
    print("Generating %s random solid solution: %dx%d x%d = %d atoms"
          % (cfg["lattice_type"], N1, N2, N3, N1 * N2 * N3))
    print("Concentrations: "
          + ", ".join("%s=%.4f" % (lab, c)
                      for lab, c in zip(cfg["labels"], cfg["conc"])))

    lattice = build_rss(cfg["cell"], cfg["conc"], seed=seed)
    pos = site_positions(cfg["cell"], a1, a2, a3)
    write_lattice(cfg["output_file"], cfg["cell"], lattice, cfg["labels"],
                  pos, cfg["lattice_parameter"], a1, a2, a3)

    # report the realised composition (differs slightly from target by rounding)
    realised = np.array([np.mean(lattice == c) for c in range(cfg["components"])])
    print("Realised:       "
          + ", ".join("%s=%.4f" % (lab, r)
                      for lab, r in zip(cfg["labels"], realised)))
    print("Lattice written to %s" % cfg["output_file"])
    return lattice


def main():
    if len(sys.argv) < 2:
        print("usage: python3 lattice_generator.py <input_file>")
        return 1
    cfg = read_input(sys.argv[1])
    generate(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
