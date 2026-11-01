# graphene-grain-boundary

> **Legacy research code (2019), not maintained.** Kept for reference.

An early Python script that builds a graphene grain-boundary structure: it relaxes Voronoi seed points with Lloyd's algorithm, inspired by the algorithm of Ophus et al., "Large-Scale Experimental and Theoretical Study of Graphene Grain Boundary Structures" (2015).

- `gb_shiftperiter_modify.py`: the relaxation script. It reads its starting seed coordinates from a text file (`edge1coord150k.txt`), which is not included, and writes the relaxed coordinates to a new file.
- `1.png`, `GB200k.png`: the initial grain-boundary structure and the structure after 200,000 steps of Lloyd's algorithm.

Requires Python 3 with NumPy, SciPy, Matplotlib and Shapely.
