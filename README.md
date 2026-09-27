# Exact Diagonalization of the Hubbard Model

I decided to re-explore a [challenging problem](http://qmcchem.ups-tlse.fr/files/caffarel/Hub_Inf_PRL_1994.pdf) first encountered in graduate school as an excuse to experiment with modern web technologies, specifically [WebAssembly](https://webassembly.org/) and [Serverless](https://serverless.com).

Because a large part of the challenge is memory, I started by implementing a version of [this stackoverflow answer](https://stackoverflow.com/a/36345790/8479938) in [AssemblyScript](https://github.com/AssemblyScript/assemblyscript).  I wasn't impressed by the compiler output so I hand tuned the [.wat](http://webassembly.github.io/spec/core/text/index.html) directly because why not.

That experiment, and a later Julia port, now live in [legacy/](legacy/).  The current implementation is a Rust command line application.

## Usage

```sh
cargo run --release -- --sites 16 --up 8 --down 8 -t 1 -U 4
```

```
lattice    16 sites, tilt [4, 0]
basis      1297076 states, 128 symmetries, 19.8 MiB for 2 vectors
lanczos    79 iterations, residual 8.7e-9
energy     -13.621854821163 (-0.851365926323 per site)
```

`--sites` must be a sum of two squares $u^2 + v^2 \le 32$: the lattice is the square with sides $(u, v)$ and $(-v, u)$ with periodic boundary conditions.  `--eigenvector` also builds the ground state to check it and report double occupancy.

The search is restricted to one symmetry sector, chosen with `--momentum gamma|m` (crystal momentum $(0, 0)$ or $(\pi, \pi)$) and `--irrep a1|a2|b1|b2` (point group irrep, B₁ being $d_{x^2-y^2}$); the default is the fully symmetric `gamma a1`.  `--all-sectors` tries each of them in turn and reports the lowest, and `--full-basis` ignores symmetries altogether.

## How memory is kept small

* The Hamiltonian is never stored.  Each row of $y \leftarrow Hx - \beta y$ is generated on the fly and only touches its own entry of $y$, so plain Lanczos needs just two vectors (three to replay the recurrence for the eigenvector).
* States are symmetrized under translations and the point group (rotations always, reflections when the cluster is not chiral).  Symmetries act on the spin up configuration first, so bookkeeping scales with $\binom{N}{N_\uparrow}$ rather than with the dimension.
* Configurations are ranked with the combinatorial number system split into two lookup tables (H. Q. Lin, Phys. Rev. B **42**, 6561).

## Tests

```sh
cargo test
```

compares against dense diagonalization of small clusters, both in the full basis and projected onto every symmetry sector, against exact noninteracting and atomic limits, and against published exact energies of 2×2, 3×3 and 4×4 clusters (H. Shi and S. Zhang, Phys. Rev. B **88**, 125132).
