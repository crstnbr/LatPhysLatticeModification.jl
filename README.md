# LatPhysLatticeModification.jl

Lattice modification module of [`LatticePhysics.jl`](https://github.com/janattig/LatticePhysics.jl)



## Contents

Modification of lattices and unitcells. Currently implemented:
1.  Relabeling to bipartite
2.  Relabeling to integer / string
3.  Adding and removing sites


## Installation

You can install the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to first install the unregistered
dependencies manually with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysBase.jl"
```
to finally install the main package with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysLatticeModification.jl"
```
