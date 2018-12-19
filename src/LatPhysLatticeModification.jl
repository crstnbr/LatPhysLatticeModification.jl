################################################################################
#
#   module LatticePhysics_LatticeModification
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> MODIFICATION OF LATTICES COMPONENTS
#       - Adding / Removing Connections
#       - Adding / Removing sites
#
#   --> RELABELING BONDS / SITES
#       - relabeling bonds
#       - relabeling sites
#
#   --> SPATIAL MODIFICATIONS
#       - Scaling  in space
#       - Rotating in space
#       - Shifting in space
#       - Lattice vector optimization (relabeling unitcells)
#
################################################################################


# Start of module
module LatPhysLatticeModification


# Using libraries
using LatPhysBase




# ADDING / REMOVING SITES
include("site_add_remove.jl")

# TRANSFORMATIONS IN REAL SPACE
include("real_space_transformation.jl")



# RELABELING

# setting labels of sites to bipartite
include("labeling_bipartite.jl")

# relabel all bonds to string
include("labeling_to_string.jl")

# relabel all sties to int
include("labeling_to_int.jl")


# End of module
end
