function relabelBondsToString(
            lattice :: L
        ) :: Lattice{S,Bond{String,N},U} where {LS,LB,D,N,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U,L<:AbstractLattice{S,B,U}}

    # create a list of new bonds
    bondlist = Bond{String,N}[
        newBond(Bond{String,N}, from(b), to(b), string(label(b)), wrap(b))
        for b in bonds(lattice)
    ]

    # return a new lattice object
    return newLattice(
        Lattice{S,Bond{String,N},U},
        deepcopy(latticeVectors(lattice)),
        deepcopy(sites(lattice)),
        bondlist,
        deepcopy(unitcell(lattice))
    )
end

# export the relabeling function
export relabelBondsToString
