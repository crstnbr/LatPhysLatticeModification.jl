function getEmbeddedIn3D(
            unitcell :: U
        ) :: Unitcell{Site{LS,3}, B} where {LS,B,S<:AbstractSite{LS,2},U<:AbstractUnitcell{S,B}}

    # build a new list of sites
    sites_new = Site{LS,3}[
        newSite(Site{LS,3}, [point(s)[1], point(s)[2], 0.0], label(s)) for s in sites(unitcell)
    ]

    # create and return a new Unitcell
    return newUnitcell(
        Unitcell{Site{LS,3}, B},
        latticeVectors(unitcell),
        sites_new,
        bonds(unitcell)
    )
end

function getEmbeddedIn3D(
            lattice :: L
        ) :: Lattice{Site{LS,3}, B, U} where {LS,B,S<:AbstractSite{LS,2},U,L<:AbstractLattice{S,B,U}}

    # build a new list of sites
    sites_new = Site{LS,3}[
        newSite(Site{LS,3}, [point(s)[1], point(s)[2], 0.0], label(s)) for s in sites(lattice)
    ]

    # create and return a new Unitcell
    return newLattice(
        Lattice{Site{LS,3}, B, U},
        latticeVectors(lattice),
        sites_new,
        bonds(lattice),
        unitcell(lattice)
    )
end

export getEmbeddedIn3D
