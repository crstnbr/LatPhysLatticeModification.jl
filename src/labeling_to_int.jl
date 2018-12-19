function relabelSitesWithIntIndex(
            unitcell :: U
        ) :: Unitcell{Site{Int64,D},B} where {LS,LB,D,N,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # create a list of new sites
    sitelist = Site{Int64,D}[
        newSite(Site{Int64,D}, point(site(unitcell,i)), i)
        for i in 1:numSites(unitcell)
    ]

    # return a new lattice object
    return newUnitcell(
        Unitcell{Site{Int64,D},B},
        deepcopy(latticeVectors(unitcell)),
        sitelist,
        deepcopy(bonds(unitcell))
    )
end

# export the relabeling function
export relabelSitesWithIntIndex
