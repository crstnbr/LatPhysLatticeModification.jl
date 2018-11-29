#################################################################################
#
#   FUNCTIONS FOR ADDING / REMOVING SITES IN LATTICES AND UNITCELLS
#
#   STRUCTURE OF THE FILE:
#   - removing single sites
#       - from unitcells
#       - from lattices
#   - removing multiple sites
#       - from unitcells
#       - from lattices
#
#   - adding single sites
#       - to unitcells
#       - to lattices
#   - adding multiple sites
#       - to unitcells
#       - to lattices
#
################################################################################






################################################################################
#
#   REMOVING SINGLE SITES
#       - from unitcells
#       - from lattices
#
################################################################################


# Remove a single site from a unitcell (and return that)
function removeSite!(
        unitcell :: U,
        index    :: Integer
    ) :: S where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # store the site to be removed
    s = sites(unitcell)[index]

    # remove the site from the list and pass the new list again to the unitcell
    sites!(unitcell, deleteat!(sites(unitcell), index))

    # get all bonds that are not connected to the removed site
    bond_list = filter(b -> to(b) != index && from(b) != index, bonds(unitcell))
    # relabeling of bond indices
    for b in bond_list
        # set the new indices
        from!(b, from(b)>index ? from(b)-1 : from(b))
        to!(b,   to(b)>index   ? to(b)-1   : to(b))
    end
    # set the bond list in the unitcell object
    bonds!(unitcell, bond_list)
    # return the removed site
    return s
end

# Remove a single site from a lattice (and return that)
function removeSite!(
        lattice :: L,
        index   :: Integer
    ) :: S where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B},
        DL,NL,LSL,LBL,SL<:AbstractSite{LSL,DL}, BL<:AbstractBond{LBL,NL}, L<:AbstractLattice{SL,BL,U}}

    # store the site to be removed
    s = sites(lattice)[index]

    # remove the site from the list and pass the new list again to the lattice
    sites!(lattice, deleteat!(sites(lattice), index))

    # get all bonds that are not connected to the removed site
    bond_list = filter(b -> to(b) != index && from(b) != index, bonds(lattice))
    # relabeling of bond indices
    for b in bond_list
        # set the new indices
        from!(b, from(b)>index ? from(b)-1 : from(b))
        to!(b,   to(b)>index   ? to(b)-1   : to(b))
    end
    # set the bond list in the lattice object
    bonds!(lattice, bond_list)
    # return the removed site
    return s
end


# export the function
export removeSite!
