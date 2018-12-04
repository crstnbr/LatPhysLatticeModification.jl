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





################################################################################
#
#   REMOVING MULTIPLE SITES
#       - from unitcells
#       - from lattices
#
################################################################################


# Remove multiple sites from a unitcell (and return them)
function removeSite!(
        unitcell :: U,
        index    :: Integer ...
    ) :: Vector{S} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # get the list of sites
    sites_to_return = S[site(unitcell,i) for i in index]

    # get all sites of the unitcell
    site_list = sites(unitcell)
    # build a index list
    indices_new = collect(1:length(site_list))
    # fill the new index list at positions of removed sites
    for i in index
        indices_new[i] = -1
    end
    # fill the new list successively from smallest to largest
    if indices_new[1] == -1
        indices_new[1] = 0
    end
    for i in 2:length(indices_new)
        if indices_new[i] == -1
            indices_new[i] = indices_new[i-1]
        else
            indices_new[i] = indices_new[i-1] + 1
        end
    end
    # get a new site list
    site_list = S[site_list[i] for i in 1:length(site_list) if !(i in index)]
    # register new list of sites within the unitcell object
    sites!(unitcell, site_list)


    # get all bonds that are not connected to any removed site
    bond_list = filter(b -> !(to(b) in index || from(b) in index), bonds(unitcell))
    # relabeling of bond indices
    for b in bond_list
        # set the new indices
        from!(b, indices_new[from(b)])
        to!(b,   indices_new[to(b)])
    end
    # set the bond list in the unitcell object
    bonds!(unitcell, bond_list)


    # return the removed sites
    return sites_to_return
end

# Remove multiple sites from a lattice (and return them)
function removeSite!(
        lattice :: L,
        index   :: Integer ...
    ) :: Vector{S} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B},
        DL,NL,LSL,LBL,SL<:AbstractSite{LSL,DL}, BL<:AbstractBond{LBL,NL}, L<:AbstractLattice{SL,BL,U}}

    # get the list of sites
    sites_to_return = S[site(lattice,i) for i in index]

    # get all sites of the lattice
    site_list = sites(lattice)
    # build a index list
    indices_new = collect(1:length(site_list))
    # fill the new index list at positions of removed sites
    for i in index
        indices_new[i] = -1
    end
    # fill the new list successively from smallest to largest
    if indices_new[1] == -1
        indices_new[1] = 0
    end
    for i in 2:length(indices_new)
        if indices_new[i] == -1
            indices_new[i] = indices_new[i-1]
        else
            indices_new[i] = indices_new[i-1] + 1
        end
    end
    # get a new site list
    site_list = SL[site_list[i] for i in 1:length(site_list) if !(i in index)]
    # register new list of sites within the lattice object
    sites!(lattice, site_list)


    # get all bonds that are not connected to any removed site
    bond_list = filter(b -> !(to(b) in index || from(b) in index), bonds(lattice))
    # relabeling of bond indices
    for b in bond_list
        # set the new indices
        from!(b, indices_new[from(b)])
        to!(b,   indices_new[to(b)])
    end
    # set the bond list in the lattice object
    bonds!(lattice, bond_list)


    # return the removed sites
    return sites_to_return
end








################################################################################
#
#   ADDING SINGLE SITES
#       - to unitcells
#       - to lattices
#
################################################################################

# Add a single site to a unitcell
function addSite!(
            unitcell :: U,
            position :: Vector{<:Real},
            label    :: L = getDefaultLabel(L)
        ) :: S where {L,D,B,S<:AbstractSite{L,D},U<:AbstractUnitcell{S,B}}

    # create a new site
    site = newSite(S, position, label)

    # obtain the list of sites in the unitcell
    site_list = sites(unitcell)
    # push to the list of sites
    push!(site_list, site)
    # set the list of sites in the unitcell
    sites!(unitcell, sites_list)

    # return the new site
    return site
end


# Add a single site to a lattice
function addSite!(
            lattice  :: LA,
            position :: Vector{<:Real},
            label    :: L = getDefaultLabel(L)
        ) :: S where {L,D,B,S<:AbstractSite{L,D},U,LA<:AbstractLattice{S,B,U}}

    # create a new site
    site = newSite(S, position, label)

    # obtain the list of sites in the lattice
    site_list = sites(lattice)
    # push to the list of sites
    push!(site_list, site)
    # set the list of sites in the lattice
    sites!(lattice, sites_list)

    # return the new site
    return site
end


# export the function
export addSite!
