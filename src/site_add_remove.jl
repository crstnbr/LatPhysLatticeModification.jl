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
#   - removing disconnected sites
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
    ) :: SL where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B},
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
        index    :: Integer...
    ) :: Vector{S} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # return the respective function
    return removeSite!(unitcell, unique(index))
end

# Remove multiple sites from a unitcell (and return them)
function removeSite!(
        unitcell :: U,
        index    :: Vector{<:Integer}
    ) :: Vector{S} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # get the list of sites
    sites_to_return = S[site(unitcell,i) for i in index]

    # get all sites of the unitcell
    site_list = sites(unitcell)
    # build a index list
    indices_new = collect(1:length(site_list))
    index_gets_removed = Bool[false for i in 1:numSites(unitcell)]
    # fill the new index list at positions of removed sites
    for i in index
        indices_new[i] = -1
        index_gets_removed[i] = true
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
    site_list = S[site_list[i] for i in 1:length(site_list) if !index_gets_removed[i]]
    # register new list of sites within the unitcell object
    sites!(unitcell, site_list)


    # get all bonds that are not connected to any removed site
    bond_list = filter(b -> !(index_gets_removed[to(b)] || index_gets_removed[from(b)]), bonds(unitcell))
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
    ) :: Vector{SL} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B},
        DL,NL,LSL,LBL,SL<:AbstractSite{LSL,DL}, BL<:AbstractBond{LBL,NL}, L<:AbstractLattice{SL,BL,U}}

    # call the respective funtion
    return removeSite!(lattice, unique(index))
end

# Remove multiple sites from a lattice (and return them)
function removeSite!(
        lattice :: L,
        index   :: Vector{<:Integer}
    ) :: Vector{SL} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B},
        DL,NL,LSL,LBL,SL<:AbstractSite{LSL,DL}, BL<:AbstractBond{LBL,NL}, L<:AbstractLattice{SL,BL,U}}

    # get the list of sites
    sites_to_return = SL[site(lattice,i) for i in index]

    # get all sites of the lattice
    site_list = sites(lattice)
    # build a index list
    indices_new = collect(1:length(site_list))
    index_gets_removed = Bool[false for i in 1:numSites(lattice)]
    # fill the new index list at positions of removed sites
    for i in index
        indices_new[i] = -1
        index_gets_removed[i] = true
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
    site_list = SL[site_list[i] for i in 1:length(site_list) if !index_gets_removed[i]]
    # register new list of sites within the lattice object
    sites!(lattice, site_list)


    # get all bonds that are not connected to any removed site
    bond_list = filter(b -> !(index_gets_removed[to(b)] || index_gets_removed[from(b)]), bonds(lattice))
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
#   REMOVING DISCONNECTED SITES
#       - from unitcells
#       - from lattices
#
################################################################################

# Remove multiple disconnected sites (from site with index) from a unitcell (and return them)
function removeDisconnectedSites!(
        unitcell :: U,
        index    :: Integer = 1
    ) :: Vector{S} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    #######################################
    # label all sites by Hoshen Kopelman
    #######################################

    # label variables
    current_label = 2
    site_labels = zeros(Int64, numSites(unitcell)) .- 1
    site_labels[index] = 1

    # the list of real labels
    real_labels = zeros(Int64, numSites(unitcell)) .- 1
    function reallabel(l :: Int64) :: Int64
        if real_labels[l] == -1
            return l
        else
            return reallabel(real_labels[l])
        end
    end

    # the list of bonds (organized)
    bond_list = organizedBondsFrom(unitcell)


    # iterate over all sites
    for s in 1:numSites(unitcell)
        # skip the preset one
        if s == index
            continue
        end
        # check neighbors
        for nb in bond_list[s]
            # check if a label is already set
            if site_labels[to(nb)] > 0
                # set the current label
                if site_labels[s] > 0 && reallabel(site_labels[s])==reallabel(site_labels[to(nb)])
                    # do nothing
                elseif site_labels[s] > 0
                    # make new connection
                    real_labels[max(site_labels[s], site_labels[to(nb)])] = min(site_labels[s], site_labels[to(nb)])
                    # label accordingly
                    site_labels[s] = min(site_labels[s], site_labels[to(nb)])
                else
                    # just use the label
                    site_labels[s] = reallabel(site_labels[to(nb)])
                end
            end
        end
        # if still no label, create new one
        if site_labels[s] < 0
            site_labels[s] = current_label
            current_label = current_label + 1
        end
    end

    # replace all labels by their real counterpart
    site_labels = map(s -> reallabel(s), site_labels)

    # remove all labels that are not 1
    sites_to_remove = [i for i in 1:numSites(unitcell) if site_labels[i]!=1]
    if length(sites_to_remove) == 0
        return S[]
    else
        return removeSite!(unitcell, sites_to_remove)
    end
end

# Remove multiple disconnected sites (from site with index) from a lattice (and return them)
function removeDisconnectedSites!(
        lattice :: L,
        index   :: Integer = 1
    ) :: Vector{SL} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B},
        DL,NL,LSL,LBL,SL<:AbstractSite{LSL,DL}, BL<:AbstractBond{LBL,NL}, L<:AbstractLattice{SL,BL,U}}

    #######################################
    # label all sites by Hoshen Kopelman
    #######################################

    # label variables
    current_label = 2
    site_labels = zeros(Int64, numSites(lattice)) .- 1
    site_labels[index] = 1

    # the list of real labels
    real_labels = zeros(Int64, numSites(lattice)) .- 1
    function reallabel(l :: Int64) :: Int64
        if real_labels[l] == -1
            return l
        else
            return reallabel(real_labels[l])
        end
    end

    # the list of bonds (organized)
    bond_list = organizedBondsFrom(lattice)


    # iterate over all sites
    for s in 1:numSites(lattice)
        # skip the preset one
        if s == index
            continue
        end
        # check neighbors
        for nb in bond_list[s]
            # check if a label is already set
            if site_labels[to(nb)] > 0
                # set the current label
                if site_labels[s] > 0 && reallabel(site_labels[s])==reallabel(site_labels[to(nb)])
                    # do nothing
                elseif site_labels[s] > 0
                    # make new connection
                    real_labels[max(site_labels[s], site_labels[to(nb)])] = min(site_labels[s], site_labels[to(nb)])
                    # label accordingly
                    site_labels[s] = min(site_labels[s], site_labels[to(nb)])
                else
                    # just use the label
                    site_labels[s] = reallabel(site_labels[to(nb)])
                end
            end
        end
        # if still no label, create new one
        if site_labels[s] < 0
            site_labels[s] = current_label
            current_label = current_label + 1
        end
    end

    # replace all labels by their real counterpart
    site_labels = map(s -> reallabel(s), site_labels)

    # remove all labels that are not 1
    sites_to_remove = [i for i in 1:numSites(lattice) if site_labels[i]!=1]
    if length(sites_to_remove) == 0
        return S[]
    else
        return removeSite!(lattice, sites_to_remove)
    end
end

# export functions
export removeDisconnectedSites!




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
