################################################################################
#
#   ADDING SINGLE BONDS
#       - to unitcells
#       - to lattices
#
################################################################################

# Add a single bond to a unitcell and return its index in the list
function addBond!(
            unitcell :: U,
            from  :: Integer,
            to    :: Integer,
            label :: L = getDefaultLabel(L),
            wrap  :: NTuple{N,Int64} = NTuple{N,Int64}(zeros(N)),
            returning_bond :: Bool = true
        ) :: Int64 where {L,N,B<:AbstractBond{L,N},S,U<:AbstractUnitcell{S,B}}

    # create a new site
    bond = newBond(B, from, to, label, wrap)

    # obtain the list of bonds in the unitcell
    bond_list = bonds(unitcell)
    # maybe push the returning bond into the list first
    if returning_bond
        # build return wrap
        return_wrap = NTuple{N,Int64}([-w for w in wrap])
        # push to list
        push!(bond_list, newBond(B, to, from, label, return_wrap))
    end
    # push to the list of bonds
    push!(bond_list, bond)
    # set the list of bonds in the unitcell
    bonds!(unitcell, bond_list)

    # return the new bond index
    return length(bond_list)
end


# Add a single bond to a lattice and return its index in the list
function addBond!(
            lattice  :: LA,
            from  :: Integer,
            to    :: Integer,
            label :: L = getDefaultLabel(L),
            wrap  :: NTuple{N,Int64} = NTuple{N,Int64}(zeros(N)),
            returning_bond :: Bool = true
        ) :: Int64 where {L,N,B<:AbstractBond{L,N},S,U,LA<:AbstractLattice{S,B,U}}

    # create a new site
    bond = newBond(B, from, to, label, wrap)

    # obtain the list of bonds in the lattice
    bond_list = bonds(lattice)
    # maybe push the returning bond into the list first
    if returning_bond
        # build return wrap
        return_wrap = NTuple{N,Int64}([-w for w in wrap])
        # push to list
        push!(bond_list, newBond(B, to, from, label, return_wrap))
    end
    # push to the list of bonds
    push!(bond_list, bond)
    # set the list of bonds in the lattice
    bonds!(lattice, bond_list)

    # return the new bond index
    return length(bond_list)
end

# export the function
export addBond!
