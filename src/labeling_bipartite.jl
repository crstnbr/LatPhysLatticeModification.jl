# Function to relabel a lattice to bipartite labeling
function relabelBipartite!(
            lattice :: L,
            label_A :: LS,
            label_B :: LS,
            seed_site_A :: Int64 = 1
        ) where {LS,D,S<:AbstractSite{LS,D},B,U,L<:AbstractLattice{S,B,U}}


    #-------------
    # PREPARATION
    #-------------

    # indicate that all labels are not set yet
    relabeled  = Bool[false for i in 1:numSites(lattice)]
    # save the sublattice as an integer
    sublattice = Int64[-1 for i in 1:numSites(lattice)]

    # get the ordered connection list of the lattice
    bonds_organized = organizedBondsFrom(lattice)

    # set the first label to be sublattice A
    label!(site(lattice, 1), deepcopy(label_A))
    # set the sublattice to be 1
    sublattice[1] = 1
    # indicate that this site has indeed been relabeled
    relabeled[1] = true


    #------------
    # RELABELING
    #------------

    # iterate while there are unlabeled labels left
    while false in relabeled
        # look through all sites
        for i in 1:numSites(lattice)
            # continue if label already set
            if relabeled[i]
                continue
            end
            # check the neighbors for their sublattice
            neighbor_sublattice = -1
            for b in bonds_organized[i]
                if relabeled[to(b)]
                    # check if the sublattice has to be set
                    if neighbor_sublattice < 0
                        neighbor_sublattice = sublattice[to(b)]
                    elseif neighbor_sublattice != sublattice[to(b)]
                        # ERROR - lattice is not bipartite
                        error("The given lattice cannot be labeled to bipartite")
                    end
                end
            end
            # check if neighbors sublattice could be identified
            if neighbor_sublattice > 0
                # set the label depending on the neighbor sublattice
                if neighbor_sublattice == 1
                    # set the sublattice to 2
                    sublattice[i] = 2
                    # set the label B
                    label!(site(lattice, i), deepcopy(label_B))
                    # indicate that the label has been found
                    relabeled[i] = true
                    # set all neighbor sublattices
                    for b in bonds_organized[i]
                        # set the sublattice to 1
                        sublattice[to(b)] = 1
                        # set the label A
                        label!(site(lattice, to(b)), deepcopy(label_A))
                        # indicate that the label has been found
                        relabeled[to(b)] = true
                    end
                elseif neighbor_sublattice == 2
                    # set the sublattice to 1
                    sublattice[i] = 1
                    # set the label A
                    label!(site(lattice, i), deepcopy(label_A))
                    # indicate that the label has been found
                    relabeled[i] = true
                    # set all neighbor sublattices
                    for b in bonds_organized[i]
                        # set the sublattice to 2
                        sublattice[to(b)] = 2
                        # set the label B
                        label!(site(lattice, to(b)), deepcopy(label_B))
                        # indicate that the label has been found
                        relabeled[to(b)] = true
                    end
                else
                    # give an error
                    error("somewhere the neighbor sublattice took the value "*string(neighbor_sublattice))
                end
            end
        end
    end
end

# wrapper using the default labels A and B
function relabelBipartite!(
        lattice :: L,
        seed_site_A :: Int64 = 1
    ) where {LS,D,S<:AbstractSite{LS,D},B,U,L<:AbstractLattice{S,B,U}}

    # call the general function
    return relabelBipartite!(lattice, getDefaultLabelA(LS), getDefaultLabelB(LS), seed_site_A)
end


# export the function
export relabelBipartite!
