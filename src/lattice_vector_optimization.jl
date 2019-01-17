
# Function to optimize lattice vectors
function optimizeLatticeVectors2D!(unitcell::Unitcell; verbose::Bool=false)
    # init the loop by assuming one refolded a vector
    refolded_one_vector = true
    # as long as one vector was refolded, repeat this loop
    while refolded_one_vector
        # init inside the loop by assuming one did not refold a vector
        refolded_one_vector = false
        # try a1
        if dot(unitcell.lattice_vectors[1],unitcell.lattice_vectors[1]) > dot(unitcell.lattice_vectors[1].-unitcell.lattice_vectors[2], unitcell.lattice_vectors[1].-unitcell.lattice_vectors[2])
            # print
            if verbose
                println("a1 --> a1 - a2")
            end
            # change the lattice vector
            unitcell.lattice_vectors[1] = unitcell.lattice_vectors[1] .- unitcell.lattice_vectors[2]
            # change all connections
            for c in unitcell.connections
                if c[4][1] != 0
                    c[4] = (c[4][1], c[4][2]+c[4][1])
                end
            end
            # did refold
            refolded_one_vector = true
        end
        if dot(unitcell.lattice_vectors[1],unitcell.lattice_vectors[1]) > dot(unitcell.lattice_vectors[1].+unitcell.lattice_vectors[2], unitcell.lattice_vectors[1].+unitcell.lattice_vectors[2])
            # print
            if verbose
                println("a1 --> a1 + a2")
            end
            # change the lattice vector
            unitcell.lattice_vectors[1] = unitcell.lattice_vectors[1] .+ unitcell.lattice_vectors[2]
            # change all connections
            for c in unitcell.connections
                if c[4][1] != 0
                    c[4] = (c[4][1], c[4][2]-c[4][1])
                end
            end
            # did refold
            refolded_one_vector = true
        end
        # try a2
        if dot(unitcell.lattice_vectors[2],unitcell.lattice_vectors[2]) > dot(unitcell.lattice_vectors[2].-unitcell.lattice_vectors[1], unitcell.lattice_vectors[2].-unitcell.lattice_vectors[1])
            # print
            if verbose
                println("a2 --> a2 - a1")
            end
            # change the lattice vector
            unitcell.lattice_vectors[2] = unitcell.lattice_vectors[2] .- unitcell.lattice_vectors[1]
            # change all connections
            for c in unitcell.connections
                if c[4][2] != 0
                    c[4] = (c[4][1]+c[4][2], c[4][2])
                end
            end
            # did refold
            refolded_one_vector = true
        end
        if dot(unitcell.lattice_vectors[2],unitcell.lattice_vectors[2]) > dot(unitcell.lattice_vectors[2].+unitcell.lattice_vectors[1], unitcell.lattice_vectors[2].+unitcell.lattice_vectors[1])
            # print
            if verbose
                println("a2 --> a2 + a1")
            end
            # change the lattice vector
            unitcell.lattice_vectors[2] = unitcell.lattice_vectors[2] .+ unitcell.lattice_vectors[1]
            # change all connections
            for c in unitcell.connections
                if c[4][2] != 0
                    c[4] = (c[4][1]-c[4][2], c[4][2])
                end
            end
            # did refold
            refolded_one_vector = true
        end
    end
end

# Function to optimize lattice vectors
function optimizeLatticeVectors3D!(unitcell::Unitcell; verbose::Bool=false)
    # init the loop by assuming one refolded a vector
    refolded_one_vector = true
    # as long as one vector was refolded, repeat this loop
    while refolded_one_vector
        # init inside the loop by assuming one did not refold a vector
        refolded_one_vector = false
        # try all combinations of the remaining two lattice vectors
        for comb in Array{Int64,1}[[1,0], [-1,0], [0,1], [0,-1], [1,1], [1,-1], [-1,1], [-1,-1]]
            # try a1
            a1p = unitcell.lattice_vectors[1] .+ comb[1].*unitcell.lattice_vectors[2] .+ comb[2].*unitcell.lattice_vectors[3]
            # check if this is closer
            if dot(unitcell.lattice_vectors[1],unitcell.lattice_vectors[1]) > dot(a1p, a1p)
                # print
                if verbose
                    println("a1 --> a1 + $(comb[1])*a2 + $(comb[2])*a3")
                end
                # change the lattice vector
                unitcell.lattice_vectors[1] = a1p
                # change all connections
                for c in unitcell.connections
                    c[4] = (c[4][1], c[4][2]-c[4][1]*comb[1], c[4][3]-c[4][1]*comb[2])
                end
                # did refold
                refolded_one_vector = true
            end
            # try a2
            a2p = unitcell.lattice_vectors[2] .+ comb[1].*unitcell.lattice_vectors[1] .+ comb[2].*unitcell.lattice_vectors[3]
            # check if this is closer
            if dot(unitcell.lattice_vectors[2],unitcell.lattice_vectors[2]) > dot(a2p, a2p)
                # print
                if verbose
                    println("a2 --> a2 + $(comb[1])*a1 + $(comb[2])*a3")
                end
                # change the lattice vector
                unitcell.lattice_vectors[2] = a2p
                # change all connections
                for c in unitcell.connections
                    c[4] = (c[4][1]-c[4][2]*comb[1], c[4][2], c[4][3]-c[4][2]*comb[2])
                end
                # did refold
                refolded_one_vector = true
            end
            # try a3
            a3p = unitcell.lattice_vectors[3] .+ comb[1].*unitcell.lattice_vectors[1] .+ comb[2].*unitcell.lattice_vectors[2]
            # check if this is closer
            if dot(unitcell.lattice_vectors[3],unitcell.lattice_vectors[3]) > dot(a3p, a3p)
                # print
                if verbose
                    println("a3 --> a3 + $(comb[1])*a1 + $(comb[2])*a2")
                end
                # change the lattice vector
                unitcell.lattice_vectors[3] = a3p
                # change all connections
                for c in unitcell.connections
                    c[4] = (c[4][1]-c[4][3]*comb[1], c[4][2]-c[4][3]*comb[2], c[4][3])
                end
                # did refold
                refolded_one_vector = true
            end
        end
    end
end
