# Parareal Thermal Radiation Diffusion Mesh Generation

module Mesh

# mutable struct Params
#     dx::Vector{Float64}
#     dt::Float64
#     nx::Int64
#     nt::Int64
#     ngroups::Int64
#     c::Float64
#     a::Float64
#     T_src::Vector{Float64}
#     T_init::Float64
#     x_nodes::Vector{Float64}
#     x_centers::Vector{Float64}
#     t::Vector{Float64}
#     rho::Vector{Float64}
#     cV::Vector{Float64}
#     sigma_a::Matrix{Float64}
#     beta::Vector{Float64}
#     f::Array{Float64,2}
#     T::Array{Float64,2}
#     D_centers::Matrix{Float64}
#     D_edges::Matrix{Float64}
#     E::Array{Float64,3}
#     E_m::Array{Float64,2}
#     F::Array{Float64,2}
#     S::Array{Float64,3}
# end

mutable struct Params
    dx::Vector{Float64}
    dt::Float64
    nx::Int
    nt::Int
    ngroups::Int
    c::Float64
    a::Float64
    T_src::Vector{Float64}
    Tm_init::Float64
    Tr_init::Float64
    x_nodes::Vector{Float64}
    x_centers::Vector{Float64}
    t::Vector{Float64}
    rho::Vector{Float64}
    cV::Vector{Float64}
    sigma_a
    beta::Vector{Float64}
    f
    T
    D_centers
    D_edges
    E
    E_m
    F
    S
    N_coarse::Int
end

mutable struct ParamsMC
    nx::Int
    nt::Int
    ngroups::Int
    x_nodes::Vector{Float64}
    dx::Vector{Float64}
    x_centers::Vector{Float64}
    t::Vector{Float64}
    dt::Float64
    Tm_init::Float64
    Tr_init::Float64
    c::Float64
    a::Float64
    T_src::Vector{Float64}
    rho::Vector{Float64}
    cV::Vector{Float64}
    sigma_a::Vector{Float64}
    beta::Vector{Float64}
    f::Matrix{Float64}         
    T::Matrix{Float64}
    D_centers::Matrix{Float64}  
    D_edges::Matrix{Float64}    
    S::Matrix{Float64}
    E_m::Matrix{Float64}
    E::Matrix{Float64}
    E_source::Matrix{Float64}
    E_emitted::Matrix{Float64}
    E_rad_tally::Vector{Float64}
    E_mat_tally::Vector{Float64}
    E_escaped::Float64
    E_total::Vector{Float64}
    F::Matrix{Float64}
    f_plus::Vector{Float64}
    f_minus::Vector{Float64}
    P_plus::Vector{Float64}
    P_minus::Vector{Float64}
    P_c::Vector{Float64}
    P_a::Vector{Float64}
    matrix_diag::Vector{Float64}
    N_source::Vector{Int64}
    N_particles::Int
    N_coarse::Int
end



function mesh_generation(inputs)
        """
        Function to generate mesh quantities
        Parameters: 
        inputs: Dict - Dictionary of input parameters
        Returns:
        mesh: MeshStruct - Mesh object
        """

        print("Generating mesh quantities... \n")

        #geometry = uppercase(inputs["GEOMETRY"])
        #meshtype = uppercase(inputs["MESHTYPE"])

        dx_temp = parse(Float64, inputs["dx"])
        dt = parse(Float64, inputs["dt"])
        nx =  parse(Int64, inputs["nx"])
        xsize = dx_temp * nx
        nt = parse(Int64, inputs["nt"])
        ngroups = parse(Int64, inputs["ngroups"])
        dx = fill(dx_temp, nx)
        x_nodes = collect(LinRange(0, xsize, nx+1))
        x_centers = collect(LinRange(dx_temp/2, xsize - dx_temp/2, nx))
        t = collect(LinRange(0, dt*nt, nt+1))
        c = parse(Float64, inputs["c"]) # Speed of light
        a = parse(Float64, inputs["a"]) # Radiation constant
        # Handle T_src as either a scalar or array
        if isa(inputs["T_src"], String)
            T_src = [parse(Float64, inputs["T_src"]), 0.0] # Source Temperature (left boundary, right is vacuum)
        else
            T_src = convert(Vector{Float64}, inputs["T_src"]) # Already an array from input parser
        end
        Tm_init = parse(Float64, inputs["Tm_init"]) # Initial Material Temperature
        Tr_init = parse(Float64, inputs["Tr_init"]) # Initial Radiation Temperature
        rho = parse(Float64, inputs["rho"]) .* ones(nx) # Density
        cV = parse(Float64, inputs["cV"]) .* ones(nx) # Specific Heat
       
        if ngroups > 1
            E = (a*c*Tr_init^4)*ones(nx, nt+1, ngroups) # Radiation Energy
            S = zeros(nx, nt+1, ngroups) # Radiation Source
            sigma_a = ones(nx,ngroups) # Absorption Opacity
            for g in 1:ngroups
                sigma_a[:,g] .= inputs["sigma_a"][g]
                for i in 1:nt+1
                    S[:, i, g] .= region_joiner("1D", collect(inputs["S_regs"]), collect(inputs["S_vals"]), x_nodes, nx)
                end
            end
          
        else
            E = (a*c*Tr_init^4)*ones(nx, nt+1) # Radiation Energy
            S = zeros(nx, nt+1) # Radiation Source
            sigma_a = ones(nx) * parse(Float64, inputs["sigma_a"]) 
            for i in 1:nt+1
                S[:, i] .= region_joiner("1D", collect(inputs["S_regs"]), collect(inputs["S_vals"]), x_nodes, nx)
            end
        end
        beta = ones(nx) # Radiation - Material Energy Coupling Coefficient 
        f = ones(nx, ngroups) # Fleck factor
        T = Tm_init*ones(nx, nt+1) # Temperature
        D_centers = ones(nx,ngroups)./(3*sigma_a) # Diffusion Coefficient - Cell Centered
        D_edges = ones(nx-1,ngroups) # Diffusion Coeffecient - Cell Edges
        E_m = zeros(nx, nt+1) # Material Energy
    
        F = zeros(nx+1, nt+1) # Radiation Flux
        N_coarse = parse(Int64, inputs["N_coarse"]) # Coarse Time Step Factor
       
        #S = zeros(nx, nt+1, ngroups) # Radiation Source
        #S[1:25, :, :] .= 1.0 # Initial Source in first 50 cells

        E_m[:,1] = rho.*cV.*T[:,1] # Initial Material Energy

        D_edges[:,:] = 2*dx[1:end-1].*(D_centers[1:end-1,:].*D_centers[2:end,:])./(dx[2:end].*D_centers[1:end-1,:] + dx[1:end-1].*D_centers[2:end,:])


        if uppercase(inputs["solver"]) == "DETERMINISTIC"
            params = Params(dx, dt, nx, nt, ngroups, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, T, D_centers, D_edges, E, E_m, F, S, N_coarse)
        elseif uppercase(inputs["solver"]) == "MC" || uppercase(inputs["solver"]) == "HYBRID"
            N_particles = parse(Int64, inputs["N_particles"]) # Number of MC Particles
            E_source = zeros(nx, nt+1) # Radiation Energy Source
            E_emitted = zeros(nx, nt+1) # Emitted Radiation Energy
            E_rad_tally = zeros(nx) # Radiation Energy Tally
            E_mat_tally = zeros(nx) # Material Energy Deposition Tally
            E_escaped = 0.0 # Escaped Energy Tally
            E_total = zeros(nt+1) # Total Energy in System
            f_plus = zeros(nx) # Probability of moving right
            f_minus = zeros(nx) # Probability of moving left
            P_plus = zeros(nx) # Probability of moving right normalized
            P_minus = zeros(nx) # Probability of moving left normalized
            P_c = zeros(nx) # Probability of being colliding without absorbing normalized
            P_a = zeros(nx) # Probability of being absorbed normalized
            matrix_diag = zeros(nx) # Diagonal of matrix for implicit solve
            N_source = zeros(Int64, nx) # Number of source particles per cell

            params = ParamsMC(nx, nt, ngroups, x_nodes, dx, x_centers, t, dt, Tm_init, Tr_init, c, a, T_src, rho, cV, sigma_a, beta, f, T, D_centers, D_edges, S, E_m, E, E_source, E_emitted, E_rad_tally, E_mat_tally, E_escaped, E_total, F, f_plus, f_minus, P_plus, P_minus, P_c, P_a, matrix_diag, N_source, N_particles, N_coarse)
        end

        return params
    end


    function region_joiner(geometry, regions, values, nodes, Ncells)
        """ Function to join regions and values for mesh generation
            Parameters:
            geometry: String - Geometry of the problem
            regions: Array - Array of regions
            values: Array - Array of values
            nodes: Array - Array of nodal points
            Ncells: Int - Number of cells
            Returns:
            joined_array: Array - Array of joined values
        """

        joined_array = zeros(Float64, Ncells)

        region_index = 1
        if geometry == "1D"
            if length(regions) <= 1
                joined_array = fill(parse(Float64, values[1]), Ncells)
            else
                for i in eachindex(joined_array)
                    if regions[region_index] >= nodes[i+1]
                        joined_array[i] = values[region_index]
                    else
                        region_index += 1
                        joined_array[i] = values[region_index]
                    end
                end
            end
        elseif geometry == "2D"
            if length(regions[1]) <= 1
                joined_array = fill(parse(Float64, values[1]), Ncells)
            else
                for region_index in eachindex(regions[1])
                        xstart, xend = regions[1][region_index][1]
                        ystart, yend = regions[1][region_index][2]
                        xstartindex = findlast(x -> x <= xstart, nodes[1])
                        xendindex = findlast(x -> x <= xend, nodes[1])-1
                        ystartindex = findlast(y -> y <= ystart, nodes[2])
                        yendindex = findlast(y -> y <= yend, nodes[2])-1
                        joined_array[xstartindex:xendindex, ystartindex:yendindex] .= values[region_index]
                end
            end
        end
        return joined_array
    end

    function surface_definer(geometry, regions, values, perimeter, nodes, Ncells)
        """ Function to define the surface temperature
            Parameters:
            geometry: String - Geometry of the problem
            regions: Array - Array of regions
            values: Array - Array of values
            perimeter: Array - Array of perimeter values
            nodes: Array - Array of nodal points
            Ncells: Int - Number of cells
            Returns:
            surface: Array - Array of surface values
        """
        if geometry == "1D"
            surface = (values[1], values[2])
        elseif geometry == "2D"
            surface = [[], [], [], []]
            surfacevals = [[], [], [], []]
            # Format values for region joiner function
            for i in 1:4
                if length(values[1][i]) <= 1
                    surfacevals[i] = [string(values[1][i])]
                else
                    surfacevals[i] = collect(values[1][i])
                end
            end
            # Bottom Boundary
            surface[1] = region_joiner("1D", collect(regions[1][1]), surfacevals[1], nodes[1], Ncells[1])
            # Top Boundary
            surface[2] = region_joiner("1D", collect(regions[1][2]), surfacevals[2], nodes[1], Ncells[1]) 
            # Left Boundary
            surface[3] = region_joiner("1D", collect(regions[1][3]), surfacevals[3], nodes[2], Ncells[2])
            # Right Boundary
            surface[4] = region_joiner("1D", collect(regions[1][4]), surfacevals[4], nodes[2], Ncells[2])
        end

        return surface 
    end


end