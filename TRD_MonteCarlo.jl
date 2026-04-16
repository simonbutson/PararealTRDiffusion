# Monte Carlo Thermal Radiation Diffusion Solvers

module MonteCarlo

using Statistics
using Random

include("TRD_Inputs.jl")
include("TRD_Mesh.jl")
include("TRD_Solvers.jl")

function Updates(params, inputs, n)
(; nx, dt, dx, c, a, cV, rho, sigma_a, f, T, D_centers, D_edges, beta, f_plus, f_minus, P_plus, P_minus, P_c, P_a, matrix_diag) = params


    if inputs["cVmodel"] == "LINEARIZED"
        cV .= (4.0 * a * T[:, n].^3) # Update specific heat for time-Step
    end


    beta .= (4.0 * a * T[:, n].^3) ./ (rho .* cV)  # Update beta for time-Step
    f .= 1 ./(1 .+ c*dt*sigma_a.*beta) # Update fleck factor for time-Step
    
   
    if inputs["LBC"] == "MARSHAK"
        matrix_diag[1]= 1 + c*dt*f[1]*sigma_a[1] + (c*dt/dx[1])*(1.5*D_centers[1]*sigma_a[1]/(1+3*dx[1]*sigma_a[1]/4) + D_edges[2]/(dx[2]))
    elseif inputs["LBC"] == "REFLECT"
        matrix_diag[1]= 1 + c*dt*f[1]*sigma_a[1] + (c*dt/dx[1])*(D_centers[1]/(dx[2]))
    end
    matrix_diag[2:end-1] .= ones(nx-2) + c*dt*f[2:end-1].*sigma_a[2:end-1] + c*dt*(D_edges[1:end-1]./dx[2:end-1] .+ D_edges[2:end]./dx[2:end-1])./dx[2:end-1] 
    
    if inputs["RBC"] == "MARSHAK"
        matrix_diag[end]= 1 + c*dt*f[end]*sigma_a[end] + (c*dt/dx[end])*(1.5*D_centers[end]*sigma_a[end]/(1+3*dx[end]*sigma_a[end]/4) + D_edges[end-1]/(dx[end-1]))
    elseif inputs["RBC"] == "REFLECT"
        matrix_diag[end]= 1 + c*dt*f[end]*sigma_a[end] + (c*dt/dx[end])*(D_centers[end]/(dx[end-1]))
    end
    
    f_plus[1] = (c*dt / dx[1]) * D_edges[1] / dx[1] 
    f_plus[2:end-1] .= (c*dt ./ dx[2:end-1]).* (D_edges[2:end] ./ dx[2:end-1]) 
    f_plus[end] = (c*dt/dx[end])*(D_centers[end] /dx[end]) 
    
    f_minus[1] = (c*dt/dx[1]) * D_centers[1] / dx[1] 
    f_minus[2:end-1] .= (c*dt./dx[2:end-1]) .* (D_edges[1:end-1] ./ dx[2:end-1]) 
    f_minus[end] = (c*dt/dx[end])*D_edges[end] /dx[end] 

    P_plus[1:end-1] .= f_plus[1:end-1] ./ matrix_diag[1:end-1]
    if inputs["RBC"] == "MARSHAK"
        P_plus[end] = (c*dt/dx[end])*(1.5*D_centers[end]*sigma_a[end]) / ((1+3*dx[end]*sigma_a[end]/4)*matrix_diag[end])
    elseif inputs["RBC"] == "REFLECT"
        P_plus[end] = 0.0
    end
    if inputs["LBC"] == "MARSHAK"
        P_minus[1] = (c*dt/dx[1])*(1.5*D_centers[1]*sigma_a[1]) / ((1+3*dx[1]*sigma_a[1]/4)*matrix_diag[1])
    elseif inputs["LBC"] == "REFLECT"
        P_minus[1] = 0.0
    end
    P_minus[2:end] .= f_minus[2:end] ./ matrix_diag[2:end]
    P_c .= ones(nx)./matrix_diag
    P_a .= c*sigma_a.*f.*dt./matrix_diag
end


function Sourcing(params, inputs, n)
    (; nx, dt, dx, a, c, f, sigma_a, T, T_src, E, E_source, E_emitted, E_total, S, N_particles, N_source) = params
    E_emitted[:,n+1] .= c*dt.*f[:].*sigma_a[:].*a.*T[:,n].^4 
    E_source[:,n+1] .= E[:,n] .+ E_emitted[:,n+1] .+ dt*dx.*S[:,n]
    E_source[1,n+1] += (a*c*dt*T_src[1]^4)/(2*dx[1])  # Left Boundary Source
    E_source[end,n+1] += (a*c*dt*T_src[2]^4)/(2*dx[end])  # Right Boundary Source
    E_total[n+1] = sum(E_source[:,n+1])
    #print("Total energy at time ", round((n)*dt, sigdigits=4), ": ", E_total[n+1], "\n")
    

    for i in 1:nx
        N_source[i] = Int(round(N_particles*E_source[i,n+1]/E_total[n+1])+1)
    end
end

function MC(params, inputs, n, rng::AbstractRNG=Random.default_rng())
    (; nx, N_source, E_source, P_plus, P_minus, P_c, P_a, E_rad_tally, E_mat_tally, E_escaped) = params
    E_rad_tally .= 0.0 # Radiation energy tally
    E_mat_tally .= 0.0 # Material energy deposition tally
    E_escaped = 0.0 # Escaped energy tally
    for i in 1:nx
        for _ in 1:N_source[i]
            x = Int(i)
            E_particle = E_source[i,n+1]/N_source[i]
            E_start = E_particle
            while E_particle > 1e-2 * E_start
                u = rand(rng)
                E_rad_tally[x] += E_particle * P_c[x]
                E_mat_tally[x] += E_particle * P_a[x]
                E_particle *= (1 - P_c[x] - P_a[x])
                if u <= P_minus[x] / (P_minus[x] + P_plus[x])
                    x -= 1
                    if x == 0
                        E_escaped += E_particle
                        break
                    end
                else
                    x +=1
                    if x == nx+1
                        E_escaped += E_particle
                        break
                    end
                end
            end
            if x !=0 && x != nx+1
                E_rad_tally[x] += E_particle * P_c[x] / (P_c[x] + P_a[x])
                E_mat_tally[x] += E_particle * P_a[x] / (P_c[x] + P_a[x])
            end
        end
    end
end

function Tally(params, inputs, n)
    (; dx, dt,rho, cV, E_rad_tally, E_mat_tally, E_emitted, E_total, E_escaped, T, E, E_m) = params
    T[:,n+1] = T[:,n] + (E_mat_tally .- E_emitted[:,n+1]) ./ (rho.*cV) # Update Temperature
    E[:,n+1] = E_rad_tally # Update Radiation Energy 
    E_m[:,n+1] = E_m[:,n] + (E_mat_tally .- E_emitted[:,n+1]) # Update Material Energy
    
    #print("Energy check at time ", round((n)*dt, sigdigits=4), ": ", (E_total[n+1] - (sum(E_mat_tally) + sum(E_rad_tally) + E_escaped))/E_total[n+1], "\n")

    # print("Material energy sum: ", sum(params.E_m[:,n+1]), "\n")
    # print("Radiation energy sum: ", sum(params.E[:,n+1]), "\n")
    # print("Temperature sum: ", sum(params.T[:,n+1]), "\n")
end


function MC_Main(params, inputs)
    for n in 1:params.nt
        Updates(params, inputs, n)
        Sourcing(params, inputs, n)
        MC(params, inputs, n)
        Tally(params, inputs, n)
    end
    return params.E, params.E_m, params.T
end

function Smooth(v)
    v_smooth = copy(v)
    n = length(v)
    for i in 2:n-1
        v_smooth[i] = 0.25*v[i-1] + 0.5*v[i] + 0.25*v[i+1]
    end
    return v_smooth
end

function Parareal_MC_Main(params, inputs, threads, epsilon, E_ref, E_m_ref)  
    # Coarse Parameters
    (; nx, nt, dt, dx, x_nodes, x_centers, t, Tm_init, Tr_init, c, a, T_src, rho, cV, sigma_a, beta, f, T, D_centers, D_edges, S, E_m, E, E_source, E_emitted, E_rad_tally, E_mat_tally, E_escaped, E_total, F, f_plus, f_minus, P_plus, P_minus, P_c, P_a, matrix_diag, N_source, N_particles, N_coarse) = params
    E_parallel = (a*c*Tr_init^4)*ones(nx, nt+1,threads)
    E_m_parallel = zeros(nx, nt+1,threads)
    T_parallel = Tm_init*ones(nx, nt+1,threads)
    E_parallel[:,1,:] .= E[:,1]
    E_m_parallel[:,1,:] .= E_m[:,1]
    T_parallel[:,1,:] .= T[:,1]

    Tend = dt * nt
    dt_coarse = Tend / (threads * N_coarse)
    chunk_size = div(nt, threads)
    E_coarse = (a*c*Tr_init^4)*ones(nx, N_coarse*threads+1)
    E_m_coarse = zeros(nx, N_coarse*threads+1)
    T_coarse = Tm_init*ones(nx, N_coarse*threads+1)
    F_coarse = zeros(nx+1, N_coarse*threads+1)
    E_m_coarse[:,1] = E_m[:,1]
    S_coarse = zeros(nx, N_coarse*threads+1)
    for n in 1:N_coarse*threads+1
        S_coarse[:,n] = S[:,n]
    end

    if uppercase(inputs["solver"]) == "HYBRID"
        coarse_params = Mesh.Params(dx, dt_coarse, nx, threads*N_coarse, 1, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, T_coarse, D_centers, D_edges, E_coarse, E_m_coarse, F_coarse, S_coarse, N_coarse)
        E_coarse, E_m_coarse, T_coarse = Solvers.SerialTRTDiffusion(coarse_params, inputs)
    else
        coarse_params = Mesh.ParamsMC(nx, threads*N_coarse, 1, x_nodes, dx, x_centers, t, dt_coarse, Tm_init, Tr_init, c, a, T_src, rho, cV, sigma_a, beta, f, T_coarse, D_centers, D_edges, S_coarse, E_m_coarse, E_coarse, E_source, E_emitted, E_rad_tally, E_mat_tally, E_escaped, E_total, F, f_plus, f_minus, P_plus, P_minus, P_c, P_a, matrix_diag, N_source, N_particles, N_coarse)
    
        E_coarse, E_m_coarse, T_coarse = MC_Main(coarse_params, inputs)
    end

    E_coarse = E_coarse[:, 1:N_coarse:end]
    
    E_m_coarse = E_m_coarse[:, 1:N_coarse:end]
    T_coarse = T_coarse[:, 1:N_coarse:end]
    #F_coarse = F_coarse[:, 1:N_coarse:end]
    #S_coarse = S_coarse[:, 1:N_coarse:end]

    E_coarse_old = copy(E_coarse)
    E_m_coarse_old = copy(E_m_coarse)
    T_coarse_old = copy(T_coarse)

    old_E_error = 1.0
    old_E_m_error = 1.0
    E_spectral_radii = Vector{Float64}()
    E_m_spectral_radii = Vector{Float64}()
    local_E = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_E_m = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_T = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_f = [zeros(nx) for _ in 1:threads]
    local_S = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_F = [zeros(nx+1, chunk_size+1) for _ in 1:threads]    
    local_cV = [zeros(nx) for _ in 1:threads]
    local_params = Vector{Mesh.ParamsMC}(undef, threads)
    
    new_E_error = 2.0 * epsilon 
    new_E_m_error = 2.0 * epsilon
    iterations = 0 
    
    E_prev_iter = copy(E_parallel[:,:,1])
    E_m_prev_iter = copy(E_m_parallel[:,:,1])

    rng = Xoshiro(1234) # Base RNG for reproducibility, can be used to create independent RNGs for each thread if needed

    while (new_E_error > epsilon || new_E_m_error > epsilon) && iterations < threads
        iterations += 1
        
        if iterations > 1
            E_prev_iter .= E_parallel[:,:,iterations-1]
            E_m_prev_iter .= E_m_parallel[:,:,iterations-1]
        end

        Threads.@threads for p in 1:threads
            # Set independent, repeatable seed for this thread/chunk to avoid race conditions
            rng = Xoshiro(1234 + p) 

            local_E[p][:,1] = E_coarse[:,p]
            local_E_m[p][:,1] = E_m_coarse[:,p]
            local_T[p][:,1] = T_coarse[:,p]

            # IMPLEMENTATION OF INCREASING PARTICLES
            # Scale particle count by iteration number to reduce noise in later stages
            current_N_particles = N_particles * iterations
            
            # Update local parameters with new particle count
            local_params[p] = Mesh.ParamsMC(nx, chunk_size, 1, x_nodes, dx, x_centers, t, dt, Tm_init, Tr_init, c, a, T_src, rho, cV, sigma_a, copy(beta), copy(f), local_T[p], D_centers, D_edges, local_S[p], local_E_m[p], local_E[p], copy(E_source), copy(E_emitted), copy(E_rad_tally), copy(E_mat_tally), copy(E_escaped), copy(E_total), local_F[p], copy(f_plus), copy(f_minus), copy(P_plus), copy(P_minus), copy(P_c), copy(P_a), copy(matrix_diag), copy(N_source), current_N_particles, N_coarse)

            for n in 1:chunk_size
                Updates(local_params[p], inputs, n)
                Sourcing(local_params[p], inputs, n)
                MC(local_params[p], inputs, n, rng)
                Tally(local_params[p], inputs, n)
            end

        end
        
        # Combine local results into global arrays
        E_parallel[:, 2:end, iterations] = hcat([local_E[p][:,2:end] for p in 1:threads]...)
        E_m_parallel[:, 2:end, iterations] = hcat([local_E_m[p][:,2:end] for p in 1:threads]...)
        T_parallel[:, 2:end, iterations] = hcat([local_T[p][:,2:end] for p in 1:threads]...)

        # Prediction-Correction
        for n in 2:threads+1
            
            if uppercase(inputs["solver"]) == "HYBRID"
                coarse_params = Mesh.Params(dx, dt_coarse, nx, N_coarse, 1, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, hcat(T_coarse[:,n-1], zeros(nx, N_coarse)), D_centers, D_edges, hcat(E_coarse[:,n-1], zeros(nx, N_coarse)), hcat(E_m_coarse[:,n-1], zeros(nx, N_coarse)), F_coarse, S_coarse, N_coarse)
                E_coarse_new, E_m_coarse_new, T_coarse_new = Solvers.SerialTRTDiffusion(coarse_params, inputs)
            else
                coarse_params = Mesh.ParamsMC(nx, N_coarse, 1, x_nodes, dx, x_centers, t, dt_coarse, Tm_init, Tr_init, c, a, T_src, rho, cV, sigma_a, beta, f, hcat(T_coarse[:,n-1], zeros(nx, N_coarse)), D_centers, D_edges, S_coarse, hcat(E_m_coarse[:,n-1], zeros(nx, N_coarse)), hcat(E[:,n-1], zeros(nx, N_coarse)), E_source, E_emitted, E_rad_tally, E_mat_tally, E_escaped, E_total, F_coarse, f_plus, f_minus, P_plus, P_minus, P_c, P_a, matrix_diag, N_source, N_particles, N_coarse )
                E_coarse_new, E_m_coarse_new, T_coarse_new = MC_Main(coarse_params, inputs)
            end

            # Standard Parareal Update using solution from CURRENT iteration fine solver
            t_idx = Int(1+(n-1)*chunk_size)
            
            E_coarse[:,n]   = E_coarse_new[:,end] .+ Smooth(E_parallel[:, t_idx, iterations] .- E_coarse[:,n])
            E_m_coarse[:,n] = E_m_coarse_new[:,end] .+ Smooth(E_m_parallel[:, t_idx, iterations] .- E_m_coarse[:,n])
            T_coarse[:,n]   = T_coarse_new[:,end] .+ Smooth(T_parallel[:, t_idx, iterations] .- T_coarse[:,n])

            E_coarse_old[:,n] = E_coarse[:,end]
            E_m_coarse_old[:,n] = E_m_coarse[:,end]
            T_coarse_old[:,n] = T_coarse[:,end]

        end

        new_E_error = maximum(abs.(E_parallel[:,:,iterations] - E_prev_iter))
        new_E_m_error = maximum(abs.(E_m_parallel[:,:,iterations] - E_m_prev_iter))
        
        print("Max E value: ", maximum(E_parallel[:,:,iterations]), " at x = ", round(x_centers[argmax(E_parallel[:,:,iterations])[1]], sigdigits=4), " and t = ", round(t[argmax(E_parallel[:,:,iterations])[2]], sigdigits=4), "\n")
        print("Max E_m value: ", maximum(E_m_parallel[:,:,iterations]), " at x = ", round(x_centers[argmax(E_m_parallel[:,:,iterations])[1]], sigdigits=4), " and t = ", round(t[argmax(E_m_parallel[:,:,iterations])[2]], sigdigits=4), "\n")
        print("Max T value: ", maximum(T_parallel[:,:,iterations]), " at x = ", round(x_centers[argmax(T_parallel[:,:,iterations])[1]], sigdigits=4), " and t = ", round(t[argmax(T_parallel[:,:,iterations])[2]], sigdigits=4), "\n")

        print("Iteration ", iterations, ": ϵ_R = ", new_E_error, " with spectral radius: ", new_E_error / old_E_error ," at x = ", round(x_centers[argmax(abs.(E_ref - E_parallel[:,:,iterations]))[1]], sigdigits=4), " and t = ", round(t[argmax(abs.(E_ref - E_parallel[:,:,iterations]))[2]], sigdigits=4), "\n")
        print("Iteration ", iterations, ": ϵ_M = ", new_E_m_error, " with spectral radius: ", new_E_m_error / old_E_m_error , " at x = ", round(x_centers[argmax(abs.(E_m_ref - E_m_parallel[:,:,iterations]))[1]], sigdigits=4), " and t = ", round(t[argmax(abs.(E_m_ref - E_m_parallel[:,:,iterations]))[2]], sigdigits=4), "\n")
        push!(E_spectral_radii, new_E_error / old_E_error)
        push!(E_m_spectral_radii, new_E_m_error / old_E_m_error)
        old_E_error = new_E_error
        old_E_m_error = new_E_m_error
    end
    print("Parallel Scheme Converged in ", iterations, " iterations.\n")
    print("Mean Spectral Radius for E: ", mean(E_spectral_radii[2:end]), "\n")
    print("Mean Spectral Radius for E_m: ", mean(E_m_spectral_radii[2:end]), "\n")
    return E_parallel[:,:,iterations], E_m_parallel[:,:,iterations], T_parallel[:,:,iterations]

 
end



end