# Thermal Radiation Diffusion Solvers

module Solvers

using LinearAlgebra
using Statistics
using Plots

include("TRD_Inputs.jl")
include("TRD_Mesh.jl")

function SerialTRTDiffusion(params, inputs)
    # Extract fields from params struct    
    (; dx, dt, nx, nt, c, a, T_src, rho, cV, sigma_a, beta, f, T, D_centers, D_edges, E_m, E, F, S) = params

    #print("CFL Condition: ", (maximum(D_edges)*c)*dt/minimum(dx)^2, "\n")
    
    dl = zeros(nx-1) # Lower diagonal of matrix A
    d = zeros(nx) # Main diagonal of matrix A
    du = zeros(nx-1) # Upper diagonal of matrix A
    b = zeros(nx) # Source terms
    
    for n in 1:nt
        #print("Time-Step: ", n, "\n")
        F[1,n] =  (a*c*dt/4)*(T_src[1])^4 # Left Temperature Source Boundary flux
        F[end,n] =  (a*c*dt/4)*(T_src[2])^4 # Right Temperature Source Boundary flux
        if uppercase(inputs["cVmodel"]) == "LINEARIZED"
            cV = 4*a*T[:,n].^3 # Update specific heat for time-Step
        end
        beta = (4 * a * T[:, n].^3)./(rho.*cV)  # Update beta for time-Step
        f = 1 ./(1 .+ c*dt*sigma_a.*beta) # Update fleck factor for time-Step

        # Left Boundary Condition
        if inputs["LBC"] == "MARSHAK"
            d[1] = 1 + c*dt*f[1]*sigma_a[1] + (c*dt/dx[1])*1.5*D_centers[1]*sigma_a[1]/(1 + 3*dx[1]*sigma_a[1]/4) + (c*dt/dx[1])*(D_edges[1]/dx[2])
            du[1] = -(c*dt/dx[1])*(D_edges[1]/dx[2])
            b[1] = c*dt*f[1]*sigma_a[1]*a*T[1,n]^4 + E[1,n] + 4*(1.5*D_centers[1]*sigma_a[1]/(1 + 3*dx[1]*sigma_a[1]/4))*F[1,n]/(c*dx[1]) + dt*dx[1]*S[1,n]
        elseif inputs["LBC"] == "REFLECT"
            # Zero flux condition: -D * dE/dx = 0 at boundary
            # Use one-sided difference: dE/dx ≈ (E[2] - E[1])/dx[1] = 0
            # This gives E[1] = E[2], but implemented via the finite difference equation
            d[1] = 1 + c*dt*f[1]*sigma_a[1] + c*dt*D_edges[1]/(dx[1]*dx[2])
            du[1] = -c*dt*D_edges[1]/(dx[1]*dx[2])
            b[1] = c*dt*f[1]*sigma_a[1]*a*T[1,n]^4 + E[1,n] + dt*dx[1]*S[1,n]
        end

        for j in 2:nx-1
            dl[j-1] = -c*dt*D_edges[j-1]/(dx[j-1]*dx[j])
            d[j] = 1 + c*dt*f[j]*sigma_a[j] + c*dt*(D_edges[j-1]/dx[j] + D_edges[j]/dx[j+1])/dx[j]
            du[j] = -c*dt*D_edges[j]/(dx[j]*dx[j+1])
            b[j] = c*dt*f[j]*sigma_a[j]*a*T[j,n]^4 + E[j,n] + dt*dx[j]*S[j,n]
        end
        # Right Boundary Condition
        if inputs["RBC"] == "MARSHAK"
            dl[nx-1] = -(c*dt/dx[nx-1])*(D_edges[nx-1]/dx[nx])
            d[nx] = 1 + c*dt*f[nx]*sigma_a[nx] + (c*dt/dx[nx])*(1.5*D_centers[nx]*sigma_a[nx]/(1 + 3*dx[nx]*sigma_a[nx]/4)) + (c*dt/dx[nx])*(D_edges[nx-1]/dx[nx])
            b[nx] = c*dt*f[nx]*sigma_a[nx]*a*T[nx,n]^4 + E[nx,n] + 4*(1.5*D_centers[nx]*sigma_a[nx]/(1 + 3*dx[nx]*sigma_a[nx]/4))*F[nx+1,n]/(c*dx[nx]) + dt*dx[nx]*S[nx,n] # Vacuum Right Boundary
        elseif inputs["RBC"] == "REFLECT"
            dl[nx-1] = -c*dt*D_edges[nx-1]/(dx[nx-1]*dx[nx])
            d[nx] = 1 + c*dt*f[nx]*sigma_a[nx] + c*dt*D_edges[nx-1]/(dx[nx-1]*dx[nx])
            b[nx] = c*dt*f[nx]*sigma_a[nx]*a*T[nx,n]^4 + E[nx,n] + dt*dx[nx]*S[nx,n]
        end

        TD = Tridiagonal(dl, d, du)

        E[:,n+1] = TD \ b # Solve for new radiation energy

        E_m[:,n+1] = E_m[:,n] + dt*c*f.*sigma_a.*(E[:,n+1] - a.*T[:,n].^4) # Material energy update

        T[:,n+1] = T[:,n] + (E_m[:,n+1]-E_m[:,n]) ./ (rho.*cV) # Update Material Temperature

        #T[:,n+1] = (T[:,n] + (c*dt*f.*sigma_a./(rho.*cV)).*(E[:,n+1] + 3*a*T[:,n].^4)) ./ (1 .+ 4*c*dt*f.*sigma_a.*a.*(T[:,n].^3)./(rho.*cV)) # Update Material Temperature

        #E_m[:,n+1] = rho.*cV.*T[:,n+1] # Update Material Energy
    end
    return E, E_m, T
end

function SerialTRTDiffusion_MG(params, inputs)
    # Extract fields from params struct for multi-group version
    (; dx, dt, nx, nt, ngroups, c, a, T_src, rho, cV, sigma_a, beta, f, T, D_centers, D_edges, E_m, E, F, S) = params
    
    for n in 1:nt
        #print("Time-Step: ", n, "\n")
        F[1,n] = (a*c*dt/4)*T_src[1].^4 # Left Temperature Source Boundary flux
        F[end,n] = (a*c*dt/4)*T_src[2].^4 # Right Vacuum Boundary flux
        if uppercase(inputs["cVmodel"]) == "LINEARIZED"
            cV = 4*a*T[:,n].^3 # Update specific heat for time-Step
        end
        beta = (4 * a * T[:, n].^3)./(rho.*cV)  # Update beta for time-Step
    

        for g in 1:ngroups
            dl = zeros(nx-1) # Lower diagonal of matrix A
            d = zeros(nx) # Main diagonal of matrix A
            du = zeros(nx-1) # Upper diagonal of matrix A
            b = zeros(nx) # Source terms
            f[:,g] = 1 ./(1 .+ c*dt*sigma_a[:,g].*beta) # Update fleck factor for all groups
            
            # Left Boundary Condition
            if uppercase(inputs["LBC"]) == "MARSHAK"
                d[1] = 1 + c*dt*f[1,g]*sigma_a[1,g] + (c*dt/dx[1])*1.5*D_centers[1,g]*sigma_a[1,g]/(1 + 3*dx[1]*sigma_a[1,g]/4) + (c*dt/dx[1])*(D_edges[1,g]/dx[2])
                du[1] = -(c*dt/dx[1])*(D_edges[1,g]/dx[2])
                b[1] = c*dt*f[1,g]*sigma_a[1,g]*a*T[1,n]^4/ngroups + E[1,n,g] + 4*(1.5*D_centers[1,g]*sigma_a[1,g]/(1 + 3*dx[1]*sigma_a[1,g]/4))*F[1,n]/(c*dx[1]*ngroups) + dt*dx[1]*S[1,n,g]/(ngroups)
            elseif uppercase(inputs["LBC"]) == "REFLECT"
                # Zero flux condition for each group
                d[1] = 1 + c*dt*f[1,g]*sigma_a[1,g] + c*dt*D_edges[1,g]/(dx[1]*dx[2])
                du[1] = -c*dt*D_edges[1,g]/(dx[1]*dx[2])
                b[1] = c*dt*f[1,g]*sigma_a[1,g]*a*T[1,n]^4/ngroups + E[1,n,g] + dt*dx[1]*S[1,n,g]/(ngroups)
            end

            for j in 2:nx-1
                dl[j-1] = -c*dt*D_edges[j-1,g]/(dx[j-1]*dx[j])
                d[j] = 1 + c*dt*f[j,g]*sigma_a[j,g] + c*dt*(D_edges[j-1,g]/dx[j] + D_edges[j,g]/dx[j+1])/dx[j]
                du[j] = -c*dt*D_edges[j,g]/(dx[j]*dx[j+1])
                b[j] = (c*dt*f[j,g]*sigma_a[j,g]*a*T[j,n]^4)/ngroups + E[j,n,g] + dt*dx[j]*S[j,n,g]/(ngroups)
            end
            
            # Right Boundary Condition
            if uppercase(inputs["RBC"]) == "MARSHAK"
                dl[nx-1] = -(c*dt/dx[nx-1])*(D_edges[nx-1,g]/dx[nx])
                d[nx] = 1 + c*dt*f[nx,g]*sigma_a[nx,g] + (c*dt/dx[nx])*(1.5*D_centers[nx,g]*sigma_a[nx,g]/(1 + 3*dx[nx]*sigma_a[nx,g]/4)) + (c*dt/dx[nx])*(D_edges[nx-1,g]/dx[nx])
                b[nx] = c*dt*f[nx,g]*sigma_a[nx,g]*a*T[nx,n]^4/ngroups + E[nx,n,g] + 4*(1.5*D_centers[nx,g]*sigma_a[nx,g]/(1 + 3*dx[nx]*sigma_a[nx,g]/4))*F[nx+1,n]/(c*dx[nx]*ngroups) + dt*dx[nx]*S[nx,n,g]/(ngroups)
            elseif uppercase(inputs["RBC"]) == "REFLECT"
                # Zero flux condition for each group
                dl[nx-1] = -c*dt*D_edges[nx-1,g]/(dx[nx-1]*dx[nx])
                d[nx] = 1 + c*dt*f[nx,g]*sigma_a[nx,g] + c*dt*D_edges[nx-1,g]/(dx[nx-1]*dx[nx])
                b[nx] = c*dt*f[nx,g]*sigma_a[nx,g]*a*T[nx,n]^4/ngroups + E[nx,n,g] + dt*dx[nx]*S[nx,n,g]/(ngroups)
            end

            TD = Tridiagonal(dl, d, du)
            E[:,n+1,g] = TD \ b # Solve for new radiation energy for group g
        end
        # Update material energy and temperature (summing over all groups)
        E_m[:,n+1] = E_m[:,n] + dt*c*sum(f .* sigma_a .* (E[:,n+1,:] .- (a.*T[:,n].^4)/ngroups), dims=2)
        #E_m[:,n+1] = E_m[:,n] + dt*c*sum(f .* sigma_a .* E[:,n+1,:],dims=2) .- dt*c* (a.*T[:,n].^4)
        #E_m[:,n+1] = E_m[:,n] + dt*c*(f[:,1].*sigma_a[:,1] .* (E[:, n+1, 1] - 0.5*a.*T[:,n].^4)) + dt*c*(f[:,2].*sigma_a[:,2] .* (E[:, n+1, 2] - 0.5*a.*T[:,n].^4))
        
        T[:,n+1] = T[:,n] + (E_m[:,n+1]-E_m[:,n]) ./ (rho.*cV) # Update Material Temperature
    end
    return E, E_m, T
end


function PararealTRTDiffusion(params, inputs, threads, epsilon, E_ref, E_m_ref)
     # Extract fields from params struct
    (; dx, dt, nx, nt, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, T, D_centers, D_edges, E_m, E, F, S, N_coarse) = params
    
    Tend = dt * nt
    dt_coarse = Tend / (threads * N_coarse)
    chunk_size = div(nt, threads)
    E_coarse = (a*c*Tr_init^4)*ones(nx, N_coarse*threads+1)
    E_m_coarse = zeros(nx, N_coarse*threads+1)
    T_coarse = Tm_init*ones(nx, N_coarse*threads+1)
    F_coarse = zeros(nx+1, N_coarse*threads+1)
    E_m_coarse[:,1] .= E_m[:,1]
    S_coarse = zeros(nx, N_coarse*threads+1)
    for n in 1:N_coarse*threads+1
        S_coarse[:,n] .= S[:,1]
    end
  

    #print("Initial temperature check: ", T_init, "\n")
    #print("Initial radiation energy check: ", E[:,1], "\n")

    coarse_params = Mesh.Params(dx, dt_coarse, nx, threads*N_coarse, 1, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, T_coarse, D_centers, D_edges, E_coarse, E_m_coarse, F_coarse, S_coarse, N_coarse)

    E_parallel = (a*c*Tr_init^4)*ones(nx, nt+1)
    E_m_parallel = zeros(nx, nt+1)
    T_parallel = Tm_init*ones(nx, nt+1)
    E_parallel[:,1] .= E[:,1]
    E_m_parallel[:,1] .= E_m[:,1]
    T_parallel[:,1] .= T[:,1]
    # Step 1: Coarse solution
    E_coarse, E_m_coarse, T_coarse = SerialTRTDiffusion(coarse_params, inputs)
    E_coarse = E_coarse[:, 1:N_coarse:end] # Downsample coarse solution arrays
    E_m_coarse = E_m_coarse[:, 1:N_coarse:end]
    T_coarse = T_coarse[:, 1:N_coarse:end]
    #F_coarse = F_coarse[:, 1:N_coarse:end]
    #S_coarse = S_coarse[:, 1:N_coarse:end]
    E_coarse_old = copy(E_coarse)
    E_m_coarse_old = copy(E_m_coarse)
    T_coarse_old = copy(T_coarse)
    E_old = zeros(size(E_parallel))
    E_m_old = ones(size(E_m_parallel))
    old_E_error = 1.0
    old_E_m_error = 1.0
    E_spectral_radii = Vector{Float64}()
    E_m_spectral_radii = Vector{Float64}()
    local_E = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_E_m = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_T = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_f = [zeros(nx) for _ in 1:threads]
    local_F = [zeros(nx+1, chunk_size+1) for _ in 1:threads]
    local_cV = [zeros(nx) for _ in 1:threads]
    iterations = 0
    while maximum(abs.(E_ref - E_parallel)) > epsilon || maximum(abs.(E_m_ref - E_m_parallel)) > epsilon
        #E_old .= copy(E_parallel)
        #E_m_old .= copy(E_m_parallel)

        Threads.@threads for p in 1:threads
            local_E[p][:,1] .= E_coarse[:,p]
            local_E_m[p][:,1] .= E_m_coarse[:,p]
            local_T[p][:,1] .= T_coarse[:,p]

            #print("Thread ", Threads.threadid(), " starting chunk ", p, "\n")

            dl = zeros(nx-1)
            d = zeros(nx)
            du = zeros(nx-1)
            b = zeros(nx) # Source terms

            #print("Thread ", Threads.threadid(), " working on chunk ", p, "\n")

            for n in 1:chunk_size
                # Matrix assembly for TRT-Diffusion equations
                 local_F[p][1,n] = (a*c*dt/4)*(T_src[1])^4
                 local_F[p][end,n] = (a*c*dt/4)*(T_src[2])^4

                if uppercase(inputs["cVmodel"]) == "LINEARIZED"
                    local_cV[p] = 4*a*local_T[p][:, n].^3 # Update specific heat for time-step
                else
                    local_cV[p] = cV # Use original cV if not linearized
                end
                local_beta = (4 * a * local_T[p][:, n].^3)./(rho.*local_cV[p])
                local_f[p] = 1 ./(1 .+ c*dt*sigma_a.*local_beta)
                # Left boundary
                if inputs["LBC"] == "MARSHAK"
                    d[1] = 1 + c*dt*local_f[p][1]*sigma_a[1] + (c*dt/dx[1])*1.5*D_centers[1]*sigma_a[1]/(1 + 3*dx[1]*sigma_a[1]/4) + (c*dt/dx[1])*(D_edges[1]/dx[2])
                    du[1] = -(c*dt/dx[1])*(D_edges[1]/dx[2])
                    b[1] = c*dt*local_f[p][1]*sigma_a[1]*a*local_T[p][1,n]^4 + local_E[p][1,n] + 4*(1.5*D_centers[1]*sigma_a[1]/(1 + 3*dx[1]*sigma_a[1]/4))*local_F[p][1,n]/(c*dx[1]) + dt*dx[1]*S[1,n]
                elseif inputs["LBC"] == "REFLECT"
                    # Zero flux condition
                    d[1] = 1 + c*dt*local_f[p][1]*sigma_a[1] + c*dt*D_edges[1]/(dx[1]*dx[2])
                    du[1] = -c*dt*D_edges[1]/(dx[1]*dx[2])
                    b[1] = c*dt*local_f[p][1]*sigma_a[1]*a*local_T[p][1,n]^4 + local_E[p][1,n] + dt*dx[1]*S[1,n]
                end
                for j in 2:nx-1
                    dl[j-1] = -c*dt*D_edges[j-1]/(dx[j-1]*dx[j])
                    d[j] = 1 + c*dt*local_f[p][j]*sigma_a[j] + c*dt*(D_edges[j-1]/dx[j] + D_edges[j]/dx[j+1])/dx[j]
                    du[j] = -c*dt*D_edges[j]/(dx[j]*dx[j+1])
                    b[j] = c*dt*local_f[p][j]*sigma_a[j]*a*local_T[p][j,n]^4 + local_E[p][j,n] + dt*dx[j]*S[j,n]
                end
                # Right boundary
                if inputs["RBC"] == "MARSHAK"
                    dl[nx-1] = -(c*dt/dx[nx-1])*(D_edges[nx-1]/dx[nx])
                    d[nx] = 1 + c*dt*local_f[p][nx]*sigma_a[nx] + (c*dt/dx[nx])*(1.5*D_centers[nx]*sigma_a[nx]/(1 + 3*dx[nx]*sigma_a[nx]/4)) + (c*dt/dx[nx])*(D_edges[nx-1]/dx[nx])
                    b[nx] = c*dt*local_f[p][nx]*sigma_a[nx]*a*local_T[p][nx,n]^4 + local_E[p][nx,n] + 4*(1.5*D_centers[nx]*sigma_a[nx]/(1 + 3*dx[nx]*sigma_a[nx]/4))*local_F[p][nx+1,n]/(c*dx[nx]) + dt*dx[nx]*S[nx,n]
                elseif inputs["RBC"] == "REFLECT"
                    # Zero flux condition
                    dl[nx-1] = -c*dt*D_edges[nx-1]/(dx[nx-1]*dx[nx])
                    d[nx] = 1 + c*dt*local_f[p][nx]*sigma_a[nx] + c*dt*D_edges[nx-1]/(dx[nx-1]*dx[nx])
                    b[nx] = c*dt*local_f[p][nx]*sigma_a[nx]*a*local_T[p][nx,n]^4 + local_E[p][nx,n] + dt*dx[nx]*S[nx,n]
                end

                TD = Tridiagonal(dl, d, du)
                local_E[p][:,n+1] = TD \ b
                local_E_m[p][:,n+1] = local_E_m[p][:,n] + dt*c*local_f[p].*sigma_a.*(local_E[p][:,n+1] - a.*local_T[p][:,n].^4)
                local_T[p][:,n+1] = local_T[p][:,n] + (local_E_m[p][:,n+1]-local_E_m[p][:,n]) ./ (rho.*local_cV[p])

                #local_T[p][:,n+1] = (local_T[p][:,n] + (c*dt*f.*sigma_a./(rho.*cV)).*(local_E[p][:,n+1] + 3*a*local_T[p][:,n].^4)) ./ (1 .+ 4*c*dt*f.*sigma_a.*a.*(local_T[p][:,n].^3)./(rho.*cV)) # Update Material Temperature

                #local_E_m[p][:,n+1] = rho.*cV.*local_T[p][:,n+1] # Update Material Energy
            end
        end
        # Step 2b: Combine local results into global arrays
        E_parallel[:,2:end] = hcat([local_E[p][:,2:end] for p in 1:threads]...)
        E_m_parallel[:,2:end] = hcat([local_E_m[p][:,2:end] for p in 1:threads]...)
        T_parallel[:,2:end] = hcat([local_T[p][:,2:end] for p in 1:threads]...)
        #print(T_parallel[:,end].^4 /dx[1], "\n")
        #sleep(10)
        # Step 2c: Perform prediction and correction step to update coarse solution with parallel fine solution
        for n in 2:threads+1
            #coarse_params = (nx, N_coarse, dx, Tend/(threads*N_coarse), c, a, T_src, rho, cV, sigma_a, beta, f, hcat(T_coarse[:,n-1], zeros(nx, N_coarse)), D_centers, D_edges, hcat(E_m_coarse[:,n-1], zeros(nx, N_coarse)), hcat(E_coarse[:,n-1], zeros(nx, N_coarse)), F_coarse)
            coarse_params = Mesh.Params(dx, dt_coarse, nx, N_coarse, 1, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, hcat(T_coarse[:,n-1], zeros(nx, N_coarse)), D_centers, D_edges, hcat(E_coarse[:,n-1], zeros(nx, N_coarse)), hcat(E_m_coarse[:,n-1], zeros(nx, N_coarse)), F_coarse, S_coarse, N_coarse)
    
            E_coarse_new, E_m_coarse_new, T_coarse_new = SerialTRTDiffusion(coarse_params, inputs)

            E_coarse[:,n] = E_coarse_new[:,end] .+ (E_parallel[:, Int(1+(n-1)*chunk_size)] .- E_coarse_old[:,n])
            E_m_coarse[:,n] = E_m_coarse_new[:,end] .+ (E_m_parallel[:, Int(1+(n-1)*chunk_size)] .- E_m_coarse_old[:,n])
            T_coarse[:,n] = T_coarse_new[:,end] .+ (T_parallel[:, Int(1+(n-1)*chunk_size)] .- T_coarse_old[:,n])

            E_coarse_old[:,n] .= E_coarse_new[:,end]
            E_m_coarse_old[:,n] .= E_m_coarse_new[:,end]
            T_coarse_old[:,n] .= T_coarse_new[:,end]
        end

        #E_coarse .= E_coarse .+ (E_parallel[:,1:Int(chunk_size/N_coarse):end] .- E_coarse)
        #E_m_coarse .= E_m_parallel[:,1:Int(chunk_size/N_coarse):end]
        #T_coarse .= T_parallel[:,1:Int(chunk_size/N_coarse):end]
        iterations += 1
        new_E_error = maximum(abs.(E_ref - E_parallel))
        new_E_m_error = maximum(abs.(E_m_ref - E_m_parallel))
        print("Iteration ", iterations, ": ϵ_R = ", new_E_error, " with spectral radius: ", new_E_error / old_E_error ," at x = ", round(x_centers[argmax(abs.(E_ref - E_parallel))[1]], sigdigits=4), " and t = ", round(t[argmax(abs.(E_ref - E_parallel))[2]], sigdigits=4), "\n")
        print("Iteration ", iterations, ": ϵ_M = ", new_E_m_error, " with spectral radius: ", new_E_m_error / old_E_m_error , " at x = ", round(x_centers[argmax(abs.(E_m_ref - E_m_parallel))[1]], sigdigits=4), " and t = ", round(t[argmax(abs.(E_m_ref - E_m_parallel))[2]], sigdigits=4), "\n")
        push!(E_spectral_radii, new_E_error / old_E_error)
        push!(E_m_spectral_radii, new_E_m_error / old_E_m_error)
        old_E_error = new_E_error
        old_E_m_error = new_E_m_error
        
        #p = plot(x_centers, E_coarse[:,end].^(1/4)/a, title="Final Temperatures at Iteration $(iterations)", label= "Radiation Temperature", xlabel="x", ylabel="T_r", xscale=:log10, color=:red, minorgrid=:true)
        #plot!(x_centers, T_coarse[:,end], label="Material Temperature", color=:blue, linestyle=:dash)
        p = plot(x_centers, abs.(mean!(ones(nx), E_ref-E_parallel)), title="Absolute Average Error at Iteration $(iterations)", label= "Radiation Energy", xlabel="x", ylabel="Absolute Average Change", yscale=:log10, color=:red, minorgrid=:true)
        plot!(x_centers, abs.(mean!(ones(nx), E_m_ref - E_m_parallel)), label="Material Energy", color=:blue, linestyle=:dash)
        ylims!(p, 1e-16, 1e1)
        #display(p)
        savefig(p, "outputs/frames/frame_$(lpad(iterations, 5, '0')).png")  # Save each frame
    end
    print("Parallel Scheme Converged in ", iterations, " iterations.\n")
    print("Mean Spectral Radius for E: ", mean(E_spectral_radii[2:end]), "\n")
    print("Mean Spectral Radius for E_m: ", mean(E_m_spectral_radii[2:end]), "\n")
    return E_parallel, E_m_parallel, T_parallel
end



function PararealTRTDiffusion_MG(params, inputs, threads, epsilon, E_ref, E_m_ref)
     # Extract fields from params struct
    (; dx, dt, nx, nt, ngroups, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, T, D_centers, D_edges, E_m, E, F, S, N_coarse) = params

    Tend = dt * nt
    dt_coarse = Tend / (threads * N_coarse)
    chunk_size = div(nt, threads)
    E_coarse = (a*c*Tr_init^4)*ones(nx, N_coarse*threads+1, ngroups)
    E_m_coarse = zeros(nx, N_coarse*threads+1)
    T_coarse = Tm_init*ones(nx, N_coarse*threads+1)
    F_coarse = zeros(nx+1, N_coarse*threads+1)
    E_m_coarse[:,1] .= E_m[:,1]
    S_coarse = zeros(nx, N_coarse*threads+1, ngroups)
    for g in 1:ngroups
        for n in 1:N_coarse*threads+1
            S_coarse[:,n,g] .= S[:,1,g]
        end
    end

    #print("Initial temperature check: ", T_init, "\n")
    #print("Initial radiation energy check: ", E[:,1], "\n")



    #coarse_params = (nx, threads*N_coarse, dx, Tend/(threads*N_coarse), c, a, T_src, rho, cV, sigma_a, beta, f, T_coarse, D_centers, D_edges, E_m_coarse, E_coarse, F_coarse)
    coarse_params = Mesh.Params(dx, dt_coarse, nx, threads*N_coarse, ngroups, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, T_coarse, D_centers, D_edges, E_coarse, E_m_coarse, F_coarse, S_coarse, N_coarse)
    E_coarse = E_coarse[:, 1:N_coarse:end, :] # Downsample coarse solution arrays
    E_m_coarse = E_m_coarse[:, 1:N_coarse:end]
    T_coarse = T_coarse[:, 1:N_coarse:end]
    #F_coarse = F_coarse[:, 1:N_coarse:end]
    #S_coarse = S_coarse[:, 1:N_coarse:end, :]

    E_parallel = (a*c*Tr_init^4)*ones(nx, nt+1, ngroups)
    E_m_parallel = zeros(nx, nt+1)
    T_parallel = Tm_init*ones(nx, nt+1)
    E_parallel[:,1,:] .= E[:,1,:]
    E_m_parallel[:,1] .= E_m[:,1]
    T_parallel[:,1] .= T[:,1]
    # Step 1: Coarse solution
    E_coarse, E_m_coarse, T_coarse = SerialTRTDiffusion_MG(coarse_params, inputs)
    E_coarse_old = copy(E_coarse)
    E_m_coarse_old = copy(E_m_coarse)
    T_coarse_old = copy(T_coarse)
    E_old = zeros(size(E_parallel))
    E_m_old = zeros(size(E_m_parallel))
    old_E_error = 1.0
    old_E_m_error = 1.0
    E_spectral_radii = Vector{Float64}()
    E_m_spectral_radii = Vector{Float64}()
    local_E = [zeros(nx, chunk_size+1, ngroups) for _ in 1:threads]
    local_E_m = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_T = [zeros(nx, chunk_size+1) for _ in 1:threads]
    local_F = [zeros(nx+1, chunk_size+1) for _ in 1:threads]
    local_f = [zeros(nx, ngroups) for _ in 1:threads]
    local_cV = [zeros(nx) for _ in 1:threads]
    iterations = 0
    print("Starting Parareal Iterations...\n")
    while maximum(abs.(E_ref - E_parallel)) > epsilon || maximum(abs.(E_m_ref - E_m_parallel)) > epsilon || iterations == 0
        iterations += 1
        #E_old .= copy(E_parallel)
        #E_m_old .= copy(E_m_parallel)

        Threads.@threads for p in 1:threads
            local_E[p][:,1,:] .= E_coarse[:,p,:]
            local_E_m[p][:,1] .= E_m_coarse[:,p]
            local_T[p][:,1] .= T_coarse[:,p]

            dl = zeros(nx-1)
            d = zeros(nx)
            du = zeros(nx-1)
            b = zeros(nx) # Source terms

            for n in 1:chunk_size
                # Matrix assembly for TRT-Diffusion equations
                 local_F[p][1,n] = (a*c*dt/4)*(T_src[1])^4
                 local_F[p][end,n] = (a*c*dt/4)*(T_src[2])^4

                if uppercase(inputs["cVmodel"]) == "LINEARIZED"
                    local_cV[p] = 4*a*local_T[p][:, n].^3 # Update specific heat for time-step
                else
                    local_cV[p] = cV # Use original cV if not linearized
                end
                local_beta = (4 * a * local_T[p][:, n].^3)./(rho.*local_cV[p])
                for g in 1:ngroups
                     local_f[p][:, g] = 1 ./(1 .+ c*dt*sigma_a[:, g].*local_beta)
                    # Left boundary
                    if inputs["LBC"] == "MARSHAK"
                        d[1] = 1 + c*dt*local_f[p][1,g]*sigma_a[1,g] + (c*dt/dx[1])*1.5*D_centers[1,g]*sigma_a[1,g]/(1 + 3*dx[1]*sigma_a[1,g]/4) + (c*dt/dx[1])*(D_edges[1,g]/dx[2])
                        du[1] = -(c*dt/dx[1])*(D_edges[1,g]/dx[2])
                        b[1] = c*dt*local_f[p][1,g]*sigma_a[1,g]*a*local_T[p][1,n]^4/ngroups + local_E[p][1,n,g] + 4*(1.5*D_centers[1,g]*sigma_a[1,g]/(1 + 3*dx[1]*sigma_a[1,g]/4))*local_F[p][1,n]/(c*dx[1]*ngroups) + dt*dx[1]*S[1,n,g]/(ngroups)
                    elseif inputs["LBC"] == "REFLECT"
                        # Zero flux condition for each group
                        d[1] = 1 + c*dt*local_f[p][1,g]*sigma_a[1,g] + c*dt*D_edges[1,g]/(dx[1]*dx[2])
                        du[1] = -c*dt*D_edges[1,g]/(dx[1]*dx[2])
                        b[1] = c*dt*local_f[p][1,g]*sigma_a[1,g]*a*local_T[p][1,n]^4/ngroups + local_E[p][1,n,g] + dt*dx[1]*S[1,n,g]/(ngroups)
                    end
                    
                    for j in 2:nx-1
                        dl[j-1] = -c*dt*D_edges[j-1,g]/(dx[j-1]*dx[j])
                        d[j] = 1 + c*dt*local_f[p][j,g]*sigma_a[j,g] + c*dt*(D_edges[j-1,g]/dx[j] + D_edges[j,g]/dx[j+1])/dx[j]
                        du[j] = -c*dt*D_edges[j,g]/(dx[j]*dx[j+1])
                        b[j] = c*dt*local_f[p][j,g]*sigma_a[j,g]*a*local_T[p][j,n]^4/ngroups + local_E[p][j,n,g] + dt*dx[j]*S[j,n,g]/(ngroups)
                    end
                    # Right boundary
                    if inputs["RBC"] == "MARSHAK"
                        dl[nx-1] = -(c*dt/dx[nx-1])*(D_edges[nx-1,g]/dx[nx])
                        d[nx] = 1 + c*dt*local_f[p][nx,g]*sigma_a[nx,g] + (c*dt/dx[nx])*(1.5*D_centers[nx,g]*sigma_a[nx,g]/(1 + 3*dx[nx]*sigma_a[nx,g]/4)) + (c*dt/dx[nx])*(D_edges[nx-1,g]/dx[nx])
                        b[nx] = c*dt*local_f[p][nx,g]*sigma_a[nx,g]*a*local_T[p][nx,n]^4/ngroups + local_E[p][nx,n,g] + 4*(1.5*D_centers[nx,g]*sigma_a[nx,g]/(1 + 3*dx[nx]*sigma_a[nx,g]/4))*local_F[p][nx+1,n]/(c*dx[nx]*ngroups) + dt*dx[nx]*S[nx,n,g]/(ngroups)
                    elseif inputs["RBC"] == "REFLECT"
                        # Zero flux condition for each group
                        dl[nx-1] = -c*dt*D_edges[nx-1,g]/(dx[nx-1]*dx[nx])
                        d[nx] = 1 + c*dt*local_f[p][nx,g]*sigma_a[nx,g] + c*dt*D_edges[nx-1,g]/(dx[nx-1]*dx[nx])
                        b[nx] = c*dt*local_f[p][nx,g]*sigma_a[nx,g]*a*local_T[p][nx,n]^4/ngroups + local_E[p][nx,n,g] + dt*dx[nx]*S[nx,n,g]/(ngroups)
                    end
                    TD = Tridiagonal(dl, d, du)
                    local_E[p][:,n+1,g] = TD \ b
                end
                local_E_m[p][:,n+1] = local_E_m[p][:,n] + dt*c*sum(local_f[p] .* sigma_a .* (local_E[p][:,n+1,:] .- (a.*local_T[p][:,n].^4 / ngroups)), dims=2)
              
                local_T[p][:,n+1] = local_T[p][:,n] + (local_E_m[p][:,n+1]-local_E_m[p][:,n]) ./ (rho.*local_cV[p])

                #local_T[p][:,n+1] = (local_T[p][:,n] + (c*dt*f.*sigma_a./(rho.*cV)).*(local_E[p][:,n+1] + 3*a*local_T[p][:,n].^4)) ./ (1 .+ 4*c*dt*f.*sigma_a.*a.*(local_T[p][:,n].^3)./(rho.*cV)) # Update Material Temperature

                #local_E_m[p][:,n+1] = rho.*cV.*local_T[p][:,n+1] # Update Material Energy
            end
        end
        # Step 2b: Combine local results into global arrays
        for g in 1:ngroups
            E_parallel[:,2:end,g] = hcat([local_E[p][:,2:end,g] for p in 1:threads]...)
        end
        E_m_parallel[:,2:end] = hcat([local_E_m[p][:,2:end] for p in 1:threads]...)
        T_parallel[:,2:end] = hcat([local_T[p][:,2:end] for p in 1:threads]...)
        # Step 2c: Perform prediction and correction step to update coarse solution with parallel fine solution
        for n in 2:threads+1
            #coarse_params = (nx, N_coarse, dx, Tend/(threads*N_coarse), c, a, T_src, rho, cV, sigma_a, beta, f, hcat(T_coarse[:,n-1], zeros(nx, N_coarse)), D_centers, D_edges, hcat(E_m_coarse[:,n-1], zeros(nx, N_coarse)), hcat(E_coarse[:,n-1], zeros(nx, N_coarse)), F_coarse)
            E_coarse_temp = zeros(nx, N_coarse+1, ngroups)
            E_coarse_temp[:,1,:] .= E_coarse[:,n-1,:]
            coarse_params = Mesh.Params(dx, dt_coarse, nx, N_coarse, ngroups, c, a, T_src, Tm_init, Tr_init, x_nodes, x_centers, t, rho, cV, sigma_a, beta, f, hcat(T_coarse[:,n-1], zeros(nx, N_coarse)), D_centers, D_edges, E_coarse_temp, hcat(E_m_coarse[:,n-1], zeros(nx, N_coarse)), F_coarse, S_coarse, N_coarse)
    
            E_coarse_new, E_m_coarse_new, T_coarse_new = SerialTRTDiffusion_MG(coarse_params, inputs)

            E_coarse[:,n,:] = E_coarse_new[:,end,:] .+ (E_parallel[:, Int(1+(n-1)*chunk_size),:] .- E_coarse_old[:,n,:])
            E_m_coarse[:,n] = E_m_coarse_new[:,end] .+ (E_m_parallel[:, Int(1+(n-1)*chunk_size)] .- E_m_coarse_old[:,n])
            T_coarse[:,n] = T_coarse_new[:,end] .+ (T_parallel[:, Int(1+(n-1)*chunk_size)] .- T_coarse_old[:,n])

            if any(E_coarse[:,n,:] .< 0)
                print("Negative radiation energy detected in coarse correction at iteration ", iterations, ", time chunk ", n, "\n")
                #return E_parallel, E_m_parallel, T_parallel, iterations
            elseif any(E_m_coarse[:,n] .< 0)
                print("Negative material energy detected in coarse correction at iteration ", iterations, ", time chunk ", n, "\n")
                #return E_parallel, E_m_parallel, T_parallel, iterations
            elseif any(T_coarse[:,n] .< 0)
                print("Negative temperature detected in coarse correction at iteration ", iterations, ", time chunk ", n, "\n")
                #return E_parallel, E_m_parallel, T_parallel, iterations
            end

            E_coarse_old[:,n,:] .= E_coarse_new[:,end,:]
            E_m_coarse_old[:,n] .= E_m_coarse_new[:,end]
            T_coarse_old[:,n] .= T_coarse_new[:,end]
        end

        #E_coarse .= E_coarse .+ (E_parallel[:,1:Int(chunk_size/N_coarse):end] .- E_coarse)
        #E_m_coarse .= E_m_parallel[:,1:Int(chunk_size/N_coarse):end]
        #T_coarse .= T_parallel[:,1:Int(chunk_size/N_coarse):end]
        
        # for k in 1:threads
        #     print(k , " ", sum(local_E[k][:,:,:]), " ", sum(local_E_m[k][:,:,:]), " ", sum(local_T[k][:,:,:]), "\n")
        # end

        new_E_error = maximum(abs.(E_ref - E_parallel))
        new_E_m_error = maximum(abs.(E_m_ref - E_m_parallel))
        
        print("Iteration ", iterations, ": ϵ_R = ", new_E_error, " with spectral radius: ", new_E_error / old_E_error ," at x = ", round(x_centers[argmax(abs.(E_ref - E_parallel))[1]], sigdigits=4), " and t = ", round(t[argmax(abs.(E_ref - E_parallel))[2]], sigdigits=4), "\n")
        print("Iteration ", iterations, ": ϵ_M = ", new_E_m_error, " with spectral radius: ", new_E_m_error / old_E_m_error , " at x = ", round(x_centers[argmax(abs.(E_m_ref - E_m_parallel))[1]], sigdigits=4), " and t = ", round(t[argmax(abs.(E_m_ref - E_m_parallel))[2]], sigdigits=4), "\n")
        push!(E_spectral_radii, new_E_error / old_E_error)
        push!(E_m_spectral_radii, new_E_m_error / old_E_m_error)
        old_E_error = new_E_error
        old_E_m_error = new_E_m_error
        #p = plot(x_centers, E_coarse[:,end].^(1/4)/a, title="Final Temperatures at Iteration $(iterations)", label= "Radiation Temperature", xlabel="x", ylabel="T_r", xscale=:log10, color=:red, minorgrid=:true)
        #plot!(x_centers, T_coarse[:,end], label="Material Temperature", color=:blue, linestyle=:dash)
        p = plot(x_centers, abs.(mean!(ones(nx), E_ref[:,:,1]-E_parallel[:,:,1])), title="Absolute Average Error at Iteration $(iterations)", label= "Radiation Energy Group 1", xlabel="x", ylabel="Absolute Average Change", yscale=:log10, color=:red, minorgrid=:true)
        for g in 2:ngroups
            plot!(x_centers, abs.(mean!(ones(nx), E_ref[:,:,g]-E_parallel[:,:,g])), label= "Radiation Energy Group $g", xlabel="x", ylabel="Absolute Average Change", yscale=:log10, color=:red, minorgrid=:true)
        end
        plot!(x_centers, abs.(mean!(ones(nx), E_m_ref - E_m_parallel)), label="Material Energy", color=:blue, linestyle=:dash)
        ylims!(p, 1e-16, 1e1)
        display(p)
        savefig(p, "outputs/frames/frame_$(lpad(iterations, 5, '0')).png")  # Save each frame
    end
    print("Parallel Scheme Converged in ", iterations, " iterations.\n")
    print("Mean Spectral Radius for E: ", mean(E_spectral_radii[2:end]), "\n")
    print("Mean Spectral Radius for E_m: ", mean(E_m_spectral_radii[2:end]), "\n")
    return E_parallel, E_m_parallel, T_parallel
end


end