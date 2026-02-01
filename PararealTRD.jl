# Parareal Thermal Radiation Diffusion


module PararealTRD

include("TRD_Inputs.jl")
include("TRD_Solvers.jl")
include("TRD_Mesh.jl")
include("TRD_Outputs.jl")
include("TRD_MonteCarlo.jl")


function main(args)

    println("Running Parareal Thermal Radiation Diffusion Simulation...")
    
    if isempty(args)
        print("No input file provided, exiting... \n") # Exception if no input file is provided
        return
    else
        input_file = args[1] # Use user provided input file
    end

    # Reading Input File
    inputs = Inputs.readInputs(input_file)

    # Generating Mesh Quantities
    params = Mesh.mesh_generation(inputs)

    # TRT Diffusion Solver Call
    if uppercase(inputs["solver"]) == "DETERMINISTIC"
        
        if uppercase(inputs["mode"]) == "SERIAL"
            if params.ngroups > 1
                #E_mg = zeros(params.nx, params.ngroups)
                print("Running Serial Multi-Group TRT Diffusion Solver... \n")
                params.E, params.E_m, params.T = Solvers.SerialTRTDiffusion_MG(params, inputs)
            else
                print("Running Serial Single-Group TRT Diffusion Solver... \n")
                params.E, params.E_m, params.T = Solvers.SerialTRTDiffusion(params, inputs)
            end
        elseif uppercase(inputs["mode"]) == "PARAREAL"
            if params.ngroups > 1
                print("Running Serial Reference Multi-Group TRT Diffusion Solver... \n")
                E_ref, E_m_ref, T_ref = Solvers.SerialTRTDiffusion_MG(params, inputs)
                print("Running Parareal Multi-Group TRT Diffusion Solver... \n")
                params.E, params.E_m, params.T = Solvers.PararealTRTDiffusion_MG(params, inputs, Threads.nthreads(), parse(Float64, inputs["epsilon"]), E_ref, E_m_ref)
            else
                print("Running Serial Reference Single-Group TRT Diffusion Solver... \n")
                E_ref, E_m_ref, T_ref = Solvers.SerialTRTDiffusion(params, inputs)
                print("Running Parareal Single-Group TRT Diffusion Solver... \n")
                params.E, params.E_m, params.T = Solvers.PararealTRTDiffusion(params, inputs, Threads.nthreads(), parse(Float64, inputs["epsilon"]), E_ref, E_m_ref)
            end
        else
            print("Invalid mode specified in input file. Use 'SERIAL' or 'PARAREAL'. \n")
            return
        end
    elseif uppercase(inputs["solver"]) == "MC"
        if uppercase(inputs["mode"]) == "SERIAL"
            print("Running Serial Monte Carlo TRT Diffusion Solver... \n")
            params.E, params.E_m, params.T = MonteCarlo.MC_Main(params, inputs)
        elseif uppercase(inputs["mode"]) == "PARAREAL"
            print("Running Serial Reference Single-Group TRT Diffusion Solver... \n")
            E_ref, E_m_ref, T_ref = Solvers.SerialTRTDiffusion(params, inputs)
            print("Running Parareal Monte Carlo TRT Diffusion Solver... \n")
            params.E, params.E_m, params.T = MonteCarlo.Parareal_MC_Main(params, inputs, Threads.nthreads(), parse(Float64, inputs["epsilon"]), E_ref, E_m_ref)
        end
    end

    # Output Results
    Outputs.plotting(params, inputs)
   
    println("Simulation complete.")

end
 

main(ARGS)

end