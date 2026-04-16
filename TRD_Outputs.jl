# Output module for Thermal Radiation Diffusion results

module Outputs

using Plots
using LaTeXStrings
using FFMPEG
using JLD2

include("TRD_Inputs.jl")

function plotting(params, inputs)
    # Extract necessary parameters
    x_centers = params.x_centers
    dx = params.dx
    t = params.t
    a = params.a
    E = params.E
    E_m = params.E_m
    T = params.T

    # Create output directories if they don't exist
    if !isdir("outputs")
        mkdir("outputs")
    end
    if !isdir("outputs/frames")
        mkpath("outputs/frames")
    end

    # Get the number of time steps available
    nt_available = size(E, 2) - 1  # -1 because we have nt+1 time points

    if uppercase(inputs["output"]) == "TEMPERATURE"
        output_quantity = "Temperature"
    elseif uppercase(inputs["output"]) == "ENERGY"
        output_quantity = "Energy Density"    
    end
    
    if uppercase(inputs["mode"]) == "SERIAL"
        
        title1 = "Serial Radiation $output_quantity vs Position"
        title2 = "Serial Material $output_quantity vs Position"
    else
        title1 = "Parareal Radiation $output_quantity vs Position"
        title2 = "Parareal Material $output_quantity vs Position"
    end

    if output_quantity == "Temperature"
        E_plot = (max.(E,0) ./ a) .^ (1/4) 
        Emat_plot = T
    elseif output_quantity == "Energy Density"
        E_plot = E ./ dx
        Emat_plot = (T .^4) ./ dx
    end

    # Save the mesh and simulation variables to a JLD2 file
    filename = "outputs\\" * inputs["name"] * "_" * titlecase(inputs["mode"]) * ".jld2"
    jldsave(filename; x_centers, t, E_plot, Emat_plot)

    if params.ngroups == 1     
        p1 = plot(x_centers, E_plot[:,end],  xlabel=L"x", ylabel=L"T_r", minorgrid=:true, label="t = $(t[end])")
        if nt_available >= 1000
            plot!(x_centers, E_plot[:,1001], label="t = $(t[1001])")
        end
        plot!(x_centers, E_plot[:,101], label="t = $(t[101])")
    else
        p1 = plot(x_centers, E_plot[:,end,1],  xlabel=L"x", ylabel=L"T_r", minorgrid=:true, label="t = $(t[end]) Group 1")
        plot!(x_centers, E_plot[:,end,2], label="t = $(t[end]) Group 2", linestyle=:dash)

        if nt_available >= 1000
            plot!(x_centers, E_plot[:,1001,1], label="t = $(t[1001]) Group 1")
            plot!(x_centers, E_plot[:,1001,2], label="t = $(t[1001]) Group 2", linestyle=:dash)
        end
        plot!(x_centers, E_plot[:,1,1], label="t = $(t[1]) Group 1")
        plot!(x_centers, E_plot[:,1,2], label="t = $(t[1]) Group 2", linestyle=:dash)
    end

    display(p1)
    p2 = plot(x_centers, Emat_plot[:,end],  xlabel=L"x", ylabel=L"T_m", minorgrid=:true, label="t = $(t[end])")
    
    # Add intermediate time points if they exist
    if nt_available >= 1000
        plot!(x_centers, Emat_plot[:,1001], label="t = $(t[1001])")
    end
    plot!(x_centers, Emat_plot[:,101], label="t = $(t[101])")
    display(p2)
    if uppercase(inputs["mode"]) == "SERIAL"
        savefig(p1, "outputs/TRT_Diffusion_RadTemp_Serial.png")
        savefig(p2, "outputs/TRT_Diffusion_MatTemp_Serial.png")
    else
        savefig(p1, "outputs/TRT_Diffusion_RadTemp_Parallel.png")
        savefig(p2, "outputs/TRT_Diffusion_MatTemp_Parallel.png")
    end

   

imagesdirectory = "outputs/frames/frame"
framerate = 5
gifname = "outputs/TRT_Temperature_Iterations.gif"

# Check if frame files exist before trying to create GIF
frame_files = filter(f -> endswith(f, ".png"), readdir("outputs/frames"))
if !isempty(frame_files)
    println("Creating GIF from $(length(frame_files)) frame files...")
    
    # Step 1: Generate palette
    FFMPEG.ffmpeg_exe(`-framerate $(framerate) -i $(imagesdirectory)_%05d.png -vf palettegen -y palette.png`)

    # Step 2: Create GIF using palette
    FFMPEG.ffmpeg_exe(`-framerate $(framerate) -i $(imagesdirectory)_%05d.png -i palette.png -lavfi paletteuse -y $(gifname)`)
    
    println("GIF created successfully: $(gifname)")
else
    println("No frame files found, skipping GIF creation")
end

# Clean up frame files
if isdir("outputs/frames")
    for fname in readdir("outputs/frames")
        if endswith(fname, ".png")
            rm(joinpath("outputs/frames", fname))       
        end
    end
end


end

end