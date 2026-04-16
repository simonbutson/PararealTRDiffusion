using Plots
using Statistics
using LaTeXStrings
using JLD2

# Su-Olson Benchmarks Thermal Radiation Diffusion Benchmark Plotter

xMG = [0.0, 0.20, 0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 3.50, 5.0, 7.0, 9.0]
xSO = [0.01000, 0.10000, 0.17783, 0.31623, 0.45000, 0.50000, 0.56234, 0.75000, 1.00000, 1.33352, 1.77828, 3.16228, 5.62341, 10.0000, 17.78279]

tMG = [0.1, 0.3, 1.0, 3.0]
tSO = [0.10000, 0.31623, 1.0000, 3.16228, 10.0000, 31.6228, 100.000]

U1MG = [0.03873 0.08265 0.16857 0.30028;
      0.03649 0.07923 0.16426 0.29531;
      0.02951 0.06891 0.15134 0.28045;
      0.02402 0.06109 0.14164 0.26934;
      0.01840 0.05282 0.13115 0.25719;
      0.01027 0.03884 0.11187 0.23393;
      0.00535 0.02790 0.09478 0.21218;
      0.00214 0.01785 0.07633 0.18718;
      0.00075 0.01098 0.06081 0.16454;
      0.00023 0.00649 0.04793 0.14418;
      0.00000 0.00368 0.03735 0.12593;
      0.00000 0.00000 0.00650 0.05219;
      0.00000 0.00000 0.00071 0.01901;
      0.00000 0.00000 0.00000 0.00395;
      0.00000 0.00000 0.00000 0.00062]

U2MG = [0.04578 0.11394 0.23845 0.44678;
      0.04501 0.10820 0.22479 0.42434;
      0.03684 0.08416 0.17908 0.35415;
      0.02293 0.05926 0.14007 0.29871;
      0.00902 0.03433 0.10051 0.24141;
      0.00084 0.00994 0.04988 0.15585;
      0.00004 0.00232 0.02364 0.09963;
      0.00000 0.00028 0.00882 0.05695;
      0.00000 0.00003 0.00323 0.03321;
      0.00000 0.00001 0.00125 0.02023;
      0.00000 0.00000 0.00057 0.01312;
      0.00000 0.00000 0.00004 0.00264;
      0.00000 0.00000 0.00000 0.00070;
      0.00000 0.00000 0.00000 0.00010;
      0.00000 0.00000 0.00000 0.00001]

VMG =  [0.00452 0.03326 0.20363 0.66656;
      0.00447 0.03200 0.19265 0.63301  
      0.00380 0.02552 0.15342 0.52613
      0.00229 0.01736 0.11740 0.43986
      0.00077 0.00919 0.08102 0.35113
      0.00011 0.00258 0.03862 0.22372
      0.00003 0.00087 0.01903 0.14403
      0.00001 0.00035 0.00882 0.08608
      0.00000 0.00018 0.00489 0.05476
      0.00000 0.00009 0.00314 0.03754
      0.00000 0.00005 0.00218 0.02760
      0.00000 0.00000 0.00025 0.00778
      0.00000 0.00000 0.00002 0.00228
      0.00000 0.00000 0.00000 0.00036
      0.00000 0.00000 0.00000 0.00004]    

USO = [0.09403 0.24356 0.50359 0.95968 1.86585 0.66600 0.35365
      0.09326 0.24002 0.49716 0.95049 1.85424 0.66562 0.35360
      0.09128 0.23207 0.48302 0.93036 1.82889 0.66479 0.35347
      0.08230 0.20515 0.43743 0.86638 1.74866 0.66216 0.35309
      0.06086 0.15981 0.36656 0.76956 1.62824 0.65824 0.35252
      0.04766 0.13682 0.33271 0.72433 1.57237 0.65643 0.35225
      0.03171 0.10856 0.29029 0.66672 1.50024 0.65392 0.35188
      0.00755 0.05086 0.18879 0.51507 1.29758 0.64467 0.35051
      0.00064 0.01583 0.10150 0.35810 1.06011 0.62857 0.34809
      0.00000 0.00244 0.04060 0.21309 0.79696 0.60098 0.34382
      0.00000 0.00000 0.01011 0.10047 0.52980 0.55504 0.33636
      0.00000 0.00000 0.00003 0.00634 0.12187 0.37660 0.30185
      0.00000 0.00000 0.00000 0.00000 0.00445 0.11582 0.21453
      0.00000 0.00000 0.00000 0.00000 0.00000 0.00384 0.07351
      0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00269
]

VSO = [0.00466 0.03816 0.21859 0.75342 1.75359 0.67926 0.35554
      0.00464 0.03768 0.21565 0.74557 1.74218 0.67885 0.35548
      0.00458 0.03658 0.20913 0.72837 1.71726 0.67796 0.35536
      0.00424 0.03253 0.18765 0.67348 1.63837 0.67517 0.35497
      0.00315 0.02476 0.15298 0.58978 1.51991 0.67100 0.35438
      0.00234 0.02042 0.13590 0.55041 1.46494 0.66907 0.35411
      0.00137 0.01515 0.11468 0.50052 1.39405 0.66640 0.35374
      0.00023 0.00580 0.06746 0.37270 1.19584 0.65656 0.35235
      0.00000 0.00139 0.03173 0.24661 0.96571 0.63947 0.34988
      0.00000 0.00015 0.01063 0.13729 0.71412 0.61022 0.34555
      0.00000 0.00000 0.00210 0.05918 0.46369 0.56166 0.33797
      0.00000 0.00000 0.00000 0.00281 0.09834 0.37513 0.30294
      0.00000 0.00000 0.00000 0.00000 0.00306 0.11060 0.21452
      0.00000 0.00000 0.00000 0.00000 0.00000 0.00334 0.07269
      0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00258
]
#  # Load Simulation Data from JLD2 file
 filename = "outputs\\SuOlson_Parareal.jld2"
 @load filename x_centers t E_plot Emat_plot


# p = plot(t, mean(E_plot, dims=1)[:], minorgrid=true, xlabel=L"t", ylabel=L"\bar{T}", label=L"\bar{T}_r", xscale=:log10, linewidth = 2.5, palette=:tol_bright)
# plot!(t, mean(Emat_plot, dims=1)[:], label=L"\bar{T}_m", linewidth = 2.5)
# xlims!(p, 5e-3, maximum(t))

# display(p)

# savefig("outputs\\Equilibrium_Parareal_MeanTemperatures.png")



      # p = plot(x_centers, E_plot[:,end], label=L"T_r" * " at t="*string(t[end]), minorgrid=true, palette=:tol_bright)
      # plot!(p, x_centers, E_plot[:,1001], label=L"T_r" * " at t="*string(t[1001]))
      # plot!(p, x_centers, E_plot[:,101], label=L"T_r" * " at t="*string(t[101]))
      # plot!(p, x_centers, E_plot[:,11], label=L"T_r" * " at t="*string(t[11]))
      # plot!(p, x_centers, E_plot[:,1], label=L"T_r" * " at t="*string(t[1]))
      # ylims!(p, 0.4, 1.1)

      # display(p)

      # savefig(p, "outputs\\Equilibrium_Parareal_RadiationTemperature.png")

#     p = scatter(xMG, U1MG[:,1], markershape=:circle, minorgrid=true, xscale =:log10, label=false, palette=:tol_bright)
#     plot!(p, x_centers, E_plot[:,1001,1], label=L"U_1" * " at t="*string(t[1001]), linewidth=2.5)
#     scatter!(p, xMG, U2MG[:,1], markershape=:rect, label=false)
#     plot!(p, x_centers, E_plot[:,1001,2], label=L"U_2" * " at t="*string(t[1001]), linewidth=2.5)
#     scatter!(p, xMG, U1MG[:,3], markershape=:diamond, label=false)
#     plot!(p, x_centers, E_plot[:,end,1], label=L"U_1" * " at t="*string(t[end]), linewidth=2.5)
#     scatter!(p, xMG, U2MG[:,3], markershape=:star, label=false)
#     plot!(p, x_centers, E_plot[:,end,2], label=L"U_2" * " at t="*string(t[end]), linewidth=2.5)
#     xlims!(p, 1e-1, 10)
#     ylabel!(p, L"U_g")
#     xlabel!(p, L"x")
#     display(p)

#     savefig(p, "outputs\\MultiGroup_Parareal_RadiationEnergyDensity.png")
  
#     p = scatter(xSO[1:end-1], USO[1:end-1,1], markershape=:rect, minorgrid=true, xscale=:log10, label=false, palette=:tol_bright)
#     plot!(p, x_centers, E_plot[:,1001], label=L"U" * " at t="*string(t[1001]), linewidth=2.5)
#     scatter!(p, xSO[1:end-1], USO[1:end-1,3], markershape=:diamond, label=false)
#     plot!(p, x_centers, E_plot[:,end], label=L"U" * " at t="*string(t[end]), linewidth=2.5)
#     ylabel!(p, L"U")
#     xlabel!(p, L"x")
#     #xlims!(p, 1e-3, 10)
#     display(p)

      #savefig(p, "outputs\\SuOlson_Parareal_RadiationEnergyDensity.png")

      # p2 = plot(x_centers, Emat_plot[:,end], label=L"T_m" * " at t="*string(t[end]), minorgrid=true)
      # plot!(p2, x_centers, Emat_plot[:,1001], label=L"T_m" * " at t="*string(t[1001]))
      # plot!(p2, x_centers, Emat_plot[:,101], label=L"T_m" * " at t="*string(t[101]))
      # plot!(p2, x_centers, Emat_plot[:,11], label=L"T_m" * " at t="*string(t[11]))
      # plot!(p2, x_centers, Emat_plot[:,1], label=L"T_m" * " at t="*string(t[1]))
      # ylims!(p2, 0.4, 1.1)
      # display(p2)
      # savefig(p2, "outputs\\Equilibrium_Parareal_MaterialTemperature.png")


#     p2 = scatter(xMG, VMG[:,1], markershape=:rect, minorgrid=true, xscale=:log10, yscale=:log10, label=false, palette=:tol_bright)
#     plot!(p2, x_centers, Emat_plot[:,1001], label=L"V" * " at t="*string(t[1001]), linewidth=2.5)
#     scatter!(p2, xMG, VMG[:,3], markershape=:diamond, label=false)
#     plot!(p2, x_centers, Emat_plot[:,end], label=L"V" * " at t="*string(t[end]), linewidth=2.5)
#     ylabel!(p2, L"V")
#     xlabel!(p2, L"x")
#     xlims!(p2, 1e-1, 6)
#     ylims!(p2, 1e-3, 2)
#     display(p2)

#savefig(p2, "outputs\\MultiGroup_Parareal_MaterialEnergyDensity.png")

    p2 = scatter(xSO[1:end-1], VSO[1:end-1,1], markershape=:rect, xscale=:log10, yscale=:log10, minorgrid=true, label=false, palette=:tol_bright)
    plot!(p2, x_centers, Emat_plot[:,1001], label=L"V" * " at t="*string(t[1001]), linewidth=2.5)
    scatter!(p2, xSO[1:end-1], VSO[1:end-1,3], markershape=:diamond, label=false)
    plot!(p2, x_centers, Emat_plot[:,end], label=L"V" * " at t="*string(t[end]), linewidth=2.5)
    ylabel!(p2, L"V")
    xlabel!(p2, L"x")
    ylims!(p2, 1e-3, 2)
    display(p2)

    savefig(p2, "outputs\\SuOlson_Parareal_MaterialEnergyDensity.png")

    


# Convergence Study Speed-Up Results

N_p = [5,10,20,25,50,100]
R = [10,20,50,100]

S_eq = [1.00 0.909 1.17 1.06 1.22 1.31
        0.952 1.25 1.54 1.47 1.85 2.13
        1.11 1.35 1.72 2.00 2.63 3.13
        1.18 1.28 1.64 1.96 2.70 3.70]


S_so = [0.769 1.11 1.18 1.25 1.43 1.32
        0.952 1.25 1.82 2.00 2.13 2.13
        0.892 1.61 2.27 2.27 2.94 3.13
        0.943 1.49 2.33 2.44 2.94 3.23]        

S_mg = [0.625 0.769 1.00 1.06 1.06 1.15
        0.645 1.05 1.33 1.47 1.85 1.89
        0.746 1.16 1.72 2.00 2.63 3.13
        0.787 1.12 1.82 2.17 2.94 3.70]   

# p = plot(N_p, S_eq[1,   :], markershape=:circle, minorgrid=true, linewidth=2.5, label=L"R"*" = 10", palette=:tol_bright)
# plot!(N_p, S_eq[2,:], markershape=:square, linewidth=2.5, label=L"R"*" = 20")
# plot!(N_p, S_eq[3,:], markershape=:diamond, linewidth=2.5, label=L"R"*" = 50")
# plot!(N_p, S_eq[4,:], markershape=:utriangle, linewidth=2.5, label=L"R"*" = 100")
# xlabel!(L"N_p")
# ylabel!(L"S")

# display(p)
# savefig("outputs\\Speed_Up_Equilibrium.png")

# p2 = plot(N_p, S_so[1,   :], markershape=:circle, minorgrid=true, linewidth=2.5, label=L"R"*" = 10", palette=:tol_bright)
# plot!(N_p, S_so[2,:], markershape=:square, linewidth=2.5, label=L"R"*" = 20")
# plot!(N_p, S_so[3,:], markershape=:diamond, linewidth=2.5, label=L"R"*" = 50")
# plot!(N_p, S_so[4,:], markershape=:utriangle, linewidth=2.5, label=L"R"*" = 100")
# xlabel!(L"N_p")
# ylabel!(L"S")

# display(p2)
#  savefig("outputs\\Speed_Up_Gray_Su_Olson.png")

# p3 = plot(N_p, S_mg[1,   :], markershape=:circle, minorgrid=true, linewidth=2.5, label=L"R"*" = 10", palette=:tol_bright)
# plot!(N_p, S_mg[2,:], markershape=:square, linewidth=2.5, label=L"R"*" = 20")
# plot!(N_p, S_mg[3,:], markershape=:diamond, linewidth=2.5, label=L"R"*" = 50")
# plot!(N_p, S_mg[4,:], markershape=:utriangle, linewidth=2.5, label=L"R"*" = 100")
# xlabel!(L"N_p")
# ylabel!(L"S")

# display(p3)
# savefig("outputs\\Speed_Up_Multi_Group.png")