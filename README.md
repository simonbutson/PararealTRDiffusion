# Parareal Thermal Radiation Diffusion Solver

This code package implements a 1D multigroup thermal radiation diffusion solver using the Parareal method. Options include deterministic, Monte Carlo, and hybrid implementations of the same solver.

The deterministic solver has been validated and the example problems from the "Parareal Thermal Radiation Diffusion" paper can be recreated using the input decks in the paper-inputs section.
Available problems in the inputs folder include material/radiation equilibration, a marshak wave, as well as gray and multigroup Su Olson benchmarks. \
The Monte Carlo and hybrid implementations are still experimental and undergoing further development.

The code can be called in the following way: \
julia --threads n_threads PararealTRD.jl input_folder\file_name.txt \
where n_threads is an integer number of threads to run with and input_folder\file_name is the appropriate path to an input deck.

## Input Decks
Here is a list of input deck keywords and available options: \
name = String: give the problem a name \
mode = SERIAL or PARAREAL: choose to run in serial or parareal mode \
solver: DETERMINISTIC, MC, HYBRID: use a determinisitic, monte carlo or hybrid (deterministic coarse and MC fine) solver \
epsilon = Float: enter a floating point value for the absolute Parareal convergence tolerance \
ngroups = Int: enter the number of energy groups \
nx = Int: number of spatial cells \
dx = Float: spatial cell width \
dt = Float: fine time-step length \
nt = Float: number of fine time-steps \
c = Float: speed of light \
a = Float: radiation constant \
rho = Float: density \
cVmodel = CONSTANT or LINEARIZED:  choose a constant heat capacity or to use the linearized model from the Su Olson benchmarks\
cV = Float: heat capacity constant \
sigma_a = Float: absoprtion opacity \
T_src = [Float, Float]: temperatures of sources on left and right boundaries, use [0.0, 0.0] if no sources are present \
Tm_init = Float: Initial material temperature \
Tr_init = Float: Initial radiation temperature \
S_regs = [Float, ..., Float]: Enter values for source region cut-offs. Sources will start at the last cut-off value and will end at the next cut-off value. The start x=0.0 does not need to be included, but the end of the slab x_end = dx*nx always has to be included even if only one region is present \
S_vals = [Float, ..., Float]: Enter values for volume source strengths corresponding to the regions in S_regs \
LBC = REFLECT or MARSHAK: Reflecting or Marshak boundary conditions at left boundary \
RBC = REFLECT or MARSHAK: Reflecting or Marshak boundary conditions at right boundary \
N_coarse = Int: Ratio of number of fine over coarse time-steps. Example N_coarse = 10 means 10x as many fine time-steps as coarse. The coarse time-steps are sized as dt_coarse = N_coarse * dt \
N_particles = Int: Number of particles per time-step used by Monte Carlo solver \
output = TEMPERATURE or ENERGY: Output values as radiation/material temperatures or as energy densities 
