# Parameters for the Zelnik et al. model
using Parameters
@with_kw mutable struct params

    dir_out::String = "./"
    # Filenames
    filename_nc = dir_out*"evolution.nc"          # Full data NetCDF filename
    filename_ani_c = dir_out*"evolution_bc.gif"        # Animation filename
    filename_ani_s = dir_out*"evolution_bs.gif"        # Animation filename
    filename_final_img_b_c = dir_out*"bc_final.png"  # Final plot filename b
    filename_final_img_b_s = dir_out*"bs_final.png"  # Final plot filename b
    filename_final_img_w = dir_out*"w_final.png"  # Final plot filename w
    filename_final_nc = dir_out*"final.nc"     # Final netcdf filename
    filename_shortcuts_b=dir_out*"shortcuts_b.dat" # File where to save the shortcuts (b)
    filename_shortcuts_init_b=dir_out*"shortcuts_b.dat" # File from where to read the shortcuts (b) if flag_read_disturbance=true
    filename_shortcuts_w=dir_out*"shortcuts_w.dat" # File where to save the shortcuts (w)
    filename_shortcuts_init_w=dir_out*"shortcuts_w.dat"  #File from where to read the shortcuts (w) if flag_read_disturbance=true
    filename_init=dir_out*"initial.nc"
    animation_file_c=dir_out*"viz2D_out_bc"
    animation_file_s=dir_out*"viz2D_out_bs"
    filename_plot_shortcuts_b=dir_out*"shortcut_distances_b.png" #plot of the shortcut distances in b if flag_plot_shortcuts_b=true
    filename_plot_shortcuts_w=dir_out*"shortcut_distances_w.png" #plot of the shortcut distances in w if flag_plot_shortcuts_w=true

    loglevel = "info"  # Logging level (debug, info, warn, error)

    # Zelnik et al. parameters

    # Physics
    ν::Float64 = 1.4286            #evaporation rate
    ρ::Float64 = 0.7               #shading factor
    γ_c::Float64 = 0.4571          #water absortion rate clonal
    γ_s::Float64 = 0.4571          #water absortion rate seeds
    λ::Float64 = 1                 #ratio of seeds/clonal growth rate per water absorbed
    μ::Float64 = 1                 #ratio of seeds/clonal mortality rate
    k::Float64=1                   #ratio of seeds/clonal carrying capacities
    η_c::Float64 = 2.8             # root augmentation clonal
    η_s::Float64 = 2.8             # root augmentation seeds
    s_s::Float64 = 1               #total dispersal rate outside of point x seeds
    l_s::Float64 = 30         #seed dispersal characteristic distance
    dw::Float64 = 125          #water diffusion coefficient

    
    #precipitation rate (fixed or variable)
    #p_modes: "constant"-> precipitation fixed at p
    #"quadratic" -> precipitation drops from p_high to p_low instantly at t_drop and recovers at t_drop+tau_drought
    #"linear" -> precipitation drops from p_high to p_low starting at t_drop with rate tau_rate and begins recovery at t_drop+tau_drought
    #transition lasts tau_rate*(p_high-p_low)
    #"smoothed" -> precipitation drops from p_high to p_low with a smoothed transition (accelerating in the middle). Transition starts at
    # t_drop and each lasts 6*tau_change. tau_drought is the length of time for which it is at minimum
    p_mode="constant"
    p = 1.5 
    p_high = 1.8
    p_low = 1.5
    t_drop = 0.0
    tau_drought=100
    tau_rate=0.01
    tau_change=10
    noise_amplitude= 0.0 #noise to added at every step of the simulation
    
    # Network disturbance
    ϕb::Float64 = 0           # fraction of disturbed links
    ϕb_add::Float64 = 0       # fraction of disturbed links to add to the ones read from file
    ϕb_rem::Float64 = 0       # fraction of preset disturbed links to remove from the ones read from file
    sigma_b::Float64= 1       # standard deviation of the gaussian weighting function expressed as fraction 
                              # of maximum distance on periodic domain (sqrt(nx^2+ny^2))
                              # when sigma=1 shortcut lengths are uniformly distributed 
    ϕw::Float64 = 0           # fraction of disturbed links
    ϕw_add::Float64 = 0       # fraction of preset disturbed links 
    ϕw_rem::Float64 = 0       # fraction of preset disturbed links to remove
    sigma_w::Float64= 1    

    # Domain size
    lx::Float64 = 400
    ly::Float64 = 400  # Length of domain in dimensions x, y

    # Numerics
    numx::Int64 = 400
    numy::Int64 = 400  # Number of gridpoints in x, y

    # Initialization
    b_c_mean::Float64 = 0.4  # mean value of b
    b_s_mean::Float64 = 0.4  # mean value of b
    b_c_rand::Float64 = 0.6  # amplitude of perturbation
    b_s_rand::Float64 = 0.6  # amplitude of perturbation
    w_mean::Float64 = 1.5  # mean value of w
    w_rand::Float64 = 0.0  # amplitude of perturbation

    # Time loop
    dtstep::Float64 = 0.001  # time step. If set to 0.0 the code will estimate the time step
    nt::Int64    = 2500  # Number of adimensional times
    nout::Int64  = 1000    # how often to print stats (in time steps)
    nouta::Int64 = 5000    # how often to save animation (in time steps)
    noutf::Int64 = 1000    # how often to save netcdf file (in time steps)
end
