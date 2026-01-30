# Parameters for the Zelnik et al. model
using Parameters
@with_kw mutable struct params

    dir_out::String = "/work/users/sfilippini/"
    # Filenames
    filename_nc = dir_out*"zelnik.nc"          # Full data NetCDF filename
    filename_ani = dir_out*"zelnik.gif"        # Animation filename
    filename_final_img_b = dir_out*"b_final.png"  # Final plot filename b
    filename_final_img_w = dir_out*"w_final.png"  # Final plot filename w
    filename_final_nc = dir_out*"final.nc"     # Final netcdf filename
    filename_shortcuts_b=dir_out*"shortcuts_b.dat" # File where to save the shortcuts (b)
    filename_shortcuts_init_b=dir_out*"shortcuts_b.dat" # File from where to read the shortcuts (b) if flag_read_disturbance=true
    filename_shortcuts_w=dir_out*"shortcuts_w.dat" # File where to save the shortcuts (w)
    filename_shortcuts_init_w=dir_out*"shortcuts_w.dat"  #File from where to read the shortcuts (w) if flag_read_disturbance=true
    filename_init=dir_out*"initial.nc"
    animation_file=dir_out*"viz2D_out"

    loglevel = "info"  # Logging level (debug, info, warn, error)

    # Zelnik et al. parameters

    # Physics
    η::Float64 = 2.8             # root augmentation
    λ::Float64 = 0.4571          # soil water consumption rate
    ρ::Float64 = 0.7             # shading parameter
    ν::Float64 = 1.4286          # soil water evaporation rate
    db::Float64 = 1.0            # b diffusivity
    dw::Float64= 125.0          # w diffusivity

    p::Float64 = 1.5            # precipitation rate

    # Network disturbance
    ϕb::Float64 = 0            # fraction of disturbed links
    ϕb_set::Float64 = 0       # fraction of preset disturbed links
    ϕw::Float64 = 0            # fraction of disturbed links
    ϕw_set::Float64 = 0       # fraction of preset disturbed links

    # Domain size
    lx::Float64 = 336
    ly::Float64 = 336  # Length of domain in dimensions x, y

    # Numerics
    numx::Int64 = 400
    numy::Int64 = 400  # Number of gridpoints in x, y

    # Initialization
    b_mean::Float64 = 0.4  # mean value of b
    b_rand::Float64 = 0.6  # amplitude of perturbation
    w_mean::Float64 = 1.5  # mean value of w
    w_rand::Float64 = 0.0  # amplitude of perturbation

    # Time loop
    dtstep::Float64 = 0.001  # time step. If set to 0.0 the code will estimate the time step
    nt::Int64    = 2500  # Number of adimensional times
    nout::Int64  = 1000    # how often to print stats (in time steps)
    nouta::Int64 = 5000    # how often to save animation (in time steps)
    noutf::Int64 = 1000    # how often to save netcdf file (in time steps)
end
