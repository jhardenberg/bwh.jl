# Description: Parameters for the Zelnik et al. model
  
# Filenames
filename_nc = "zelnik.nc"    # NetCDF filename
filename_anim = "zelnik.gif" # Animation filename

# Zelnik et al. parameters

# Physics
η = 2.8             # root augmentation
λ = 0.4571428       # soil water consumption rate
ρ = 0.7             # shading parameter
ν = 1.470588        # soil water evaporation rate
db = 1.0            # b diffusivity
dw = 125.0          # w diffusivity

p = 1.55            # precipitation rate

# Domain size
lx, ly = 340, 340   # Length of domain in dimensions x, y

# Numerics
nx, ny = 340 + 2, 340 + 2   # Number of gridpoints dimensions x, y - add 2 to each dimension for ghost cells

# Initialization
b_mean = 0.4  # mean value of b
b_rand = 0.6  # amplitude of perturbation
w_mean = 1.5  # mean value of w
w_rand = 0.0  # amplitude of perturbation

# Time loop
nt    = 150000       # Number of time steps
nout  = 5000  # how often to print stats
nouta = 500  # how often to save animation
noutf = 500  # how often to save netcdf file

dtstep = 0.001  # time step. If set to 0.0 the code will estimate the time step
