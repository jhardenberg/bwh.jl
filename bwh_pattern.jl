# run with: julia --project --check-bounds=no -O3 -t 4 bwh_pattern.jl
# where -t 4 is the number of threads to use

const USE_GPU = false  # flag for GPU - does not work with Metal on Mac
const flag_netcdf = true  # flag for netcdf output
const flag_ani = false  # flag for animation

#using ImplicitGlobalGrid   # this will be for MPI
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
using ChangePrecision
using GPUArraysCore  #needed for fix for Metal
using Plots
using Printf
using NetCDF

@static if USE_GPU
    @initl_parallel_stencil(Metal, Float32, 2);
    const float_type = Float32
else
    @init_parallel_stencil(Threads, Float64, 2);
    const float_type = Float64
end

include("integrate_zelnik.jl")
include("ghost.jl")
include("plotting.jl")
include("netcdf.jl")

function to_netcdf(fname, b, w, tim; mode="a")
    data = Dict(
        "b" => Dict("data" => Array(b), "long_name" => "biomass density"),
        "w" => Dict("data" => Array(w), "long_name" => "water density")
    )
    write_to_netcdf2d(fname, data, tim, mode=mode)
end

function bwh()

    println("Running with ", Threads.nthreads(), " threads")

    @changeprecision float_type begin  # fix for Metal (GPU on Mac). All in single precision.

    # Filenames
    filename_nc = "zelnik.nc"
    filename_anim = "zelnik.gif"
    
    # Zelnik et al. parameters

    # Physics
    η = 2.8             # root augmentation
    λ = 0.4571428          # soil water consumption rate
    ρ = 0.7             # shading parameter
    ν = 1.470588          # soil water evaporation rate

    db = 1.0            # b diffusivity
    dw = 125.0        # w diffusivity

    p = 1.55           # precipitation rate

    lx, ly = 340, 340   # Length of domain in dimensions x, y

    # Numerics
    nx, ny = 340+2, 340+2   # Number of gridpoints dimensions x, y 

    # this will be for MPI
    #me, dims, nprocs, coords   = init_global_grid(nx, ny, 1, periodx=1, periody=1);

    dx, dy = lx/(nx-2), ly/(ny-2)  # Space steps in x and y dimensions

    # Array initializations
    b = 0.6 *@rand(nx, ny) .+ 0.1
    w = @zeros(nx, ny) .+ 1.5

    b2 = @zeros(nx, ny);  # Temporary array for b
    w2 = @zeros(nx, ny);  # Temporary array for w

    # Time loop

    nt         = 150000       # Number of time steps
    nout = 5000  # how often to print stats
    nouta = 500  # how often to save animation
    noutf = 500  # how often to save netcdf file

    dt = min(dx^2,dy^2)/dw/8.1   # estimate time step
    println("Estimated dt = ", dt)
    dt = 0.001
    println("Selected dt = ", dt)

    println("Simulation parameters: p=", p, " η=", η, " λ=", λ, " ρ=", ρ, " ν=", ν, " db=", db, " dw=", dw, " dt=", dt, " dx=", dx, " dy=", dy)

    # Preparation of visualisation
    if flag_ani
        ENV["GKSwstype"]="nul";
        if isdir("viz2D_out")==false mkdir("viz2D_out") end; loadpath = "./viz2D_out/";
        anim = Animation(loadpath,String[])
        println("Animation directory: $(anim.dir)")
    end

    if flag_netcdf
        if isfile(filename_nc)  # remove file if it exists
            rm(filename_nc)
        end
    end

    wtime0 = Base.time()  # start timer for performance measurement

    update_ghost_serial!(b)   # initialize ghost points
    update_ghost_serial!(w)

    for it = 1:nt
        @parallel update_b!(b2, b, w, p, η, λ, ρ, ν, db, dw, dt, dx, dy)
        @parallel update_w!(w2, b, w, p, η, λ, ρ, ν, db, dw, dt, dx, dy)

        update_ghost_serial!(b2)  # this seems to be faster than the parallel version
        update_ghost_serial!(w2)
        #update_ghost!(b2, nx, ny)
        #update_ghost!(w2, nx, ny)

        b, b2 = b2, b
        w, w2 = w2, w

        if mod(it,nout)==0
            @printf("t = %7.2f b = [%10.8f, %10.8f]  w = [%10.8f, %10.8f]\n", it*dt, minimum(b), maximum(b), minimum(w), maximum(w))
        end

        if flag_netcdf && mod(it,noutf)==0
            to_netcdf(filename_nc, b, w, it*dt)
        end

        if flag_ani && mod(it,nouta)==0
            plotb(b, dt*it, nx, ny, lx, ly); frame(anim)
        end
    end

    #finalize_global_grid();  # will be for MPI

    # Performance
    wtime    = Base.time()-wtime0
    A_eff    = (2*2)/1e9*nx*ny*sizeof(Data.Number)  # Effective main memory access per iteration [GB]
    wtime_it = wtime/nt                        # Execution time per iteration [s]
    T_eff    = A_eff/wtime_it                       # Effective memory throughput [GB/s]
    @printf("Total steps=%d, time=%1.3e sec (@ T_eff = %1.2f GB/s) \n", nt, wtime, round(T_eff, sigdigits=2))
    
    plotb(b, dt*nt, nx, ny, lx, ly)
    savefig("zelnik.png")
    to_netcdf("final.nc", b, w, nt*dt, mode="c")

    if flag_ani
        gif(anim, filename_ani, fps = 15)
    end

    # plotbwh(Array(b), Array(w), nx, ny, dt*nt)
    end #  Changeprecision
    return b, w

end # function

@allowscalar b, w = bwh(); # Fix for Meta
