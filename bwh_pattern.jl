# run with: julia --project --check-bounds=no -O3 -t 4 bwh_pattern.jl
# where -t 4 is the number of threads to use

const USE_GPU = false  # flag for GPU - does not work with Metal on Mac
const GPU_TYPE = "Metal"  # "CUDA", "AMDGPU" or "Metal"
const flag_netcdf = false  # flag for netcdf output
const flag_ani = false  # flag for animation
const flag_disturbance = false  # flag for network disturbance

#using ImplicitGlobalGrid   # this will be for MPI
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
using ChangePrecision
using GPUArraysCore  #needed for fix for Metal
using Plots
using Printf
using Logging

@static if USE_GPU
    if GPU_TYPE == "CUDA"
        @init_parallel_stencil(CUDA, Float64, 2);
        const float_type = Float64
    elseif GPU_TYPE == "AMDGPU"
        @init_parallel_stencil(AMDGPU, Float64, 2);
        const float_type = Float64
    else
        @init_parallel_stencil(Metal, Float32, 2);
        const float_type = Float32
    end
else
    @init_parallel_stencil(Threads, Float64, 2);
    const float_type = Float64
end

include("integrate_zelnik.jl")
include("ghost.jl")
include("plotting.jl")
include("netcdf.jl")
include("disturbance.jl")

function to_netcdf(fname, b, w, tim; mode="a")
    data = Dict(
        "b" => Dict("data" => Array(b), "long_name" => "biomass density"),
        "w" => Dict("data" => Array(w), "long_name" => "water density")
    )
    write_to_netcdf2d(fname, data, tim, mode=mode)
end

function bwh()

    logger = SimpleLogger(stdout, Logging.Debug)

    @info "Running with $(Threads.nthreads()) threads"

    @changeprecision float_type begin  # fix for Metal (GPU on Mac). All in single precision.

    include("params.jl")

    nx, ny = numx + 2, numy + 2   # add 2 to each dimension for ghost cells
    dx, dy = lx/numx, ly/numy  # Space steps in x and y dimensions
    dt = dtstep  # we have to copy it for crazy julia scope restrictions
    
    # Array initializations. Use @zeros, @ones or @rand
    b = b_rand *(@rand(nx, ny) .- 0.5) .+ b_mean
    w = w_rand *(@rand(nx, ny) .- 0.5) .+ w_mean

    b2 = @zeros(nx, ny);  # Temporary array for b
    w2 = @zeros(nx, ny);  # Temporary array for w

    dte = min(dx^2,dy^2)/dw/8.1   # estimate time step
    if dt == 0.0  # if dt is not set, use the estimated value
        dt = dte
    end
    @info "Estimated dt = $dte"
    @info "Selected dt = $dt"

    # this will be for MPI
    #me, dims, nprocs, coords   = init_global_grid(nx, ny, 1, periodx=1, periody=1);

    @info "Simulation parameters: p= $p, η=$η, λ=$λ ρ=$ρ ν=$ν db=$db dw=$dw dt=$dt dx=$dx dy=$dy"

    if flag_disturbance
        n_x, n_y, w_l = make_disturbed_links(nx, ny, dx, dy, ϕ)
        w_l = @ones(nx, ny)/(dx*dy)  # Do not weight with distance
    end

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

        if flag_disturbance
            @parallel (2:(nx-1), 2:(ny-1)) disturbance!(b2, b, n_x, n_y, w_l, db, dt, dx, dy)
        end

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
    @info @sprintf("Total steps=%d, time=%1.3e sec (@ T_eff = %1.2f GB/s) \n", nt, wtime, round(T_eff, sigdigits=2))
    
    plotb(b, dt*nt, nx, ny, lx, ly)
    savefig(filename_final_img)
    to_netcdf(filename_final_nc, b, w, nt*dt, mode="c")

    if flag_ani
        gif(anim, filename_ani, fps = 15)
    end

    # plotbwh(Array(b), Array(w), nx, ny, dt*nt)
    end #  Changeprecision
    return b, w

end # function

@allowscalar b, w = bwh(); # Fix for Meta
