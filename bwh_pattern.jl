# On CPU run with: julia --project --check-bounds=no -O3 -t 4 bwh_pattern.jl
# where -t 4 is the number of threads to use

# When running with Metal (on Mac Silicon) use
# julia --project -O3 -t 4 bwh_pattern.jl

const USE_GPU = false  # flag for GPU - does not work with Metal on Mac
const flag_netcdf = false  # flag for netcdf output full data (final always produced)
const flag_ani = true  # flag for animation
const flag_read_init = false # flag to read the initial condition from file
const flag_disturbance_b = false  # flag for network disturbance on biomass diffusion
const flag_disturbance_w = true # flag for network disturbance on water diffusion
const flag_read_disturbance_b =  false # flags to read from file network disturbances on biomass diffusion
const flag_read_disturbance_w = false # flag to read from file network disturbance on water diffusion
const flag_add_disturbances_b =false #flags to add new links to the ones read from file on biomass diffusion
const flag_add_disturbances_w = false #flags to add new links to the ones read from file on water diffusion
const flag_weights_b = false #flag to consider weights in network disturbances for biomass
const flag_weights_w = false #flag to consider weights in network disturbances for water
const flag_smoothing_w = false #flag to approximate the diffusion over shortcuts to a global coupling term in water density
const flag_smoothing_b = false #flag to approximate the diffusion over shortcuts to a global coupling term in water density

#using ImplicitGlobalGrid   # this will be for MPI
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
using ChangePrecision
using GPUArraysCore  #needed for fix for Metal
using Plots
using Printf
using Logging
using DelimitedFiles
using Parameters
using Statistics

@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2);
    const float_type = Float64
else
    @init_parallel_stencil(Threads, Float64, 2);
    const float_type = Float64
end

include("integrate_zelnik.jl")
include("ghost.jl")
include("plotting.jl")
include("netcdf.jl")
include("disturbance.jl")
include("util.jl")
include("params.jl")

function to_netcdf(fname, b, w, tim; mode="a")
    data = Dict(
        "b" => Dict("data" => Array(b), "long_name" => "biomass density"),
        "w" => Dict("data" => Array(w), "long_name" => "water density")
    )
    write_to_netcdf2d(fname, data, tim, mode=mode)
end

@views @allowscalar function bwh(P)

    @changeprecision float_type begin  # fix for Metal (GPU on Mac). All in single precision.


    set_loglevel(P.loglevel)

    @info "Running with $(Threads.nthreads()) threads"

    nx, ny = P.numx + 2, P.numy + 2   # add 2 to each dimension for ghost cells
    dx, dy = P.lx/P.numx, P.ly/P.numy  # Space steps in x and y dimensions
    dt = P.dtstep  # we have to copy it for crazy julia scope restrictions
    
    #set initial conditions
    if flag_read_init
        # Array initialization. Read from file.
        data=NCDataset(P.filename_init)
        b = @zeros(nx, ny);
        b.=Data.Array(data["b"])
        w = @zeros(nx, ny);
        w.=Data.Array(data["w"])
    else
        # Array initializations. Use @zeros, @ones or @rand
        b = P.b_rand *(@rand(nx, ny) .- 0.5) .+ P.b_mean
        w = P.w_rand *(@rand(nx, ny) .- 0.5) .+ P.w_mean
        #b0=mean(b)
        #w0=mean(w)
        #@info "Random initialization array <b>=$b0, <w>=$w0"
    end 

    

    b2 = @zeros(nx, ny);  # Temporary array for b
    w2 = @zeros(nx, ny);  # Temporary array for w

    dte = min(dx^2,dy^2)/P.dw/8.1   # estimate time step
    if dt == 0.0  # if dt is not set, use the estimated value
        dt = dte
    end
    @info "Estimated dt = $dte"
    @info "Selected  dt = $dt"
    nt=P.nt/dt

    # this will be for MPI
    #me, dims, nprocs, coords   = init_global_grid(nx, ny, 1, periodx=1, periody=1);

    @info "Simulation parameters: p=$(P.p), η=$(P.η), λ=$(P.λ) ρ=$(P.ρ) ν=$(P.ν) db=$(P.db) dw=$(P.dw) dt=$dt dx=$dx dy=$dy"
    
    #read and/or make shortcuts biomass
    if flag_disturbance_b
        if flag_read_disturbance_b
            shortcuts = readdlm(P.filename_shortcuts_init_b)
            n_x_b=@zeros(nx,ny)
            n_y_b=@zeros(nx,ny)
            w_l_b=@zeros(nx,ny)
            n_x_b .= Data.Array(reshape(shortcuts[:,1],nx,ny))
            n_y_b .= Data.Array(reshape(shortcuts[:,2],nx,ny))
            w_l_b .= Data.Array(reshape(shortcuts[:,3],nx,ny))
            if flag_add_disturbances_b
                n_x_b1=n_x_b
                n_y_b1=n_y_b
                n_x_b2, n_y_b2, w_l_b2 = make_disturbed_links(nx, ny, dx, dy, P.ϕ)
                n_x_b, n_y_b=add_disturbances(n_x_b1,n_y_b1, n_x_b2, n_y_b2, nx,ny)
            end
        else
            n_x_b, n_y_b, w_l_b = make_disturbed_links(nx, ny, dx, dy, P.ϕ)
            if !flag_weights_b
                w_l_b = @ones(nx, ny)/(dx*dy)  # Do not weight with distance
            end
        end
        shortcuts_save = reshape(cat(n_x_b, n_y_b, w_l_b, dims=3),nx*ny,3)
        @allowscalar writedlm(P.filename_shortcuts_b, shortcuts_save , " ")
    end
    #read and/or make shortcuts water
    if flag_disturbance_w
        if flag_read_disturbance_w
            shortcuts = readdlm(P.filename_shortcuts_init_w)
            n_x_w=@zeros(nx,ny)
            n_y_w=@zeros(nx,ny)
            w_l_w=@zeros(nx,ny)
            n_x_w .= Data.Array(reshape(shortcuts[:,1],nx,ny))
            n_y_w .= Data.Array(reshape(shortcuts[:,2],nx,ny))
            w_l_w .= Data.Array(reshape(shortcuts[:,3],nx,ny))
            if flag_add_disturbances_w
                n_x_w2, n_y_w2, w_l_w2 = make_disturbed_links(nx, ny, dx, dy, P.ϕ)
                n_x_w, n_y_w=add_disturbances(n_x_w,n_y_w, n_x_w2, n_y_w2, nx,ny)
            end
        else
            n_x_w, n_y_w, w_l_w = make_disturbed_links(nx, ny, dx, dy, P.ϕ)
            if !flag_weights_w
                w_l_w = @ones(nx, ny)/(dx*dy)  # Do not weight with distance
            end
        end
        shortcuts_save = reshape(cat(n_x_w, n_y_w, w_l_w, dims=3),nx*ny,3)
        @allowscalar writedlm(P.filename_shortcuts_w, shortcuts_save , " ")
    end

    # Preparation of visualisation
    if flag_ani
        ENV["GKSwstype"]="nul";
        if isdir(P.animation_folder)
            rm(P.animation_folder; recursive=true)
        end
        mkdir(P.animation_folder)
        loadpath = P.animation_folder*"/"
        anim = Animation(loadpath,String[])
        @info "Animation directory: $(anim.dir)"
    end
    #set up for full data file if requested 
    if flag_netcdf
        if isfile(P.filename_nc)  # remove file if it exists
            rm(P.filename_nc)
        end
    end

    wtime0 = Base.time()  # start timer for performance measurement

    update_ghost_serial!(b)   # initialize ghost points
    update_ghost_serial!(w)
    #begin simulation
    for it = 1:nt
        @parallel update_b!(b2, b, w, P.p, P.η, P.λ, P.ρ, P.ν, P.db, P.dw, dt, dx, dy)
        @parallel update_w!(w2, b, w, P.p, P.η, P.λ, P.ρ, P.ν, P.db, P.dw, dt, dx, dy)
        #add globalcoupling terms in b
        if flag_smoothing_b
            if USE_GPU
                mean_b=mean(Array(b))
            else
                mean_b = mean(vec(b))
            end
            @parallel (2:(nx-1), 2:(ny-1)) smoothing!(b2, b, mean_b, P.ϕ, P.db, dt, dx, dy)
        end
        #add global coupling terms in w
        if flag_smoothing_w
            if USE_GPU
                mean_w=mean(Array(w))
            else
                mean_w = mean(vec(w))
            end
            @parallel (2:(nx-1), 2:(ny-1)) smoothing!(w2, w, mean_w, P.ϕ, P.dw, dt, dx, dy)
        end
        #add diffusion over shortcuts in b
        if flag_disturbance_b && P.ϕ>0
            @parallel (2:(nx-1), 2:(ny-1)) disturbance!(b2, b, n_x_b, n_y_b, w_l_b, P.db, dt, dx, dy)
        end
        #add diffusion over shortcuts in w
        if flag_disturbance_w && P.ϕ>0
            @parallel (2:(nx-1), 2:(ny-1)) disturbance!(w2, w, n_x_w, n_y_w, w_l_w, P.dw, dt, dx, dy)
        end
        
        #update_ghost_serial!(b2)  # this seems to be faster than the parallel version (not working on GPU)
        #update_ghost_serial!(w2)
        update_ghost!(b2, nx, ny)
        update_ghost!(w2, nx, ny)

        b, b2 = b2, b
        w, w2 = w2, w
        #print on the log
        if mod(it,P.nout)==0
            #@printf("t = %7.2f b = [%10.8f, %10.8f]  w = [%10.8f, %10.8f]\n", it*dt, minimum(b), maximum(b), minimum(w), maximum(w))
            @printf("t = %7.2f\n", it*dt)
        end
        #write on full data file
        if flag_netcdf && mod(it,P.noutf)==0
            to_netcdf(P.filename_nc, b, w, it*dt)
        end
        #add frame to animation
        #first step added to animation to verify initial condition
        if flag_ani && it==1
            plotb(b, dt*it, nx, ny, P.lx, P.ly); frame(anim)
        end
        
        if flag_ani && mod(it,P.nouta)==0
            plotb(b, dt*it, nx, ny, P.lx, P.ly); frame(anim)
        end
    end

    #finalize_global_grid();  # will be for MPI

    # Performance
    wtime    = Base.time()-wtime0
    A_eff    = (2*2)/1e9*nx*ny*sizeof(Data.Number)  # Effective main memory access per iteration [GB]
    wtime_it = wtime/nt                        # Execution time per iteration [s]
    T_eff    = A_eff/wtime_it                       # Effective memory throughput [GB/s]
    @info @sprintf("Total steps=%d, time=%1.3e sec (@ T_eff = %1.2f GB/s) \n", nt, wtime, round(T_eff, sigdigits=2))

    plotb(b, dt*nt, nx, ny, P.lx, P.ly)
    savefig(P.filename_final_img_b)
    plotw(w, dt*nt, nx, ny, P.lx, P.ly)
    savefig(P.filename_final_img_w)
    to_netcdf(P.filename_final_nc, b, w, nt*dt, mode="c")
    #compute animation
    if flag_ani
        gif(anim, P.filename_ani, fps = 20)
        rm(P.animation_folder ,recursive=true)
    end

    # plotbwh(Array(b), Array(w), nx, ny, dt*nt)
    end #  Changeprecision
    return b, w

end # function

#@allowscalar b, w = bwh(P); # Fix for Meta
