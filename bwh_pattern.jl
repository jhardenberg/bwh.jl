# On CPU run with: julia --project --check-bounds=no -O3 -t 4 bwh_pattern.jl
# where -t 4 is the number of threads to use

# When running with Metal (on Mac Silicon) use
# julia --project -O3 -t 4 bwh_pattern.jl

const USE_GPU = true  # flag for GPU - does not work with Metal on Mac
const flag_netcdf = true  # flag for netcdf output full data (final always produced)
const flag_ani = true  # flag for animation
const flag_read_init = true # flag to read the initial condition from file
const flag_read_disturbance_b = false # flag to read from file network disturbance on biomass diffusion
const flag_read_disturbance_w = false # flag to read from file network disturbance on water diffusion
const flag_weights_b = false #flag to consider weights in network disturbances for biomass
const flag_weights_w = false #flag to consider weights in network disturbances for water
const flag_smoothing_w = false #flag to approximate the effect of disturbances in biomass to mean-field approx
const flag_smoothing_b = false #flag to approximate the effect of disturbances in water to mean-field approx
const flag_plot_shortcuts_b=false
const flag_plot_shortcuts_w=false

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
using DataStructures
using FFTW

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

function to_netcdf(fname, b_c, b_s, w, tim; mode="a")
    data = Dict(
        "b_c" => Dict("data" => Array(b_c), "long_name" => "biomass density clonal-species"),
        "b_s" => Dict("data" => Array(b_s), "long_name" => "biomass density seed-species"),
        "w" => Dict("data" => Array(w), "long_name" => "water density")
    )
    write_to_netcdf2d(fname, data, tim, mode=mode)
end

@views @allowscalar function bwh(P)

    t_start=time_ns()

    @changeprecision float_type begin  # fix for Metal (GPU on Mac). All in single precision.

    
    set_loglevel(P.loglevel)

    @info "Running with $(Threads.nthreads()) threads"

    nx, ny = P.numx + 2, P.numy + 2   # add 2 to each dimension for ghost cells
    dx, dy = P.lx/P.numx, P.ly/P.numy  # Space steps in x and y dimensions
    dt = P.dtstep  # we have to copy it for crazy julia scope restrictions
    
    if flag_read_init && isfile(P.filename_init)
        # Array initialization. Read from file.
        data=NCDataset(P.filename_init)
        b_c = @zeros(nx, ny);
        b_c.=Data.Array(data["b_c"])
        b_s = @zeros(nx, ny);
        b_s.=Data.Array(data["b_s"])
        w = @zeros(nx, ny);
        w.=Data.Array(data["w"])
    else
        # Array initializations. Use @zeros, @ones or @rand
        b_c = P.b_c_rand *(@rand(nx, ny) .- 0.5) .+ P.b_c_mean
        b_s = P.b_s_rand *(@rand(nx, ny) .- 0.5) .+ P.b_s_mean
        w = P.w_rand *(@rand(nx, ny) .- 0.5) .+ P.w_mean
        #b0=mean(b)
        #w0=mean(w)
        #@info "Random initialization array <b>=$b0, <w>=$w0"
    end 

    

    b2_c = @zeros(nx, ny);  # Temporary array for b clonal
    b2_s = @zeros(nx, ny);  # Temporary array for b seeds
    w2 = @zeros(nx, ny);  # Temporary array for w

    dte = min(dx^2,dy^2)/P.dw/8.1   # estimate time step
    if dt == 0.0  # if dt is not set, use the estimated value
        dt = dte
    end
    @info "Estimated dt = $dte"
    @info "Selected  dt = $dt"
    
    #produce Fourier transform of the distance kernel
    #g=(1/2πl^2)*exp(-(|x-x'|/l))
    xs = collect(1:nx) .- (nx+1)/2
    ys = collect(1:ny) .- (ny+1)/2
    X = reshape(xs, nx, 1) .* ones(1, ny)
    Y = ones(nx, 1) .* reshape(ys, 1, ny)
    r = sqrt.(X.^2 .+ Y.^2)
    g = exp.(-r ./ P.l_s)
    g_norm=sum(g, dims=(1,2))
    g = g ./ g_norm
    g_shift = fftshift(g)
    gf=fft(g_shift)
    if USE_GPU
        gf = CUDA.CuArray(gf)
    end


    nt=Int(round(P.nt/dt))
    t_drop=Int(round(P.t_drop/dt))
    tau_drought=Int(round(P.tau_drought/dt))
    tau_rate=P.tau_rate*dt
    tau_change=Int(round(P.tau_change/dt))

    #calculate precipitation array in case of drought model
    if P.p_mode=="constant" 
        p_array=ones(Float64,Int(nt),1)*P.p 
    elseif P.p_mode =="quadratic"
        p_array=vcat(ones(Float64,Int(t_drop),1)*P.p_high, 
        ones(Float64,Int(tau_drought),1)*P.p_low, 
        ones(Float64,Int(nt-t_drop-tau_drought),1)*P.p_high)
    elseif P.p_mode == "linear"
        length_transition=Int(round((P.p_high-P.p_low)/tau_rate))
        p_array=vcat(ones(Float64,Int(t_drop),1)*P.p_high, 
        collect(1:1:length_transition)*(-1)*tau_rate.+P.p_high, 
        ones(Float64,Int(tau_drought),1)*P.p_low, 
        collect(1:1:length_transition)*(tau_rate).+P.p_low, 
        ones(Float64,Int(nt-t_drop-tau_drought-2*length_transition),1)*P.p_high)
    elseif P.p_mode =="smoothed"
        time=collect(1:1:nt)
        midway_drop=t_drop+(3*tau_change)
        midway_rise=t_drop+(6*tau_change)+tau_drought+(6*tau_change)
        p_array=P.p_low.+((P.p_high-P.p_low)/2).*(tanh.((time.-midway_drop)./tau_change).-tanh.((time.-midway_rise)./tau_change).+2)
    end

    # this will be for MPI
    #me, dims, nprocs, coords   = init_global_grid(nx, ny, 1, periodx=1, periody=1);

    @info "Simulation parameters: p=$(P.p) γ_c=$(P.γ_c) γ_s=$(P.γ_s) η_c=$(P.η_c) η_s=$(P.η_s) λ=$(P.λ) μ=$(P.μ) ρ=$(P.ρ) ν=$(P.ν) s_c=$(P.s_c) s_s=$(P.s_s) dw=$(P.dw) dt=$dt dx=$dx dy=$dy"
    
    #read or make shortcuts
    if P.ϕb>0
        if flag_read_disturbance_b && isfile(P.filename_shortcuts_init_b)
            shortcuts = readdlm(P.filename_shortcuts_init_b)    
            n_x_b = reshape(shortcuts[:,1],nx,ny)
            n_y_b = reshape(shortcuts[:,2],nx,ny)
            w_l_b = reshape(shortcuts[:,3],nx,ny)
            if P.ϕb_add>0
                n_x_b2, n_y_b2, w_l_b2, dist_b = make_disturbed_links(nx, ny, dx, dy, P.ϕb_add, P.sigma_b)
                n_x_b, n_y_b=add_disturbances(n_x_b,n_y_b, n_x_b2, n_y_b2, nx,ny)
            end
            if P.ϕb_rem>0
                M_r=floor(Int, 2*P.ϕb_rem*nx*ny)
                n_x_b, n_y_b=remove_disturbances(n_x_b, n_y_b, nx, ny, M_r)
            end
        else
            n_x_b, n_y_b, w_l_b, dist_b = make_disturbed_links(nx, ny, dx, dy, P.ϕb, P.sigma_b)
            if !flag_weights_b
                w_l_b = @ones(nx, ny)/(dx*dy)  # Do not weight with distance
            end
        end
        if flag_plot_shortcuts_b
            plot_shortcuts(dist_b)
            savefig(P.filename_plot_shortcuts_b)
        end
        n_x_b = Data.Array(n_x_b)
        n_y_b = Data.Array(n_y_b)
        w_l_b = Data.Array(w_l_b)
        shortcuts_save = reshape(cat(n_x_b, n_y_b, w_l_b, dims=3),nx*ny,3)
        @allowscalar writedlm(P.filename_shortcuts_b, shortcuts_save , " ")
    end

    if P.ϕw>0
        if flag_read_disturbance_w && isfile(P.filename_shortcuts_init_w)
            shortcuts = readdlm(P.filename_shortcuts_init_w)
            n_x_w = reshape(shortcuts[:,1],nx,ny)
            n_y_w = reshape(shortcuts[:,2],nx,ny)
            w_l_w = reshape(shortcuts[:,3],nx,ny)
            if P.ϕw_add>0
                n_x_w2, n_y_w2, w_l_w2, dist_w = make_disturbed_links(nx, ny, dx, dy, P.ϕw_add, P.sigma_w)
                n_x_w, n_y_w=add_disturbances(n_x_w,n_y_w, n_x_w2, n_y_w2, nx,ny)
            end
            if P.ϕw_rem>0
                M_r=2*P.ϕw_rem*nx*ny
                n_x_w, n_y_w=remove_disturbances(n_x_w, n_y_w, nx, ny, M_r)
            end
        else
            n_x_w, n_y_w, w_l_w, dist_w = make_disturbed_links(nx, ny, dx, dy, P.ϕw, P.sigma_w)
            if !flag_weights_w
                w_l_w = @ones(nx, ny)/(dx*dy)  # Do not weight with distance
            end
        end
        if flag_plot_shortcuts_w
            plot_shortcuts(dist_w)
            savefig(P.filename_plot_shortcuts_w)
        end
        n_x_w = Data.Array(n_x_w)
        n_y_w = Data.Array(n_y_w)
        w_l_w = Data.Array(w_l_w)
        shortcuts_save = reshape(cat(n_x_w, n_y_w, w_l_w, dims=3),nx*ny,3)
        @allowscalar writedlm(P.filename_shortcuts_w, shortcuts_save , " ")
    end

    # Preparation of visualisation
    if flag_ani
        ENV["GKSwstype"]="nul";
        if isdir(P.animation_file_c)
            rm(P.animation_file_c; recursive=true)
        end
        if isdir(P.animation_file_s)
            rm(P.animation_file_s; recursive=true)
        end
        mkdir(P.animation_file_c)
        mkdir(P.animation_file_s)
        loadpath = P.animation_file_c*"/"
        anim_c = Animation(loadpath,String[])
        loadpath = P.animation_file_s*"/"
        anim_s = Animation(loadpath,String[])
        @info "Animation directory clonal: $(anim_c.dir)"
        @info "Animation directory seed: $(anim_s.dir)"
    end

    if flag_netcdf
        if isfile(P.filename_nc)  # remove file if it exists
            rm(P.filename_nc)
        end
    end

    wtime0 = Base.time()  # start timer for performance measurement

    update_ghost_serial!(b_c)   # initialize ghost points
    update_ghost_serial!(b_s)   # initialize ghost points
    update_ghost_serial!(w)

    for it = 1:nt
        p=p_array[Int(it)]
        if P.noise_amplitude>0
            noise_field = P.noise_amplitude *(@rand(nx, ny) .- 0.5)
            b_c=b_c+noise_field
            b_s=b_s+noise_field
        end

        bf_s=fft(fftshift(b_s))
        conv=bf_s.*gf
        b_diff=real.(ifftshift(ifft(conv)))

        @parallel update_b_clonal!(b2_c, b_c, w, P.η_c, P.γ_c, dt, dx, dy)
        @parallel update_b_seeds!(b2_s, b_s, w, b_diff, P.η_s, P.γ_s, P.λ, P.μ, P.k, P.s_s, dt, dx, dy)
        @parallel update_w!(w2, b_c, b_s, w, p, P.γ_c,P.γ_s, P.η_c, P.η_s, P.ρ, P.ν, P.dw, dt, dx, dy)

        if flag_smoothing_b
            if USE_GPU
                mean_b_c=mean(Array(b_c))
                mean_b_s=mean(Array(b_s))
            else
                mean_b_c = mean(vec(b_c))
                mean_b_s = mean(vec(b_s))
            end
            @parallel (2:(nx-1), 2:(ny-1)) smoothing!(b2_c, b_c, mean_b_c, P.ϕb, P.db, dt, dx, dy)
            @parallel (2:(nx-1), 2:(ny-1)) smoothing!(b2_s, b_s, mean_b_s, P.ϕb, P.db, dt, dx, dy)
        end

        if flag_smoothing_w
            if USE_GPU
                mean_w=mean(Array(w))
            else
                mean_w = mean(vec(w))
            end
            @parallel (2:(nx-1), 2:(ny-1)) smoothing!(w2, w, mean_w, P.ϕw, P.dw, dt, dx, dy)
        end

        if P.ϕb>0
            @parallel (2:(nx-1), 2:(ny-1)) disturbance!(b2_c, b_c, n_x_b, n_y_b, w_l_b, P.db, dt, dx, dy)
            @parallel (2:(nx-1), 2:(ny-1)) disturbance!(b2_s, b_s, n_x_b, n_y_b, w_l_b, P.db, dt, dx, dy)
        end

        if P.ϕw>0
            @parallel (2:(nx-1), 2:(ny-1)) disturbance!(w2, w, n_x_w, n_y_w, w_l_w, P.dw, dt, dx, dy)
        end

        #update_ghost_serial!(b2)  # this seems to be faster than the parallel version
        #update_ghost_serial!(w2)
        update_ghost!(b2_c, nx, ny)
        update_ghost!(b2_s, nx, ny)
        update_ghost!(w2, nx, ny)

        b_c, b2_c = b2_c, b_c
        b_s, b2_s = b2_s, b_s
        w, w2 = w2, w

        if mod(it,P.nout)==0
            #@printf("t = %7.2f b = [%10.8f, %10.8f]  w = [%10.8f, %10.8f]\n", it*dt, minimum(b), maximum(b), minimum(w), maximum(w))
            @printf("t = %7.2f\n", it*dt)
        end

        if flag_netcdf && mod(it,P.noutf)==0
            to_netcdf(P.filename_nc, b_c, b_s, w, it*dt)
        end

        if flag_ani && it==1
            plotb(b_c, dt*it, nx, ny, P.lx, P.ly, p); frame(anim_c)
            plotb(b_s, dt*it, nx, ny, P.lx, P.ly, p); frame(anim_s)
        end

        if flag_ani && mod(it,P.nouta)==0
            plotb(b_c, dt*it, nx, ny, P.lx, P.ly, p); frame(anim_c)
            plotb(b_s, dt*it, nx, ny, P.lx, P.ly, p); frame(anim_s)
        end
    end

    #finalize_global_grid();  # will be for MPI

    # Performance
    wtime    = Base.time()-wtime0
    A_eff    = (2*2)/1e9*nx*ny*sizeof(Data.Number)  # Effective main memory access per iteration [GB]
    wtime_it = wtime/nt                        # Execution time per iteration [s]
    T_eff    = A_eff/wtime_it                       # Effective memory throughput [GB/s]
    @info @sprintf("Total steps=%d, time=%1.3e sec (@ T_eff = %1.2f GB/s) \n", nt, wtime, round(T_eff, sigdigits=2))

    plotb(b_c, dt*nt, nx, ny, P.lx, P.ly, p_array[end])
    savefig(P.filename_final_img_b_c)
    plotb(b_s, dt*nt, nx, ny, P.lx, P.ly, p_array[end])
    savefig(P.filename_final_img_b_s)
    plotw(w, dt*nt, nx, ny, P.lx, P.ly, p_array[end])
    savefig(P.filename_final_img_w)
    to_netcdf(P.filename_final_nc, b_c, b_s, w, nt*dt, mode="c")

    if flag_ani
        gif(anim_c, P.filename_ani_c, fps = 20)
        rm(P.animation_file_c ,recursive=true)
        gif(anim_s, P.filename_ani_s, fps = 20)
        rm(P.animation_file_s ,recursive=true)
    end

    t_elapsed = (time_ns() - t_start)/1e9
    println("Execution time: ", t_elapsed, " s")

    # plotbwh(Array(b), Array(w), nx, ny, dt*nt)
    end #  Changeprecision
    return b_c,b_s, w

end # function

#@allowscalar b, w = bwh(P); # Fix for Meta
