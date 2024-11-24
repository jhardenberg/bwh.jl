# bwh.jl

This is a completely revisited software to integrate simple reaction-diffusion equations. 

It relies heavily on the use of the great [ParallelStencil](https://github.com/omlins/ParallelStencil.jl) package which allows to
write code which will then works in parallel using CPU (threads) or on GPU with
extremely high efficiency and scalability.

The current version is set up to integrate the equations by [Zelnick et al. 2015](https://doi.org/10.1073/pnas.1504289112).

A simple example can be run with:

```
julia --project --check-bounds=no -O3 -t 4 bwh_pattern.jl
```

This will run in parallel with 4 threads, in about 30 seconds on a Mac M2, integrating 150000 timesteps (up to time 150) integrating a domain of 340x340 points.

The code writes at the end a netcdf file and png image with the results (`zelnik.nc` and `zelnik.png`).

Optionally an animation gif file can be created by specifying `flag_ani = true` at the beginning.
Data can be saved to netcdf with `flag_netcdf = true`

Dependencies are: `NCDatasets ParallelStencil ChangePrecision GPUArraysCore Plots`
