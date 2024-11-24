using NCDatasets

function write_to_netcdf2d(filename::String, data_dict::Dict{String, Dict{String, Any}}, tim::Float64; mode="a")
    """
    Appends data from a dictionary to a NetCDF file with a time axis.

    Args:
        filename (String): The name of the NetCDF file.
        data_dict (Dict): A dictionary where keys are variable names and
                          values are dictionaries with "data" and attribute keys.
        tim (Float64): The time value for the new data.
        mode (String): The mode to open the NetCDF file. Default is "a" (append or create if not existing).
                       Use "c" to overwriute even if existing.
    """
    
    # Open or create the NetCDF file
    if mode == "c"
        if isfile(filename)  # Remove the file if it exists
            rm(filename)
        end
    else
        mode = isfile(filename) ? "a" : "c"
    end
    
    ds = NCDataset(filename, mode)

    try
        # Define spatial dimensions
        spatial_dims = ("x", "y")
        for (var_name, var_info) in data_dict
            data = var_info["data"]
            for (dim_name, size) in zip(spatial_dims, size(data))
                if !(dim_name in keys(ds.dim))
                    defDim(ds, dim_name, size)
                end
            end
        end

        # Define the time dimension
        if haskey(ds.dim, "time")
            time_var = ds["time"]
            current_size = size(time_var, 1) 
            time_var[current_size + 1] = tim  
        else
            defDim(ds, "time", Inf)  
            time_var = defVar(ds, "time", Float64, ("time",), attrib=Dict("units"=>"model time units"))
            time_var[1] = tim
        end

        # Handle each variable in the data dictionary
        for (var_name, var_info) in data_dict
            data = var_info["data"]
            
            attrs = Dict(k => v for (k, v) in var_info if k != "data")
              
            # Combine "time" and spatial dimensions
            all_dims = ("x", "y", "time")

            # Define or update the variable
            if haskey(ds, var_name)
                var = ds[var_name]
                current_size = size(var, 3)
                # Correctly index the variable when appending data
                var[:, :, current_size] = data  
            else
                var = defVar(ds, var_name, eltype(data), all_dims, attrib = attrs)
                var[:, :, 1] = data
            end
        end
    finally
        close(ds)
    end
end