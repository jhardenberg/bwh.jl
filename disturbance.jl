using Random

# Function for ParallelStencil to add disturbance to the field
@parallel_indices (ix, iy) function disturbance!(b2::Data.Array, b::Data.Array, n_x::Data.Array, n_y::Data.Array, ww::Data.Array, dd::Data.Number, dt::Data.Number, dx::Data.Number, dy::Data.Number)
    """
    Adds a disturbance to the field b2 at the indices ix and iy.
    Uses ParallelStencil index notation.

    Args:
        b2: The destination field to add the disturbance to.
        b: The original field at the previous timestep
        n_x: The x indices of the disturbed links.
        n_y: The y indices of the disturbed links.
        ww: The weights of the disturbed links based on distance
        dd: Diffusion parameter.
        dt: The time step.
        dx: The spatial step between nodes in the x direction.
        dy: The spatial step between nodes in the y direction.
    """
    iix = Int(n_x[ix, iy])
    iiy = Int(n_y[ix, iy])
    b2[ix, iy] = b2[ix, iy] + dt*dd*(b[iix, iiy] - b[ix, iy])*(ww[ix, iy])  # weights inverted and dx*dy included in w_t 
    return
end

function make_disturbed_links(nx, ny, dx, dy, ϕ)
    """
    Creates a set of disturbed links in the network.

    Args:
        nx: The number of nodes in the x direction.
        ny: The number of nodes in the y direction.
        dx: The spatial step between nodes in the x direction.
        dy: The spatial step between nodes in the y direction.
        ϕ: The fraction of disturbed links.

    Returns as ParallelStencil.Data.Array:
        n_x: The x indices of the disturbed links.
        n_y: The y indices of the disturbed links.
        w  : The weights of the disturbed links.
    """

    M =floor(Int, ϕ*nx*ny)  # number of disturbed links
    @info "ϕ = $ϕ --> Number of disturbed links: $M out of $(nx*ny)"

    xx = float_type.(rand(1:(nx-2), nx, ny)) # Random x index for disturbance
    yy = float_type.(rand(1:(ny-2), nx, ny)) # Random y index for disturbance
    n_x = reshape(1:nx, :, 1) .* ones(1, ny)  # x index
    n_y = reshape(1:ny, 1, :) .* ones(ny, 1)  # y index
    n_x0 = copy(n_x)
    n_y0 = copy(n_y)
    
    assign_random_elements!(n_x, n_y, xx, yy, M) # Select only M links

    distx = abs.(n_x - n_x0)
    disty = abs.(n_y - n_y0)
    distx = min.(distx, nx .- distx).^2
    disty = min.(disty, ny .- disty).^2
    ww = distx .+ disty  # link weights
    ww[ww .== 0.0] .= 1.0   # avoid division by zero
    ww = 1.0 ./ (ww*dx*dy)  # weights inverted and dx*dy included in d_t

    ww = Data.Array(ww)  # Convert to a parallel data array
    n_x = Data.Array(n_x)  # Convert to a parallel data array
    n_y = Data.Array(n_y)  # Convert to a parallel data array 

    return n_x, n_y, ww
end

function assign_random_elements!(dd1, dd2, ss1, ss2, M)
  """
  Assigns M random elements of matrices ss1 and ss2 to the corresponding 
  values in matrices dd1 and dd2 respectively. 

  Args:
    ss1: The first source matrix.
    ss2: The second source matrix.
    dd1: The first destination matrix.
    dd2: The second destination matrix.
    M: The number of random elements to assign.
  """
  
  # Check if the matrices have compatible dimensions
  if  size(ss1) != size(dd1) || size(ss1) != size(ss2) || size(dd1) != size(dd2)
    error("Matrices must have the same dimensions for element-wise assignment.")
  end

  # Get the dimensions of the matrices
  rows_ss1, cols_ss1 = size(ss1)

  # Generate M unique random linear indices
  indices = randperm(rows_ss1 * cols_ss1)[1:M]

  # Assign values from ss1 and ss2 to dd1 and dd2 at the random indices
  dd1[indices] .= ss1[indices]
  dd2[indices] .= ss2[indices]

  return
end
