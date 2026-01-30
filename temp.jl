using Random
using Statistics
using StatsBase

function remove_disturbances!(n_x, n_y, nx, ny, M_r)

  #create normal indexes
  n_x0 = reshape(1:nx, :, 1) .* ones(1, ny)  # x index
  n_y0 = reshape(1:ny, 1, :) .* ones(ny, 1)  # y index

  #flatten arrays
  f_n_x=reshape(n_x, nx*ny)
  f_n_y=reshape(n_y, nx*ny)
  f_n_x0=reshape(n_x0, nx*ny)
  f_n_y0=reshape(n_y0,  nx*ny)

  #collect shortcut positions
  bool=collect(f_n_x.!=f_n_x0)
  randoms=findall(x->x==1, bool)
  
  #sample M_r shortcuts to remove
  r_1=sample(randoms, M_r, replace=false)
  #find what each shortcut is connect to 
  r_2=f_n_x[r_1].+((f_n_y[r_1].-1)*nx)

  #remove possible duplicates
  r=vcat(r_1,r_2)
  r=unique(r)

  #reset and reshape
  f_n_x[r]=f_n_x0[r]
  f_n_y[r]=f_n_y0[r]
  n_x=reshape(f_n_x, (nx,ny))
  n_y=reshape(f_n_y, (nx,ny))

  return
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
  dd1_0=copy(dd1)
  dd2_0=copy(dd2)
  # Generate M unique random linear indices
  indices = randperm(rows_ss1 * cols_ss1)[1:M]

  dd1[CartesianIndex.(Int.(ss1[indices].-1), Int.(ss2[indices].-1))] .= dd1_0[indices]
  dd1[indices] .= ss1[indices]
  f=dd2[CartesianIndex.(Int.(ss1[indices]), Int.(ss2[indices]))]
  dd2[CartesianIndex.(Int.(ss1[indices].-1), Int.(ss2[indices].-1))] .= dd2_0[indices]
  dd2[indices] .= ss2[indices]

  return
end

