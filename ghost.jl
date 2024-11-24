# Update ghost points

@parallel_indices (ix) function update_ghost_x!(a::Data.Array)
    # Update ghost points loopin in x direction
    a[ix,end] = a[ix,2]
    a[ix,1] = a[ix,end-1]
    return
end

@parallel_indices (iy) function update_ghost_y!(a::Data.Array)
    # Update ghost points looping in y direction
    a[end,iy] = a[2,iy]
    a[1,iy] = a[end-1,iy]
    return
end

function update_ghost!(b::Data.Array, nx::Data.Index, ny::Data.Index)
    # Update ghost points

    @parallel (1:nx) update_ghost_x!(b);
    @parallel (1:ny) update_ghost_y!(b);

    return
end

function update_ghost_serial!(a::Data.Array)
    # Update ghost points
    a[end,:] = a[2,:]
    a[1,:] = a[end-1,:]
    a[:,end] = a[:,2]
    a[:,1] = a[:,end-1]
    return
end

