# Plotting functions

function plotbwh(b, w, nx, ny, ttot)
    l = @layout [a b]	
    h1 = heatmap(b,  aspect_ratio=:equal, xlims = (0, nx), ylims = (0, ny), title=@sprintf("b - t=%3.2f", ttot) )
    h2 = heatmap(w,  aspect_ratio=:equal, xlims = (0, nx), ylims = (0, ny), title=@sprintf("w - t=%3.2f", ttot)  )
    display(plot(h1, h2, layout = l))
end

function plotb(b, t, nx, ny, lx, ly)
    X, Y      = -lx/2:lx/(nx-2):lx/2, -ly/2:ly/(ny-2):ly/2
    heatmap(X, Y, Array(b)[2:(end), 2:(end)]', aspect_ratio=1, clim=(0, 0.65), xlims=(X[1],X[end]), ylims=(Y[1],Y[end]), c=:viridis, title="b field - time = $t")
end

