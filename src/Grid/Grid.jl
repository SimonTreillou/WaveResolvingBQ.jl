struct Grid
    # Space
    Δx::Union{Float64, Nothing}            # Resolution
    Lx::Union{Float64, Nothing}            # Domain length
    Nx::Union{Integer, Nothing}            # Number of grid points
    x::Union{Vector{Float64}, Nothing}

    # Depth
    # ....
end

function Grid(;Δx=nothing,Lx=10.0,Nx=nothing,x=nothing)
    # Space
    if Δx !== nothing
        Nx = floor(Int, Lx / Δx)
        Lx = Nx * Δx
        x = collect(0:Δx:Lx-Δx)
    elseif Nx !== nothing
        Δx = Lx / Nx
        x = collect(0:Δx:Lx-Δx)
    else
        error("Vous devez fournir au moins une résolution (Δx) ou un nombre de points de grille (Nx) !")
    end
    return Grid(Δx,Lx,Nx,x)
end