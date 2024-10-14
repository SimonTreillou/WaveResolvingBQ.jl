struct Grid
    # Space
    Δx::Union{Float64, Nothing}            # Resolution
    Lx::Union{Float64, Nothing}            # Domain length
    Nx::Union{Integer, Nothing}            # Number of grid points
    x::Union{Vector{Float64}, Nothing}

    # Depth
    h::Union{Vector{Float64}, Float64, Nothing}
end

function Grid(;Δx=nothing,Lx=nothing,Nx=nothing,x=nothing,h=nothing)
    if isnothing(Lx)
        throw("You need to have at least the domain size (wtf is wrong with you?)")
    end

    # Space
    if isnothing(Nx) && !isnothing(Δx)
        x = collect(0:Δx:Lx-Δx)
        Nx = length(x)
    elseif isnothing(Δx) && !isnothing(Nx)
        Δx = Lx / Nx
        x = collect(0:Δx:Lx-Δx)
    elseif isnothing(Δx) && isnothing(Nx)
        throw("Vous devez fournir au moins une résolution (Δx) ou un nombre de points de grille (Nx) !")
    else
        x = collect(0:Δx:Lx-Δx)
    end

    # Depth 
    if isnothing(h)
        h=ones(Nx)
    elseif size(h)==()
        h=ones(Nx).*h
    end
    return Grid(Δx,Lx,Nx,x,h)
end