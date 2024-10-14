"""
    SpatialScheme is an abstract type defining the numerical scheme used for spatial discretization.

    Currently implemented:
        - UP1
        - WENO5
"""
abstract type SpatialScheme end


# Include definitions for flux depending on chosen spatial schemes.
include("UP1.jl")
include("WENO5.jl")


# Define main discretize function
function discretize(S::SpatialScheme,var,G::Grid,order,velocity::Float64)
    tmp = zeros(G.Nx)
    s = floor(Int,sign(velocity)<0) * floor(Int,order!==2) # 1 if velocity <0, 0 elsewhere

    for i=S.nghost:(G.Nx-S.nghost)
        tmp[i]=(flux(S,var,i+s,G,order) - flux(S,var,i-1+s,G,order))/G.Δx
    end

    return tmp
end

function discretize(S::SpatialScheme,var,G::Grid,order,velocity::Vector{Float64})
    tmp = zeros(G.Nx)
    s = floor.(Int,sign.(velocity).<0) .* floor(Int,order!==2)

    for i=S.nghost:(G.Nx-S.nghost)
        tmp[i]=(flux(S,var,i+s[i],G,order) - flux(S,var,i-1+s[i],G,order))/G.Δx
    end

    return tmp
end