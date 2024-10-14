"""
    UP1 is the upwind scheme of order 1 defined as:

        - (uᵢⁿ - uᵢⁿ⁻¹)/Δx for positive velocities
        - (uᵢⁿ⁺¹ - uᵢⁿ)/Δx for negative velocities
"""
struct UP1{T} <: SpatialScheme
    nghost::Union{T,Nothing}
end

function UP1(; nghost=nothing)
    if isnothing(nghost)
        nghost=2
    end
    return UP1(nghost)
end

function flux(S::UP1,var,i,G::Grid,order)
    if order==1
        return var[i]
    elseif order==2
        return (var[i+1]-var[i])/G.Δx
    end
end