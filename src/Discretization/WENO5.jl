
struct WENO5{T} <: SpatialScheme
    nghost::Union{T,Nothing}
end

function WENO5(;nghost=nothing)
    if isnothing(nghost)
        nghost=4
    end
    return WENO5(nghost)
end

function weno5(s, i, sign)
    # Constants
    k1 = 13. / 12.
    k2 = 1. / 4.
    eps = 1e-6

    g1 = 0.1
    g2 = 0.6
    g3 = 0.3

    # Get values to use, depending on sign
    s_mm = s[i - sign * 2]
    s_m = s[i - sign * 1]
    s_0 = s[i]
    s_p = s[i + sign * 1]
    s_pp = s[i + sign * 2]

    # Compute stencils
    qi1 = 1. / 3. * s_mm - 7. / 6. * s_m + 11. / 6. * s_0
    qi2 = -1. / 6. * s_m + 5. / 6. * s_0 + 1. / 3. * s_p
    qi3 = 1. / 3. * s_0 + 5. / 6. * s_p - 1. / 6. * s_pp

    # Compute indicators of smoothness
    beta1 = k1 * (s_mm - 2 * s_m + s_0) ^ 2 + k2 * (s_mm - 4 * s_m + 3 * s_0) ^ 2
    beta2 = k1 * (s_m - 2 * s_0 + s_p) ^ 2 + k2 * (s_m - s_p) ^ 2
    beta3 = k1 * (s_0 - 2 * s_p + s_pp) ^ 2 + k2 * (3 * s_0 - 4 * s_p + s_pp) ^ 2

    # Compute weights
    w1 = g1 / (beta1 + eps) ^ 2
    w2 = g2 / (beta2 + eps) ^ 2
    w3 = g3 / (beta3 + eps) ^ 2

    return (w1 * qi1 + w2 * qi2 + w3 * qi3) / (w1 + w2 + w3)
end


function flux_weno5(s_lst, u, i, dt, dx)
    """
    Main flux function for WENO5.
    """
    #s = s_lst[end]
    s = s_lst

    if u[i] > 0
        return u[i] * weno5(s, i, 1)
    else
        return u[i] * weno5(s, i + 1, -1)
    end
end

function flux(::WENO5,var,i,G::Grid,order)
    if order==1
        return weno5(var, i, 1)
    elseif order==2
        return (weno5(var, i+1, 1)-weno5(var, i, 1))/G.Δx
    end
end

function weno5_first_derivative(var,G::Grid,c)
    dvar = zeros(G.Nx)

    # Loop through interior points (avoiding boundaries for now)
    for i in 4:G.Nx-4
        if c[i]>0
            dvar[i] = (weno5(var, i, 1) - weno5(var, i-1, 1)) / G.Δx
        elseif c[i]<0
            dvar[i] = (weno5(var, i+1, -1) - weno5(var, i, -1)) / G.Δx
        end
    end
    dvar[3] = (var[4]-var[3])/G.Δx
    dvar[2] = (var[3]-var[2])/G.Δx
    dvar[1] = (var[2]-var[1])/G.Δx

    dvar[end-3] = (var[end-3]-var[end-4])/G.Δx
    dvar[end-2] = (var[end-2]-var[end-3])/G.Δx
    dvar[end-1] = (var[end-1]-var[end-2])/G.Δx
    dvar[end] = (var[end]-var[end-1])/G.Δx

    return dvar
end


function weno5_second_derivative(var,G::Grid,c)
    d2var = zeros(G.Nx)

    # Loop through interior points (avoiding boundaries for now)
    for i in 4:G.Nx-4
        S=-1
        tmp1 =  (weno5(var, i+1, S) + weno5(var, i-1, S) - 2*weno5(var, i, S) ) / G.Δx^2 * 0.5
        S=1
        tmp2 = 0.5*((weno5(var, i+1, S) + weno5(var, i-1, S) - 2*weno5(var, i, S) ) / G.Δx^2)
        S=-1
        d2var[i] = tmp1+tmp2
    end


    d2var[3] = (var[4]-2*var[3]+var[2])/G.Δx^2
    d2var[2] = (var[3]-2*var[2]+var[1])/G.Δx^2
    d2var[1] = d2var[2]

    d2var[end-3] = (var[end-2]-2*var[end-3]+var[end-4])/G.Δx^2
    d2var[end-2] =  (var[end-1]-2*var[end-2]+var[end-3])/G.Δx^2
    d2var[end-1] = (var[end]-2*var[end-1]+var[end-2])/G.Δx^2
    d2var[end] =  d2var[end-1]
    return d2var
end


"""
function discretize(::WENO5,var,G::Grid,order;velocity=var*0 .+ 1.0)
    if order==1
        tmp= weno5_first_derivative(var,G,velocity)
    elseif order==2
        tmp= weno5_second_derivative(var,G,velocity)
    else
        throw("Order isn't correct. It must be 1 or 2, depending on chosen derivative.")
    end
    return tmp
end
"""