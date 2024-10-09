
struct WENO5 <: SpatialScheme
    # Paramètres spécifiques à RK4 si nécessaire
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

function flux_weno52(s_lst, u, i, dt, dx)
    """
    Main flux function for WENO5.
    """
    s = s_lst[end]

    if u[i] > 0
        return u[i] * weno5(s, i, 1)
    else
        return u[i] * weno5(s, i + 1, -1)
    end
end

function weno5_first_derivative(η, Δx,c)
    n = length(η)
    #d_eta_dx = Union{Float64, Nothing}[nothing for _ in 1:n]  # To store the first derivative
    d_eta_dx = zeros(n)
    M = floor(Int,n/2)

    #S=-floor.(Int,sign.(c))

    # Loop through interior points (avoiding boundaries for now)
    for i in 4:n-4
        # Compute WENO5 derivative at each point using 5-point stencil
        if c[i]>0
        #if true
        #if S==1
            S=1
            #d_eta_dx[i] =  (-weno5(η[i-3:i+1], 3, S)+weno5(η[i-2:i+2], 3, S))/Δx
            d_eta_dx[i] = (weno5(η, i, 1) - weno5(η, i-1, 1)) / Δx
            #d_eta_dx[i] =  (weno5(η[i-1:i+3], 3, S)-weno5(η[i-2:i+2], 3, S))/Δx
        elseif c[i]<0
        #elseif false
            S=-1
        #elseif S==-1
            #d_eta_dx[i] =  (weno5(η[i-1:i+3], 3, -1)-weno5(η[i-2:i+2], 3, -1))/Δx
            #d_eta_dx[i] =  (weno5(η[i-2:i+2], 3, S)-weno5(η[i-1:i+3], 3, S))/Δx
            d_eta_dx[i] = (weno5(η, i+1, -1) - weno5(η, i, -1)) / Δx
        end
    end
    d_eta_dx[3] = (η[4]-η[3])/Δx
    d_eta_dx[2] = (η[3]-η[2])/Δx
    d_eta_dx[1] = (η[2]-η[1])/Δx

    d_eta_dx[n-3] = (η[n-3]-η[n-4])/Δx
    d_eta_dx[n-2] = (η[n-2]-η[n-3])/Δx
    d_eta_dx[n-1] = (η[n-1]-η[n-2])/Δx
    d_eta_dx[n] = (η[n]-η[n-1])/Δx
    #d_eta_dx[1:3].=
   # d_eta_dx[n-3:n].=NaN
    
    """
    # Handle boundary conditions (if needed)
    d_eta_dx[1]=(weno5([η[n], η[1], η[2], η[3], η[4]], 3, S)-weno5([η[n-1], η[n], η[1], η[2], η[3]], 3, S))/Δx
    d_eta_dx[2]=(weno5([η[1], η[2], η[3], η[4], η[5]], 3, S)-weno5([η[n], η[1], η[2], η[3], η[4]], 3, S))/Δx
    d_eta_dx[3]=(weno5([η[2], η[3], η[4], η[5], η[6]], 3, S)-weno5([η[1], η[2], η[3], η[4], η[5]], 3, S))/Δx

    d_eta_dx[3]=d_eta_dx[4]
    d_eta_dx[2]=d_eta_dx[4]
    d_eta_dx[1]=d_eta_dx[4]

    d_eta_dx[n]=(weno5([η[n-1], η[n], η[1], η[2], η[3]], 3, S)-weno5([η[n-2], η[n-1], η[n], η[1], η[2]], 3, S))/Δx
    d_eta_dx[n-1]=(weno5([η[n-2], η[n-1], η[n], η[1], η[2]], 3, S)-weno5([η[n-3], η[n-2], η[n-1], η[n], η[1]], 3, S))/Δx
    d_eta_dx[n-2]=(weno5([η[n-3], η[n-2], η[n-1], η[n], η[1]], 3, S)-weno5([η[n-4], η[n-3], η[n-2], η[n-1], η[n]], 3, S))/Δx
    d_eta_dx[n-3]=(weno5([η[n-4], η[n-3], η[n-2], η[n-1], η[n]], 3, S)-weno5([η[n-5], η[n-4], η[n-3], η[n-2], η[n-1]], 3, S))/Δx

    d_eta_dx[n]=d_eta_dx[n-4]
    d_eta_dx[n-1]=d_eta_dx[n-4]
    d_eta_dx[n-2]=d_eta_dx[n-4]
    d_eta_dx[n-3]=d_eta_dx[n-4]
    """
    return d_eta_dx
end


function weno5_second_derivative3(η, Δx, c)
    n = length(η)
    d_eta_dx = zeros(n)  # To store the first derivative
    M = floor(Int,n/2)

    S=1

    # Loop through interior points (avoiding boundaries for now)
    for i in 4:n-4
        # Compute WENO5 derivative at each point using 5-point stencil
        if true
        #if S==1
            S=-1
            d_eta_dx[i] =  (WaveResolvingBQ.weno5(η[i-3:i+3], 3, S)-2*WaveResolvingBQ.weno5(η[i-2:i+2], 3, S) + WaveResolvingBQ.weno5(η[i-1:i+3], 3, S))/Δx^2
            d_eta_dx[i] =  (weno5(η[i-1:i+3], 3, S)-2*weno5(η[i-2:i+2], 3, S) + weno5(η[i-3:i+1], 3, S))/Δx^2
        elseif false
            S=-1
        #elseif S==-1
            d_eta_dx[i] =  (WaveResolvingBQ.weno5(η[i-2:i+4], 3, S)-WaveResolvingBQ.weno5(η[i-3:i+3], 3, S))/Δx
        end
    end

    d_eta_dx[3] = (η[5]-2*η[4]+η[3])/Δx^2
    d_eta_dx[2] = (η[4]-2*η[3]+η[2])/Δx^2
    d_eta_dx[1] = (η[3]-2*η[2]+η[1])/Δx^2

    d_eta_dx[n-3] = (η[n-3]-2*η[n-4]+η[n-5])/Δx^2
    d_eta_dx[n-2] =  (η[n-2]-2*η[n-3]+η[n-4])/Δx^2
    d_eta_dx[n-1] = (η[n-1]-2*η[n-2]+η[n-3])/Δx^2
    d_eta_dx[n] =  (η[n]-2*η[n-1]+η[n-2])/Δx^2
    
    #weno5_stencil([η[n-1], η[n], η[1], η[2], η[3]], Δx)
    # Handle boundary conditions (if needed)
    #d_eta_dx[1] = d_eta_dx[3]  # Example boundary handling
    #d_eta_dx[2] = d_eta_dx[3]  # Example boundary handling
    #d_eta_dx[n] = d_eta_dx[n-3]  # Example boundary handling
    #d_eta_dx[n-1] =d_eta_dx[n-3]  # Example boundary handling  
    #d_eta_dx[n-2] =d_eta_dx[n-3]  # Example boundary handling   
    # For point 1, use a second-order backward finite difference
    
    """
    d_eta_dx[1]=(weno5([η[n], η[1], η[2], η[3], η[4]], 3, S)-weno5([η[n-1], η[n], η[1], η[2], η[3]], 3, S))/Δx
    d_eta_dx[2]=(weno5([η[1], η[2], η[3], η[4], η[5]], 3, S)-weno5([η[n], η[1], η[2], η[3], η[4]], 3, S))/Δx
    d_eta_dx[3]=(weno5([η[2], η[3], η[4], η[5], η[6]], 3, S)-weno5([η[1], η[2], η[3], η[4], η[5]], 3, S))/Δx

    d_eta_dx[n]=(weno5([η[n-1], η[n], η[1], η[2], η[3]], 3, S)-weno5([η[n-2], η[n-1], η[n], η[1], η[2]], 3, S))/Δx
    d_eta_dx[n-1]=(weno5([η[n-2], η[n-1], η[n], η[1], η[2]], 3, S)-weno5([η[n-3], η[n-2], η[n-1], η[n], η[1]], 3, S))/Δx
    d_eta_dx[n-2]=(weno5([η[n-3], η[n-2], η[n-1], η[n], η[1]], 3, S)-weno5([η[n-4], η[n-3], η[n-2], η[n-1], η[n]], 3, S))/Δx
    d_eta_dx[n-3]=(weno5([η[n-4], η[n-3], η[n-2], η[n-1], η[n]], 3, S)-weno5([η[n-5], η[n-4], η[n-3], η[n-2], η[n-1]], 3, S))/Δx
    """

    """
    d_eta_dx[3]=d_eta_dx[4]
    d_eta_dx[2]=d_eta_dx[4]
    d_eta_dx[1]=d_eta_dx[4]

    d_eta_dx[n]=d_eta_dx[n-4]
    d_eta_dx[n-1]=d_eta_dx[n-4]
    d_eta_dx[n-2]=d_eta_dx[n-4]
    d_eta_dx[n-3]=d_eta_dx[n-4]
    """
    return d_eta_dx
end


function weno5_second_derivative(η, Δx,c)
    n = length(η)
    #d_eta_dx = Union{Float64, Nothing}[nothing for _ in 1:n]  # To store the first derivative
    d_eta_dx = zeros(n)
    M = floor(Int,n/2)

    #S=-floor.(Int,sign.(c))

    # Loop through interior points (avoiding boundaries for now)
    for i in 4:n-4
        # Compute WENO5 derivative at each point using 5-point stencil
        S=-1
        tmp1 =  (weno5(η, i+1, S) + weno5(η, i-1, S) - 2*weno5(η, i, S) ) / Δx^2 * 0.5
        S=1
        tmp2 = 0.5*((weno5(η, i+1, S) + weno5(η, i-1, S) - 2*weno5(η, i, S) ) / Δx^2)
        S=-1
        #d_eta_dx[i] = d_eta_dx[i]*0.5 + 0.5*((weno5(η, i+1, S) + weno5(η, i-1, S) - 2*weno5(η, i, S) ) / Δx^2)
        d_eta_dx[i] = tmp1+tmp2
        #d_eta_dx[i] =  (weno5(η[i-3:i+3], 3, S)-2*weno5(η[i-2:i+2], 3, S) + weno5(η[i-1:i+3], 3, S))/Δx^2
        #(weno5(η[i-3:i+1], 3, S)+weno5(η[i-1:i+3], 3, S)-2*weno5(η[i-2:i+2], 3, S))/Δx^2
    end


    d_eta_dx[3] = (η[4]-2*η[3]+η[2])/Δx^2
    d_eta_dx[2] = (η[3]-2*η[2]+η[1])/Δx^2
    d_eta_dx[1] = d_eta_dx[2]

    d_eta_dx[n-3] = (η[n-2]-2*η[n-3]+η[n-4])/Δx^2
    d_eta_dx[n-2] =  (η[n-1]-2*η[n-2]+η[n-3])/Δx^2
    d_eta_dx[n-1] = (η[n]-2*η[n-1]+η[n-2])/Δx^2
    d_eta_dx[n] =  d_eta_dx[n-1]
    return d_eta_dx
end


function discretize(::WENO5,var,G::Grid,order;velocity=var*0 .+ 1.0)
    
    if order==1
        tmp= weno5_first_derivative(var,G.Δx,velocity) # dépend de quel type d'advection...
    elseif order==2
        tmp= weno5_second_derivative(var,G.Δx,velocity)
        #tmp= weno5_first_derivative(tmp,G.Δx,velocity)
    end
    return tmp
end