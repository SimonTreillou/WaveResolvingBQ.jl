
struct UP1 <: SpatialScheme
end


function discretize(::UP1,var,G::Grid,order;velocity=var*0 .+ 1.0)
    #n = length(var)
    dudx2 = zeros(G.Nx)  # To store the second derivative
    
    if order==1
        for i=2:(G.Nx-1)
            dudx2[i]=(var[i]-var[i-1])/G.Δx
        end
    elseif order==2
        for i=2:(G.Nx-1)
            dudx2[i]=(var[i+1]-2*var[i]+var[i-1])/G.Δx^2
        end
    end

    return dudx2
end
