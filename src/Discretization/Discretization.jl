# 1. Définir un type abstrait pour les méthodes de time-stepping
abstract type SpatialScheme end

# 2. Définir des sous-types concrets pour chaque méthode
struct UP1 <: SpatialScheme
    # Paramètres spécifiques à Euler si nécessaire
end

struct WENO5 <: SpatialScheme
    # Paramètres spécifiques à RK4 si nécessaire
end


function discretize(::UP1,var, Δx, c, a, b)
    n = length(var)
    dudx2 = zeros(n)  # To store the second derivative
    
    dudx2[3:n-2] = (circshift(var[3:n-2], -1) - 2 * var[3:n-2] + circshift(var[3:n-2], 1)) / Δx^2
    
    # Handling boundary conditions (you may need to adjust this based on your problem)
    dudx2[1] = (var[1] - 2*var[1] + var[2]) / Δx^2  # Example boundary handling
    dudx2[2] = (var[2] - 2*var[1] + var[2]) / Δx^2  # Example boundary handling
    dudx2[n-1] = (var[n] - 2*var[n-1] + var[n]) / Δx^2  # Example boundary handling
    dudx2[n] = (var[n] - 2*var[n-1] + var[n]) / Δx^2  # Example boundary handling

    dudx2 = (circshift(var, 0) - circshift(var, 1)) / (Δx)
    return dudx2
end

include("WENO5.jl")