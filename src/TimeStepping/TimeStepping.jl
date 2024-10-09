# 1. Définir un type abstrait pour les méthodes de time-stepping
abstract type TimeStepper end

# 2. Définir des sous-types concrets pour chaque méthode
#struct Euler <: TimeStepper
    # Paramètres spécifiques à Euler si nécessaire
#end

include("Euler.jl")

#struct RK3 <: TimeStepper
    # Paramètres spécifiques à RK4 si nécessaire
#end
include("RK3.jl")

#struct RK4 <: TimeStepper
    # Paramètres spécifiques à RK4 si nécessaire
#end
include("RK4.jl")

# 3. Définir la fonction générique `timestepping` avec dispatch multiple
"""
# Méthode pour Euler
function timestepping(::Euler, var, Δt, Δx, c, a, b, rhs)
    # Étape d'Euler
    return var .+ Δt .* rhs(var,Δx,c,a,b)
end

function timestepping(::RK4, var, Δt, Δx, c, a, b, rhs)
    k1 = Δt * rhs(var, Δx, c, a, b)
    k2 = Δt * rhs(var .+ 0.5 * k1, Δx, c, a, b)
    k3 = Δt * rhs(var .+ 0.5 * k2, Δx, c, a, b)
    k4 = Δt * rhs(var .+ k3, Δx, c, a, b)
    return var .+ (k1 .+ 2*k2 .+ 2*k3 .+ k4) / 6
end
"""