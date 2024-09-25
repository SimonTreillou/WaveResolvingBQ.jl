# 1. Définir un type abstrait pour les méthodes de time-stepping
abstract type Model end

# 2. Définir des sous-types concrets pour chaque méthode
struct Advection <: Model

end

function modelize(::Advection,S::SpatialScheme,η, Δx, c, a, b)
    η_x = discretize(S,η, Δx, c, a, b)
    return  c * η_x #.+ 0.01 * η_xx
end
