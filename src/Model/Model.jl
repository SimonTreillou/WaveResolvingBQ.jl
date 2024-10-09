# 1. Définir un type abstrait pour les méthodes de time-stepping
abstract type Model end

# 2. Définir des sous-types concrets pour chaque méthode
struct Advection{T} <: Model
    c::T
end

struct Diffusion{T} <: Model
    ν::T
end

struct Burgers{T} <: Model
    ν::T
end



function modelize(M::Advection,S::SpatialScheme,G::Grid,B::Boundaries,η)
    function rhs(η)
        η_x = discretize(S,η,G,1,velocity=ones(G.Nx)*M.c)
        tmp = - M.c .* η_x
        apply_boundaries!(B,G,tmp)
        return tmp
    end
    return rhs
end



function modelize(M::Burgers,S::SpatialScheme,G::Grid,B::Periodic1D,η)
    function rhs(η)
        η_x = discretize(S,η,G,1,velocity=η)
        η_xx =  discretize(S,η,G,2,velocity=η)
       # tmp = η .* η_x - M.ν * η_xx
        tmp = -η .* η_x +  M.ν * η_xx
        apply_boundaries!(B,G,tmp)
        return tmp
    end
    return rhs
end

function modelize(M::Diffusion,S::SpatialScheme,G::Grid,B::Boundaries,η)
    function rhs(η)
        η_xx = discretize(S,η,G,2,velocity=η)
        tmp = M.ν * η_xx
        apply_boundaries!(B,G,tmp)
        return tmp
    end
    return rhs
end
