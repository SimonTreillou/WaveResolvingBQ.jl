abstract type Boundaries end


struct Periodic1D{T} <: Boundaries
    num_ghost::T
end

struct Null{T} <: Boundaries
    num_ghost::T
end

struct ReflectiveBoundary1D <: Boundaries
end

# Implémentation spécifique pour Periodic
function apply_boundaries!(B::Periodic1D,G::Grid, var)
    num = B.num_ghost
    n = length(var) - 2 * num  # Supposons que u inclut les cellules fantômes

    # Copier les valeurs du bord droit vers les cellules fantômes de gauche
    var[1:num] .= var[n+1:n+num]

    # Copier les valeurs du bord gauche vers les cellules fantômes de droite
    var[n+num+1:end] .= var[num+1:num+num]
end

function apply_boundaries!(B::Null,G::Grid, var)
    num = B.num_ghost
    n = length(var) - 2 * num  # Supposons que u inclut les cellules fantômes

    # Copier les valeurs du bord droit vers les cellules fantômes de gauche
    var[1:num] .= 0.0

    # Copier les valeurs du bord gauche vers les cellules fantômes de droite
    var[n+num+1:end] .= 0.0
end
