module Equations

export compute_deta_dt, compute_du_dt

using ..Discretization

function compute_deta_dt(η, u, h, dx)
    print("lol")
    dηdx = Discretization.first_derivative(u, dx)
    du_eta_dx = Discretization.first_derivative(η, dx)
    return - ((h .+ η) .* dηdx .+ u .* du_eta_dx)
end

function compute_du_dt(η, u, h, g, dx)
    dudx = Discretization.first_derivative(u, dx)
    deta_dx = Discretization.first_derivative(η, dx)
    d3udx3 = Discretization.third_derivative(u, dx)
    return - u .* dudx - g .* deta_dx + (h^2 / 3) .* d3udx3
end

end # module