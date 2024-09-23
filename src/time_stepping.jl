module TimeStepping

export rk4_step

using ..Equations
using ..Discretization

function rk4_step(η, u, dt, dx, h, g)
    # k1
    dη1 = Equations.compute_deta_dt(η, u, h, dx)
    du1 = Equations.compute_du_dt(η, u, h, g, dx)
    
    # k2
    dη2 = Equations.compute_deta_dt(η .+ 0.5 * dt * dη1, u .+ 0.5 * dt * du1, h, dx)
    du2 = Equations.compute_du_dt(η .+ 0.5 * dt * dη1, u .+ 0.5 * dt * du1, h, g, dx)
    
    # k3
    dη3 = Equations.compute_deta_dt(η .+ 0.5 * dt * dη2, u .+ 0.5 * dt * du2, h, dx)
    du3 = Equations.compute_du_dt(η .+ 0.5 * dt * dη2, u .+ 0.5 * dt * du2, h, g, dx)
    
    # k4
    dη4 = Equations.compute_deta_dt(η .+ dt * dη3, u .+ dt * du3, h, dx)
    du4 = Equations.compute_du_dt(η .+ dt * dη3, u .+ dt * du3, h, g, dx)
    
    # Update
    η_new = η .+ (dt / 6) .* (dη1 + 2*dη2 + 2*dη3 + dη4)
    u_new = u .+ (dt / 6) .* (du1 + 2*du2 + 2*du3 + du4)
    
    return η_new, u_new
end

end # module