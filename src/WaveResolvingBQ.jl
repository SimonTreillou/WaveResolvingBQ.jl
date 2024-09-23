module WaveResolvingBQ

# Include submodule files
include("equations.jl")
include("discretization.jl")
include("time_stepping.jl")
include("boundary_conditions.jl")
include("utils.jl")

export run_simulation, BoussinesqState

using .Utils
using .Equations
using .Discretization
using .TimeStepping
using .BoundaryConditions


# Define a structure to hold the simulation state
struct BoussinesqState
    x::Vector{Float64}
    η::Vector{Float64}
    u::Vector{Float64}
    t::Float64
end

# Main simulation function
function run_simulation(L::Float64 = 1000.0, Nx::Int = 1000, h::Float64 = 10.0, g::Float64 = 9.81,
                       T::Float64 = 10.0, dt::Float64 = 0.01)
    dx = L / (Nx - 1)
    x = LinRange(0, L, Nx)
    η, u = Utils.initialize_conditions(x, L)
    
    Nt = Int(T / dt)
    
    state = BoussinesqState(x, η, u, 0.0)
    
    for n in 1:Nt
        state.η, state.u = TimeStepping.rk4_step(state.η, state.u, dt, dx, h, g)
        BoundaryConditions.apply_boundary_conditions!(state.η, state.u)
        state.t += dt
        
        # Optional: Print progress or handle output
        if n % (Nt ÷ 10) == 0
            @info "Simulation progress: $(100 * n / Nt)%"
        end
    end
    
    return state
end


end # module WaveResolvingBQ

