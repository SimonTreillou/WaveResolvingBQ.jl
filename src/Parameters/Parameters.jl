@kwdef struct SimulationParameters
    # Time
    Δt::Union{Float64, Nothing} = nothing         # Timestep
    T::Union{Float64, Nothing} = nothing          # Simulated time
    Nt::Union{Int, Nothing} = nothing             # Number of timesteps
    t::Union{Vector, Nothing} = nothing

    # Space
    Δx::Union{Float64, Nothing} = nothing         # Resolution
    Lx::Union{Float64, Nothing} = nothing         # Domain length
    Nx::Union{Int, Nothing} = nothing             # Number of grid points
    x::Union{Vector, Nothing} = nothing

    # Velocity
    c::Union{Float64, Nothing} = nothing          # Vitesse d'advection ou autre coefficient
end


function SimulationParameters(; Δt=nothing, Nt=nothing, Δx=nothing, Nx=nothing, Lx=100.0, x=nothing, T=1000.0, t=nothing, c=0.1)
    #param_int = promote(aₛₕ, aₛᵥ, maxiter)
    #param_float = promote(Δt,Δx,T,c)

    # Time
    if Δt !== nothing
        Nt = floor(Int,T/Δt)
        T = Nt*Δt
        t = collect(0:Δt:T-Δt)
    elseif Nt !== nothing
        Δt = T/Nt
        t = collect(0:Δt:T-Δt)
    else
        println("You must provide at least a time step or number of time steps!")
    end

    # Space
    if Δx !== nothing
        Nx = floor(Int,Lx/Δx)
        Lx = Nx*Δx
        x = collect(0:Δx:Lx-Δx)
    elseif Nx !== nothing
        Δx = Lx/Nx
        x = collect(0:Δx:Lx-Δx)
    else
        println("You must provide at least a resolution or number of grid points!")
    end

    SimulationParameters(Δt,T,Nt,t, Δx,Lx,Nx,x, c)
end


function check_parameters(P::SimulationParameters)
    println("Your CFL is "*string(P.c*P.Δt/P.Δx))
end