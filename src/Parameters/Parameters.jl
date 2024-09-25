struct SimulationParameters
    # Time
    Δt::Union{Float64, Nothing}         # Timestep
    T::Union{Float64, Nothing}         # Simulated time
    Nt::Union{Integer, Nothing}            # Number of timesteps
    t::Union{Vector{Float64}, Nothing}

    # Space
    Δx::Union{Float64, Nothing}        # Resolution
    Lx::Union{Float64, Nothing}        # Domain length
    Nx::Union{Integer, Nothing}            # Number of grid points
    x::Union{Vector{Float64}, Nothing}

    # Velocity
    c::Union{Float64, Nothing}         # Vitesse d'advection ou autre coefficient
end

function SimulationParameters(; Δt=nothing, Nt=nothing, Δx=nothing, Nx=nothing, Lx=100.0, x=nothing, T=1000.0, t=nothing, c=0.1)
    # Time
    if Δt !== nothing
        Nt = floor(Int, T / Δt)
        T = Nt * Δt
        t = collect(0:Δt:T-Δt)
    elseif Nt !== nothing
        Δt = T / Nt
        t = collect(0:Δt:T-Δt)
    else
        error("Vous devez fournir au moins un pas de temps (Δt) ou un nombre de pas de temps (Nt) !")
    end

    # Space
    if Δx !== nothing
        Nx = floor(Int, Lx / Δx)
        Lx = Nx * Δx
        x = collect(0:Δx:Lx-Δx)
    elseif Nx !== nothing
        Δx = Lx / Nx
        x = collect(0:Δx:Lx-Δx)
    else
        error("Vous devez fournir au moins une résolution (Δx) ou un nombre de points de grille (Nx) !")
    end

    return SimulationParameters(Δt,T,Nt,t, Δx,Lx,Nx,x, c)
end


function check_parameters(P::SimulationParameters)
    println("Your CFL is "*string(P.c*P.Δt/P.Δx))
end