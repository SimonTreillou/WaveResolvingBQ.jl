struct SimulationParameters
    # Time
    Δt::Union{Float64, Nothing}         # Timestep
    T::Union{Float64, Nothing}         # Simulated time
    Nt::Union{Integer, Nothing}            # Number of timesteps
    t::Union{Vector{Float64}, Nothing}
end

function SimulationParameters(; Δt=nothing, Nt=nothing,T=1000.0, t=nothing)
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


    return SimulationParameters(Δt,T,Nt,t)
end


function check_parameters(P::SimulationParameters)
    println("Your CFL is "*string(P.c*P.Δt/P.Δx))
end