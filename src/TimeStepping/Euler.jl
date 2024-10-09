struct Euler <: TimeStepper

end

function timestepping(::Euler,P::SimulationParameters, var, rhs)
    return var .+ P.Δt .* rhs(var)
end
