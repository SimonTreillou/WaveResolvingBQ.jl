struct RK4 <: TimeStepper

end



function timestepping(::RK4, P::SimulationParameters, var, rhs)
    k1 = P.Δt * rhs(var)
    k2 = P.Δt * rhs(var .+ 0.5 * k1)
    k3 = P.Δt * rhs(var .+ 0.5 * k2)
    k4 = P.Δt * rhs(var .+ k3)
    return var .+ (k1 .+ 2*k2 .+ 2*k3 .+ k4) / 6
end