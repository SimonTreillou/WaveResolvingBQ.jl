struct RK3 <: TimeStepper

end


function timestepping(::RK3, var, Δt, Δx, c, a, b, rhs)
    k1 = Δt * rhs(var, Δx, c, a, b)
    k2 = Δt * rhs(var .+ 0.5 * k1, Δx, c, a, b)
    k3 = Δt * rhs(var .- k1 .+ 2*k2, Δx, c, a, b)
    return var .+ (k1 .+ 4*k2 .+ k3) / 6
end
