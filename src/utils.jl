module Utils

export initialize_conditions

function initialize_conditions(x, L)
    η = exp.(-((x .- L/2) / 50.0).^2)
    u = zeros(length(x))
    return η, u
end

end # module