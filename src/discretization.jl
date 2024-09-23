module Discretization

export first_derivative, third_derivative

function first_derivative(f, dx)
    df = zeros(size(f))
    df[2:end-1] = (f[3:end] - f[1:end-2]) / (2 * dx)
    df[1] = (f[2] - f[1]) / dx
    df[end] = (f[end] - f[end-1]) / dx
    return df
end

function third_derivative(f, dx)
    d3f = zeros(size(f))
    d3f[3:end-2] = (f[5:end] - 2f[4:end-1] + 2f[2:end-3] - f[1:end-4]) / (2 * dx^3)
    d3f[1:2] .= 0.0
    d3f[end-1:end] .= 0.0
    return d3f
end

end # module