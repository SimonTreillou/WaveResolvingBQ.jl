@kwdef struct Constants
    g::Float64 = 9.81
end

function window(x,perc,A)
    n = length(x)
    tmp = zeros(n)
    imin = floor(Int,n*0.5 - n*perc*0.5)
    imax = floor(Int,n*0.5 + n*perc*0.5)
    tmp[imin:imax] .= A
    return tmp
end


