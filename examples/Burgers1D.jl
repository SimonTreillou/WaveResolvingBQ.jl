using Plots
using WaveResolvingBQ
using Statistics
using BenchmarkTools
using LaTeXStrings

# Parameters
ν = 0.07
P=WaveResolvingBQ.SimulationParameters(Δt=2*π/400 * ν,Nt=100,T=100*2*π/100 * ν)
G=WaveResolvingBQ.Grid(Lx=2*π,Δx=2*π/400,Nx=101)

#check_parameters(P)

# Init: initial function and its analytical solution
function Φ(x,t)
    return exp.(-(x .- 4*t).^2 / (4*ν*(t+1))) + exp.(-(x .- 4*t .- 2*π).^2 / (4*ν*(t+1)))
end

function dΦdx(x,t)
    term1 = -((x .- 4t)) ./ (2 * ν * (t + 1)) .* exp.(-((x .- 4t).^2) / (4 * ν * (t + 1)))
    term2 = -((x .- 4t .- 2π)) ./ (2 * ν * (t + 1)) .* exp.(-((x .- 4t .- 2π).^2) / (4 * ν * (t + 1)))
    return term1 .+ term2
end

u_init = - 2*ν ./Φ(G.x,0.0) .* dΦdx(G.x,0.0) .+ 4
exact_solution = - 2*ν ./Φ(G.x,P.T) .* dΦdx(G.x,P.T) .+ 4
plot(G.x, u_init, title=L"\textbf{Burgers \enspace equation:} \partial_t u + u \partial_x u = \nu \partial_{xx} u", xlabel="x (m)", ylabel="u (ppm)", line=(:black, 4), label="T=0")

# Model
M=WaveResolvingBQ.Burgers(ν)
# Boundaries
B=WaveResolvingBQ.Periodic1D(4)

# Euler + UP1
T=WaveResolvingBQ.Euler()
S=WaveResolvingBQ.UP1()
uUP1,TV=  WaveResolvingBQ.run(M,T,S,P,B,G,u_init)

# RK4 + WENO5
T=WaveResolvingBQ.RK4()
S=WaveResolvingBQ.WENO5()
uWENO,TV=  WaveResolvingBQ.run(M,T,S,P,B,G,u_init)


plot!(G.x,uUP1,label="Euler+UP1",line=(4,:blue,0.5))
plot!(G.x,uWENO,label="RK4+WENO5",line=(4,:green,0.5))
plot!(G.x, exact_solution, line=(2,:dash, :red),linestyle=:dot, label="Exact solution")