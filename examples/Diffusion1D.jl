using Plots
using WaveResolvingBQ
using Statistics
using BenchmarkTools

# Parameters
P=WaveResolvingBQ.SimulationParameters(Δt=0.02,T=100.0)
G=WaveResolvingBQ.Grid(Lx=100.0,Δx=1.0)

# Initialization
μ = 50.0
D = 1.0
Mₜ = 10.0
σ = 10.0
function exact_solution(x,t)
    return exp.( -(x.-μ).^2 / (2 * σ^2 + 4 .*D .*t)) .* Mₜ * σ ./ sqrt.(σ^2 + 2*D*t)
end
C_init = exact_solution(G.x,0.0)
plot(G.x, C_init, title="1D Diffusion toy case", xlabel="x (m)", ylabel="C (ppm)", line=(:black, 4, :dotted), label="T=0")

# Model
M=WaveResolvingBQ.Diffusion(D)
B=WaveResolvingBQ.Null(2)


# Timestepping & Spatial scheme: Euler+UP1
T=WaveResolvingBQ.Euler()
S=WaveResolvingBQ.UP1()
C_UP1,TV=  WaveResolvingBQ.run(M,T,S,P,B,G,C_init)

# Timestepping & Spatial scheme: RK4+WENO5
T=WaveResolvingBQ.RK4()
S=WaveResolvingBQ.WENO5()
C_WENO,TV=  WaveResolvingBQ.run(M,T,S,P,B,G,C_init)

plot!(G.x,C_UP1,label="Euler+UP1",line=(4,:blue,0.5))
plot!(G.x,C_WENO,label="RK4+WENO5",line=(4,:green,0.5))
plot!(G.x, exact_solution(G.x,P.T), line=(2,:dash, :red),linestyle=:dot, label="Exact solution")