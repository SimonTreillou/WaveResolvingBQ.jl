using Plots
using WaveResolvingBQ
using Statistics
using BenchmarkTools

# Parameters
P=WaveResolvingBQ.SimulationParameters(Δt=0.01,T=5.0)
G=WaveResolvingBQ.Grid(Lx=6,Δx=0.01)

# Init
C_init = window(G.x,0.1,.5)
C_init = zeros(G.Nx)
for i in 1:G.Nx
    if (G.x[i]>0.5) & (G.x[i]<1.0)
        C_init[i]=1.0
    end
end

plot(G.x, C_init, title="1D Advection toy case", xlabel="x (m)", ylabel="C (ppm)", line=(:black, 4, :dotted), label="T=0")

# Model
M=WaveResolvingBQ.Advection(0.5)
B=WaveResolvingBQ.Null(2)

# Euler+UP1
T=WaveResolvingBQ.Euler()
S=WaveResolvingBQ.UP1()
C_UP1,TV=  WaveResolvingBQ.run(M,T,S,P,B,G,C_init)

# RK4+WENO5
T=WaveResolvingBQ.RK4()
S=WaveResolvingBQ.WENO5()
C_WENO,TV=  WaveResolvingBQ.run(M,T,S,P,B,G,C_init)

plot!(G.x,C_UP1,label="Euler+UP1",line=(4,:blue,0.7))
plot!(G.x,C_WENO,label="RK4+WENO5",line=(4,:red,0.7))
