using Test
using WaveResolvingBQ

@testset "Window with Euler+UP1 at CFL=1" begin     
    # Parameters
    P=WaveResolvingBQ.SimulationParameters(Δt=0.01,T=4.0)
    G=WaveResolvingBQ.Grid(Lx=6,Δx=0.01)

    # Init
    C_init = zeros(G.Nx)
    for i in 1:G.Nx
        if (G.x[i]>0.5) & (G.x[i]<1.0)
            C_init[i]=1.0
        end
    end

    # Model
    M=WaveResolvingBQ.Advection(1.0)
    B=WaveResolvingBQ.Null(2)

    # Euler+UP1
    T=WaveResolvingBQ.Euler()
    S=WaveResolvingBQ.UP1()
    C_UP1,TV=  WaveResolvingBQ.run(M,T,S,P,B,G,C_init)

    # Solution
    C_sol = zeros(G.Nx)
    for i in 1:G.Nx
        if (G.x[i]>4.5) & (G.x[i]<5.0)
            C_sol[i]=1.0
        end
    end

    @test C_sol ≈ C_UP1
end