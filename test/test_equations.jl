using Test
using WaveResolvingBQ

@testset "Window with Euler+UP1 at CFL=1" begin     
    # Parameters
    P=WaveResolvingBQ.SimulationParameters(Δt=0.1,T=5000.0,Δx=0.1,Lx=100.0,c=1.0)

    # Init
    η_init = window(P.x,0.1,1.0)

    # Model
    M=WaveResolvingBQ.Advection()
    # Timestepping
    T=WaveResolvingBQ.Euler()
    # Spatial scheme
    S=WaveResolvingBQ.UP1()

    η_end,TV= WaveResolvingBQ.run(M,T,S,P,η_init)

    @test η_init ≈ η_end
end