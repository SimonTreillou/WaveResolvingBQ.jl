using Test
using WaveResolvingBQ
using Documenter

@testset "Testing WaveResolvingBQ" begin

    @testset "Window with Euler+UP1 at CFL=1" begin
        include("test_equations.jl")
    end
end