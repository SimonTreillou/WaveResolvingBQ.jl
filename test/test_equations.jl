using Test
using WaveResolvingBQ

@testset "Equations Tests" begin
    # Simple test for first derivative
    f = [0.0, 1.0, 2.0, 3.0, 4.0]
    dx = 1.0
    df_expected = [1.0, 1.0, 1.0, 1.0, 1.0]
    df_computed = BoussinesqModel.first_derivative(f, dx)
    @test df_computed ≈ df_expected
    
    # Test third derivative (example with zeros)
    d3f_expected = [0.0, 0.0, 0.0, 0.0, 0.0]
    d3f_computed = BoussinesqModel.third_derivative(f, dx)
    @test d3f_computed ≈ d3f_expected
end