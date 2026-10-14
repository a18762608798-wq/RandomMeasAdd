using Test
using RandomMeasAdd

@testset "RandomMeasAdd" begin
    include("unit/test_robust_calibration.jl")
    include("unit/test_data_io.jl")
    include("integration/test_robust_hamming.jl")
end
