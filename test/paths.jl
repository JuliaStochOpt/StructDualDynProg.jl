@testset "Path Sampler" begin
    @test_throws AssertionError StructDualDynProg.SDDP._samplepaths(4, [0.8, 0.3], false, false)
    @test_throws AssertionError StructDualDynProg.SDDP._samplepaths(4, [0.5, 0.3], true, false)
    pmf = [0.5, 0.5]
    @test [2, 2] == StructDualDynProg.SDDP._samplepaths(4, pmf, true, false)
    @test [0.5, 0.5] == pmf # It shouldn't have been modified
    @test [2, 3] == sort(StructDualDynProg.SDDP._samplepaths(5, pmf, true, false))
end
