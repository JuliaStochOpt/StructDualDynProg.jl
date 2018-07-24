@testset "Comparison" begin
    for mod in (StructProg, SDDP)
        @test mod._isapproxzero(0, 1e-10)
        @test !mod._isapproxzero(1, 1e-10)
        @test mod._isapproxzero(1e-5, 1e-4)
        @test !mod._isapproxzero(1e-3, 1e-4)
        @test mod._isapproxzero([1e-5, 1e-5], 1e-4)
        @test !mod._isapproxzero([1e-5, 1e-3], 1e-4)
        @test mod._isapprox(0, 0, 1e-3)
        @test mod._isapprox(1, 1, 1e-3)
        @test !mod._isapprox(0, 1, 1e-3)
        @test mod._lt(0, 1, 1e-3)
        @test !mod._lt(1, 1, 1e-3)
        @test !mod._lt(1, 0, 1e-3)
        @test !mod._lt(1, BigInt(1), 1e-3)
    end
end
