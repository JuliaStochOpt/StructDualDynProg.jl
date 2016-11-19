@testset "2 stages" begin
    include("data.jl")

    numScen = 2
    M1 = StructuredModel(num_scenarios=numScen)

    @variable M1 x1[1:n] >= 0
    @variable M1 v1[1:n] >= 0
    @constraints M1 begin
      x1 .== v1
    end
    @objective(M1, Min, dot(I, v1))

    for s in 1:numScen
        M2 = StructuredModel(parent=M1, prob=p2[s])
        @variable(M2, y2[1:n, 1:m] >= 0)
        @constraints M2 begin
          demand[j=1:m], sum(y2[:,j]) == D2[j,s]
          ylim[i=1:n], sum(y2[i,:]) <= x1[i]
        end
        @objective(M2, Min, dot(C, y2 * T))
    end

    fulltest(M1, 2, 340315.52, [5085,1311,3919,854], 334687.754566,15869.996575)
    # root = model2lattice(M1, 2, solver, cutmode)
    # sol = SDDP(root, 2, cutmode, :All, verbose)
    #
    # v11value = sol.sol[1:4]
    # @test sol.status == :Optimal
    # @test abs(sol.objval - 340315.52) < 0.1
    # @test norm(v11value - [5085,1311,3919,854]) < 0.1
    # SDDPclear(M1)
end
