@testset "3 stages with $solver" for solver in lp_solvers
    include("prob5.2_data.jl")

    numScen = 2
    M = StructuredModel(num_scenarios=numScen)

    x = Vector{Vector{JuMP.Variable}}(undef, 2)
    v = Vector{Vector{JuMP.Variable}}(undef, 2)
    y = Vector{Matrix{JuMP.Variable}}(undef, 2)

    x[1] = @variable(M, [1:n], lowerbound=0)
    v[1] = @variable(M, [1:n], lowerbound=0)
    @constraints M begin
        x[1] .== v[1]
    end
    @objective(M, Min, dot(ic, v[1]))


    for s in 1:numScen
        M2 = StructuredModel(parent=M, prob=p2[s], id=s)
        y[1] = @variable(M2, [1:n, 1:m], lowerbound=0)
        x[2] = @variable(M2, [1:n], lowerbound=0)
        v[2] = @variable(M2, [1:n], lowerbound=0)
        @constraints M2 begin
            x[2] .== x[1] + v[2]
            demand[j=1:m], sum(y[1][:,j]) == D2[j,s]
            ylim[i=1:n], sum(y[1][i,:]) <= x[1][i]
        end
        @objective(M2, Min, dot(ic, v[2]) + dot(C, y[1] * T))
        for S in 1:numScen
            M3 = StructuredModel(parent=M2, prob=p2[S], id=S)
            y[2] = @variable(M3, [1:n, 1:m], lowerbound=0)
            @constraints M3 begin
                demand[j=1:m], sum(y[2][:,j]) == D2[j,S]
                ylim[i=1:n], sum(y[2][i,:]) <= x[2][i]
            end
            @objective(M3, Min, dot(C, y[2] * T))
        end
    end

    if !isdefined(@__MODULE__, :testniter_3stages)
        include("prob5.2_3stages_niter.jl")
    end
    fulltest(M, 3, 406712.49, [2986,0,7329,854], 402593.71614, 17319.095064, testniter_3stages, solver)
    # root = model2lattice(M, 3, solver, cutmode)
    # sol = SDDP(root, 3, cutmode, :All, verbose)
    #
    # v11value = sol.sol[1:4]
    # @test sol.status == :Optimal
    # @test abs(sol.objval - 406712.49) < 0.1
    # @test norm(v11value - [2986,0,7329,854]) < 0.1
    # SDDPclear(M)
end
