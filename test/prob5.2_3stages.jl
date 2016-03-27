include("data.jl")

numScen = 2
M = StructuredModel(solver=GurobiSolver(), num_scenarios=numScen)

@defVar M x1[1:n] >= 0
@defVar M v1[1:n] >= 0
@addConstraints M begin
  x1 .== v1
end
@setObjective(M, Min, dot(I, v1))

for s in 1:numScen
    M2 = StructuredModel(parent=M, prob=p2[s])
    @defVar(M2, y2[1:n, 1:m] >= 0)
    @defVar(M2, x2[1:n] >= 0)
    @defVar(M2, v2[1:n] >= 0)
    @addConstraints M2 begin
      x2 .== x1 + v2
      demand[j=1:m], sum(y2[:,j]) == D2[j,s]
      ylim[i=1:n], sum(y2[i,:]) <= x1[i]
    end
    @setObjective(M2, Min, dot(I, v2) + dot(C, y2 * T))
    for S in 1:numScen
      M3 = StructuredModel(parent=M2, prob=p2[S])
      @defVar(M3, y3[1:n, 1:m] >= 0)
      @addConstraints M3 begin
        demand[j=1:m], sum(y3[:,j]) == D2[j,s]
        ylim[i=1:n], sum(y3[i,:]) <= x2[i]
      end
      @setObjective(M3, Min, dot(C, y3 * T))
    end
end

for cutmode in [:AveragedCut, :MultiCut]
  status, objval, sol = SDDP(M, 3, socp_solver, cutmode)

  v11value = sol[1:4]
  @show status
  @show objval
  @show v11value
  @test status == :Optimal
  @test abs(objval - 340315.52) < 0.1
  @test norm(v11value - [5085,1311,3919,854]) < 0.1
  SDDPclear(M)
end
