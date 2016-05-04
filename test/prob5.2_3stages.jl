include("data.jl")

numScen = 2
M = StructuredModel(num_scenarios=numScen)

@variable M x1[1:n] >= 0
@variable M v1[1:n] >= 0
@constraints M begin
  x1 .== v1
end
@objective(M, Min, dot(I, v1))

for s in 1:numScen
    M2 = StructuredModel(parent=M, prob=p2[s])
    @variable(M2, y2[1:n, 1:m] >= 0)
    @variable(M2, x2[1:n] >= 0)
    @variable(M2, v2[1:n] >= 0)
    @constraints M2 begin
      x2 .== x1 + v2
      demand[j=1:m], sum(y2[:,j]) == D2[j,s]
      ylim[i=1:n], sum(y2[i,:]) <= x1[i]
    end
    @objective(M2, Min, dot(I, v2) + dot(C, y2 * T))
    for S in 1:numScen
      M3 = StructuredModel(parent=M2, prob=p2[S])
      @variable(M3, y3[1:n, 1:m] >= 0)
      @constraints M3 begin
        demand[j=1:m], sum(y3[:,j]) == D2[j,s]
        ylim[i=1:n], sum(y3[i,:]) <= x2[i]
      end
      @objective(M3, Min, dot(C, y3 * T))
    end
end

for cutmode in [:AveragedCut, :MultiCut]
  status, objval, sol = SDDP(M, 3, solver, cutmode)

  v11value = sol[1:4]
  @show status
  @show objval
  @show v11value
  @test status == :Optimal
  @test abs(objval - 405969.63) < 0.1
  @test norm(v11value - [2986,0,7329,854]) < 0.1
  SDDPclear(M)
end
