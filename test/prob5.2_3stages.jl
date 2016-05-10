include("data.jl")

numScen = 2
M = StructuredModel(num_scenarios=numScen)

x = Vector{Vector{JuMP.Variable}}(2)
v = Vector{Vector{JuMP.Variable}}(2)
y = Vector{Matrix{JuMP.Variable}}(2)

@variable M x[1][1:n] >= 0
@variable M v[1][1:n] >= 0
@constraints M begin
  x[1] .== v[1]
end
@objective(M, Min, dot(I, v[1]))


for s in 1:numScen
  M2 = StructuredModel(parent=M, prob=p2[s])
  @variable(M2, y[1][1:n, 1:m] >= 0)
  @variable(M2, x[2][1:n] >= 0)
  @variable(M2, v[2][1:n] >= 0)
  @constraints M2 begin
    x[2] .== x[1] + v[2]
    demand[j=1:m], sum(y[1][:,j]) == D2[j,s]
    ylim[i=1:n], sum(y[1][i,:]) <= x[1][i]
  end
  @objective(M2, Min, dot(I, v[2]) + dot(C, y[1] * T))
  for S in 1:numScen
    M3 = StructuredModel(parent=M2, prob=p2[S])
    @variable(M3, y[2][1:n, 1:m] >= 0)
    @constraints M3 begin
      demand[j=1:m], sum(y[2][:,j]) == D2[j,S]
      ylim[i=1:n], sum(y[2][i,:]) <= x[2][i]
    end
    @objective(M3, Min, dot(C, y[2] * T))
  end
end

for cutmode in [:AveragedCut, :MultiCut]
  root = model2lattice(M, 3, solver, cutmode)
  sol = SDDP(root, 3, cutmode, :All)

  v11value = sol.sol[1:4]
  @show sol.status
  @show sol.objval
  @show v11value
  @test sol.status == :Optimal
  @test abs(sol.objval - 406712.49) < 0.1
  @test norm(v11value - [2986,0,7329,854]) < 0.1
  SDDPclear(M)
end
