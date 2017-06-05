using ECOS
using JuMP
using StructJuMP
# I have started a package for SDDP
# The solution of problem 5.2 and 5.3 can be seen as a particular case of this work
using StructDualDynProg

# Telling ECOS solver to shut up :)
solver = ECOS.ECOSSolver(verbose=false)
# You can replace it with another solver if you wish

include("data.jl")

# Modelisation with the package StructJuMP
numScen = 2
M1 = StructuredModel(num_scenarios=numScen)

@variable M1 x1[1:n] >= 0
# change the false to true to see the effect of adding the variable v1
if false
    @variable M1 v1[1:n] >= 0
    @constraints M1 begin
        x1 .== v1
    end
    @objective(M1, Min, dot(ic, v1))
else
    @objective(M1, Min, dot(ic, x1))
end

for s in 1:numScen
    M2 = StructuredModel(parent=M1, prob=p2[s])
    @variable(M2, y2[1:n, 1:m] >= 0)
    @constraints M2 begin
        demand[j=1:m], sum(y2[:,j]) == D2[j,s]
        ylim[i=1:n], sum(y2[i,:]) <= x1[i]
    end
    @objective(M2, Min, dot(C, y2 * T))
end

# Solving with my package SDDP

for cutmode in [:AveragedCut, :MultiCut]
    status, objval, sol = SDDP(M1, 2, solver, cutmode)

    v11value = sol[1:4]
    @show status
    @show objval
    @show v11value
    SDDPclear(M1)
end
