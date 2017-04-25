# See https://web.stanford.edu/~lcambier/fast/tuto.php
@testset "Hydro Thermal Scheduling with $solver" for solver in lp_solvers
    num_stages = 5
    numScen = 2
    C = 5
    V = 8
    d = 6
    r = [2, 10]
    x = Matrix{JuMP.Variable}(num_stages, numScen)
    y = Matrix{JuMP.Variable}(num_stages, numScen)
    p = Matrix{JuMP.Variable}(num_stages, numScen)

    models = Vector{JuMP.Model}(num_stages)
	for s in 1:num_stages
        for ξ in 1:(s == 1 ? 1 : numScen)
			if s == 1
				model = StructuredModel(num_scenarios=numScen)
			else
				model = StructuredModel(parent=models[s-1], prob=1/2, same_children_as=(ξ == 1 ? nothing : models[s]), id=ξ, num_scenarios=(s == num_stages ? 0 : numScen))
			end
			x[s, ξ] = @variable(model, lowerbound=0, upperbound=V)
			y[s, ξ] = @variable(model, lowerbound=0)
			p[s, ξ] = @variable(model, lowerbound=0)
			if s > 1
				@constraint(model, x[s, ξ] <= x[s-1, 1] + r[ξ] - y[s, ξ])
			else
                @constraint(model, x[s, ξ] <= mean(r) - y[s, ξ])
			end
			@constraint(model, p[s, ξ] + y[s, ξ] >= d)
			@objective(model, Min, C * p[s, ξ])
			if ξ == 1
				models[s] = model
			end
		end
	end

    for cutmode in [:MultiCut, :AveragedCut]
        for detectlb in [false, true]
            lattice = model2lattice(models[1], num_stages, solver, AvgCutPruningAlgo(-1), cutmode, detectlb)
            sol = SDDP(lattice, num_stages, K = 16, stopcrit = Pereira(2, 0.5) | IterLimit(10), verbose = 0)

            @test sol.attrs[:niter] == (detectlb ? 3 : 5)
            @test sol.status == :Optimal
            @test sol.objval == (detectlb ? 15 : 18.75)

            SDDPclear(models[1])
        end
    end
end
