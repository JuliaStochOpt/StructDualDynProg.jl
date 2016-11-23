#================================================================
This package solves the following problem parallely

min  c_0^T x + \sum{i = 1}^S c_i^Ty_i
s.t. b_0 - A_0 x           \in K_0,
b_i - A_i x - B_i y_i \in K_i, \forall i = 1,...,S
x           \in C_0,
x_I         \in Z
y_i \in C_i, \forall i = 1,...,S

where input to the Benders engine is:
c_all is an array of c_0, c_1, c_2, ... , c_S objective coefficients
A_all is an array of A_0, A_1, A_2, ... , A_S constraint matrices of master variables
B_all is an array of      B_1, B_2, ... , B_S constraint matrices of scenario variables (note that scenario variables do not appear at master constraints)
    b_all is an array of b_0, b_1, b_2, ... , b_S right hand side of constraints
    K_all is an array of K_0, K_1, K_2, ... , K_S constraint cones
    C_all is an array of C_0, C_1, C_2, ... , C_S variable cones
    v     is an array of master variable types (i.e. :Bin, :Int, :Cont)

    call with: julia -p <num_threads> <myscript>
    ================================================================#

    import MathProgBase

    # this function loads and solves a conic problem and returns its dual
    function loadAndSolveConicProblem(c, A, b, K, C, solver)

        # load conic model
        model = MathProgBase.ConicModel(solver)
        MathProgBase.loadproblem!(model, c, A, b, K, C)

        #println("process id $(myid()) started")

        # solve conic model
        MathProgBase.optimize!(model)
        status = MathProgBase.status(model)

        # return status and dual
        #println("process id $(myid()) status $(status)")
        return status, MathProgBase.getobjval(model), MathProgBase.getsolution(model), MathProgBase.getdual(model)
    end

    function addCuttingPlanes(prob, A_all, b_all, output, θ, separator, cutmode, TOL)
        cut_added = false
        infeasible_master = false
        xsize = size(A_all[1], 2)
        if cutmode == :AveragedCut
            averaged_optimality_cut_coef = zeros(Float64, xsize)
            averaged_optimality_cut_rhs = .0
        end
        #cuts_e = Vector{Float64}(0)
        #cuts_E = Matrix{Float64}(0, xsize)
        #cuts_d = Vector{Float64}(0)
        #cuts_D = Matrix{Float64}(0, xsize)
        # add cutting planes, one per scenario
        for i = 1:length(prob)
            coef = vec(output[i][4]' * A_all[i+1])
            rhs = vecdot(output[i][4], b_all[i+1])
            # add an infeasibility cut
            if output[i][1] == :Infeasible
                coef = vec(output[i][4]' * A_all[i+1])
                rhs = vecdot(output[i][4], b_all[i+1])
                # output[i][2] is a ray
                # so alpha * output[i][2] is also valid for any alpha >= 0.
                # Hence output[i][2] might have very large coefficients and alter
                # the numerial accuracy of the master's solver.
                # We scale it to avoid this issue
                scaling = abs(rhs)
                if scaling == 0
                    scaling = maximum(coef)
                end
                cut_added = true
                infeasible_master = true
                pushfeasibilitycut!(master_node, coef/scaling, sign(rhs))
                #cuts_d = [cuts_d; sign(rhs)]
                #cuts_d = [cuts_d; rhs/scaling]
                #cuts_D = [cuts_D; coef'/scaling]
                # add an optimality cut
            else
                if !infeasible_master
                    if cutmode == :AveragedCut
                        averaged_optimality_cut_coef += prob[i] * coef
                        averaged_optimality_cut_rhs += prob[i] * rhs
                    else
                        if !infeasible_master && θ[i] < (dot(coef, separator) - rhs)[1] - TOL
                            #@addConstraint(master_model, dot(coef, x) - θ[i] <= rhs)
                            cuts_e = [cuts_e; rhs]
                            cuts_E = [cuts_E'; coef]
                            cut_added = true
                        end
                    end
                end
            end
        end
        if cutmode == :AveragedCut && !infeasible_master && θ < (dot(averaged_optimality_cut_coef, separator) - averaged_optimality_cut_rhs)[1] - TOL
            #cuts_e = [cuts_e; averaged_optimality_cut_rhs]
            #cuts_E = [cuts_E; averaged_optimality_cut_coef']
            @show averaged_optimality_cut_rhs
            @show averaged_optimality_cut_coef'
            pushoptimalitycut!(master_node, averaged_optimality_cut_coef, averaged_optimality_cut_rhs)
            cut_added = true
        end
        return (cut_added, cuts_D, cuts_d, cuts_e, cuts_E)
    end

    function Benders_pmap(prob, c_all, A_all, B_all, b_all, K_all, C_all, v, master_solver, sub_solver, cutmode, TOL=1e-5)

        #println("Benders pmap started")
        num_master_var = length(c_all[1])
        num_scen = length(prob)

        num_bins = zeros(Int, num_scen)
        for i = 1:num_scen
            num_bins[i] = length(b_all[i+1])
        end

        #println("Load master problem")
        cut_added = true
        separator = zeros(num_master_var)
        objval = Inf
        status = :Infeasible
        niter = 0
        xsize = size(A_all[1], 2)
        cuts_e = Vector{Float64}(0)
        cuts_E = Matrix{Float64}(0, xsize)
        cuts_d = Vector{Float64}(0)
        cuts_D = Matrix{Float64}(0, xsize)
        while cut_added
            niter += 1
            #println("Iteration started")
            bigA = [A_all[1]; cuts_D]
            bigb = [b_all[1]; cuts_d]
            bigA = [bigA spzeros(size(bigA, 1), 1)]
            bigA = [bigA; cuts_E -ones(size(cuts_E, 1), 1)]
            bigb = [bigb; cuts_e]
            if size(bigA, 1) > size(A_all[1], 1)
                bigK = [K_all[1]; (:NonNeg, Vector{Int}((size(A_all[1], 1)+1):size(bigA, 1)))]
            else
                bigK = K_all[1]
            end
            bigC = [C_all[1]; (:NonNeg, [size(A_all[1], 2)+1])]
            bigc = [c_all[1]; 1]
            status, objval, primal, dual = loadAndSolveConicProblem(bigc, bigA, bigb, bigK, bigC, master_solver)
            separator = primal[1:size(A_all[1],2)]
            θ = primal[1+size(A_all[1],2)]
            if status == :Infeasible
                break
            end
            @show objval
            @show separator
            @show θ
            new_rhs = [zeros(num_bins[i]) for i in 1:num_scen]
            for i = 1:num_scen
                new_rhs[i] = b_all[i+1] - A_all[i+1] * separator
            end

            output = pmap((a1,a2,a3,a4,a5,a6)->loadAndSolveConicProblem(a1,a2,a3,a4,a5,a6), [c_all[i+1] for i = 1:num_scen], [B_all[i] for i = 1:num_scen], [new_rhs[i] for i = 1:num_scen], [K_all[i+1] for i = 1:num_scen], [C_all[i+1] for i = 1:num_scen], [sub_solver for i = 1:num_scen])

            (cut_added, new_cuts_D, new_cuts_d, new_cuts_e, new_cuts_E) = addCuttingPlanes(prob, A_all, b_all, output, θ, separator, cutmode, TOL)
            cuts_D = [cuts_D; new_cuts_D]
            cuts_d = [cuts_d; new_cuts_d]
            cuts_E = [cuts_E; new_cuts_E]
            cuts_e = [cuts_e; new_cuts_e]
        end
        @show niter
        return status, objval, separator
    end
