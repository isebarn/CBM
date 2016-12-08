@everywhere module PCBM  # defined on all workers
    type Model
        rxns
        mets
        genes
        S
        lb
        ub
        c
        b
        csense
        rxn_gene_mat
        rxn_name
        rxn_rules
        rxn_subsystem
        rxn_extra
        met_formula
        met_name
        met_extra
        gene_name
        gene_extra
        description
    end

    solvers = map(x -> in(x, readdir(Pkg.dir())) ? x : "" , ["GLPK", "CPLEX"])
    
    function set_solver(solver)
        Solver = ""
        if solver == "glpk"
            Solver = "GLPK"
        elseif solver == "cplex"
            Solver = "CPLEX"
        else 
            warn("choose one of the following solvers:")
            warn("glpk")
            warn("cplex")
            return
        end 

        global settings_default_lp = solver
    end

    if in("GLPK", solvers)
        using GLPK
        include(Pkg.dir() * "/CBM/src/core/solvers/GLPK.jl")
        global settings_default_lp = "glpk"
    end 

    if in("CPLEX", solvers)
        using CPLEX
        include(Pkg.dir() * "/CBM/src/core/solvers/CPLEX.jl")
    end

    if !isdefined(:settings_default_lp)
        if in("GLPK", solvers)
            global settings_default_lp = "glpk"
        elseif in("CPLEX", solvers)
            global settings_default_lp = "cplex"
        else 
            warn("Neither GLPK or CPLEX is installed")
        end 
    end 

    function global_setup_lp{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b_lower::Array{T}, b_upper::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
        global lp = PCBM.setup_lp(c,lb,ub,b_lower, b_upper, S, objective)
    end

    function setup_lp{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b_lower::Array{T}, b_upper::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
        if settings_default_lp == "glpk"
            return PCBM.setup_glpk(c,lb,ub,b_lower, b_upper, S, objective)
        elseif settings_default_lp == "cplex"
            return PCBM.setup_cplex(c,lb,ub,b_lower, b_upper, S, objective)
        else 
            warn("No solver")
        end 
    end

    ###############################################################################
    # stuff for FVA in Parallel 
    ###############################################################################
    # THINK setup_workers is perhaps uneccesary check when you have time
    function setup_workers(model)
        @sync for pid in workers()
            @async remotecall_fetch(PCBM.global_setup_lp, pid, 
                zeros(Float64, length(model.rxns)),
                    model.lb,
                    model.ub,
                    model.b,
                    model.b,
                    model.S)
        end 
    end 

    function fast_fva_objective(x)
        res = Array(Float64, length(x))

        for (i,v) in enumerate(x) 
            PCBM.set_objective(lp, v, 1.0)

            PCBM.solve(lp)
            res[i] = PCBM.get_objective_value(lp)

            PCBM.set_objective(lp, v, 0.0)
        end 

        return res
    end 

    function fast_fva_solution(x)
        res = Array(Float64, length(get_variable_bounds(lp)[1]), length(x))

        for (i,v) in enumerate(x) 
            PCBM.set_objective(lp, v, 1.0)

            PCBM.solve(lp)
            res[:,i] = PCBM.get_solution(lp)

            PCBM.set_objective(lp, v, 0.0)
        end 

        return res
    end 

    function fast_fva(model, optPercentage::Number = 1, flux_matrix::Bool = false)
        nrxns = length(model.rxns)
        x = collect(1:nrxns)


        lpmin = setup_lp(model.c, model.lb, model.ub, model.b,model.b, model.S)
        solve(lpmin)
        wt_min = get_objective_value(lpmin) 

        lb = deepcopy(model.lb)
        lb[find(model.c)] *= wt_min

        @sync for pid in workers()
            @async remotecall_fetch(PCBM.global_setup_lp, pid, 
                zeros(Float64, length(model.rxns)),
                    lb,
                    model.ub,
                    model.b,
                    model.b,
                    model.S)
        end         

        stepsize = div(nrxns, 3*length(workers()))
        checklist = map(x -> collect(x:x+stepsize-1), 1:stepsize:nrxns)
        checklist[end] = checklist[end][checklist[end] .<= nrxns]


        maxF = flux_matrix ?  pmap(PCBM.fast_fva_solution, checklist) : pmap(PCBM.fast_fva_objective, checklist)

        map(x -> remotecall_fetch(() -> PCBM.set_direction(lp, "min"), x), workers())

        minF = flux_matrix ?  pmap(PCBM.fast_fva_solution, checklist) : pmap(PCBM.fast_fva_objective, checklist)

        if flux_matrix
            maxV = hcat(maxF...)
            minV = hcat(minF...)
            maxF = diag(maxV)
            minF = diag(minV)
            return minF, maxF, minV, maxV
        end 

        return vcat(minF...), vcat(maxF...)
    end 
    ###############################################################################
    ###############################################################################

    ###############################################################################
    # stuff for gene deletion in Parallel 
    ###############################################################################
    function setup_workers_for_gene_deletion(model)
        @sync for pid in workers()
            @async remotecall_fetch(PCBM.global_setup_lp, pid, 
                model.c,
                model.lb,
                model.ub,
                model.b,
                model.b,
                model.S)
        end 
    end 

    function pgd_analysis(rxns)
        #progress_meter = ProgressMeter.Progress(length(rxns), 0.5, "\n\u1b[1FFVA:", 50)
        combos = sort(rxns, lt=lexless)
        num_combos = length(combos)

        solve(lp)
        wt = get_objective_value(lp)

        lb, ub = get_variable_bounds(lp)

        flow = Array(Float64, num_combos)

        ####
        for d in combos[1] 
            set_col_bounds(lp, d, 0.0, 0.0)
        end 

        solve(lp)
        flow[1] = get_objective_value(lp)

        disable = []
        for i in 2:num_combos
            enable = setdiff(combos[i-1], combos[i])
            disable = setdiff(combos[i], combos[i-1])

            for e in enable 
                set_col_bounds(lp, e, lb[e], ub[e])
            end
            
            for d in disable 
                set_col_bounds(lp, d, 0.0, 0.0)
            end 
            solve(lp)
            flow[i] = get_objective_value(lp)
            #next!(progress_meter)
        end

        for e in disable 
            set_col_bounds(lp, e, lb[e], ub[e])
        end     

        flow /= wt
        flow = round(flow, 16)

        result = Dict(zip(combos, flow))

        map(x -> set_col_bounds(lp, x, lb[x], ub[x]), 1:length(lb))

        return result
    end 

    function pgd(model, rxns)
        listlen = length(rxns)
        setup_workers_for_gene_deletion(model)

        stepsize = div(listlen, 3*length(workers()))
        checklist = map(x -> collect(x:x+stepsize-1), 1:stepsize:listlen)
        checklist[end] = checklist[end][checklist[end] .<= listlen]

        rxnlist = map(x -> rxns[x], checklist)

        res = pmap(PCBM.pgd_analysis, rxnlist)
    end 
    ###############################################################################
    ###############################################################################

end
