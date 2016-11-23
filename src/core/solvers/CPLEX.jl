function Base.show(io::IO, m::CPLEX.Model)

    @printf "%0s \n" "CPLEX Model" 
    @printf "   %-15s : %10s \n" "sense " CPLEX.get_sense(m) == :Max ? "maximize" : "minimize"
    @printf "   %-15s = %10s \n" "variables" CPLEX.num_var(m)
    @printf "   %-15s = %10s \n" "constraints" CPLEX.num_constr(m)
    @printf "   %-15s = %10s \n" "coefficients" CPLEX.get_nnz(m)
end

# To change settings for CPLEX, find the parameter you wish 
# to change, uncomment its line, and change the X to the desired
# value.
# For quick (necessarily) reference for parameters, see
# http://www-eio.upc.edu/lceio/manuals/cplex-11/pdf/refparameterscplex.pdf
# for detailed specification of parameters 
#
# Note: You can also, instead of numbers, use fx 
#
# CPLEX.set_param!(env, "CPX_PARAM_SCRIND", CPX_ON)
#
# CPLEX.set_param!(env, "CPX_PARAM_SCRIND", 1)
#
function settings_cplex()
    env = CPLEX.Env()

    CPLEX.set_param!(env, "CPX_PARAM_LPMETHOD", 1) # 0=automatic
    CPLEX.set_param!(env, "CPX_PARAM_SCRIND", 1) # 1=display messages
    CPLEX.set_param!(env, "CPX_PARAM_SIMDISPLAY", 0) # 0=simplex msg off

    #
        #CPLEX.set_param!(env, CPX_PARAM_ADVIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_AGGCUTLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_AGGFILL, X)
        #CPLEX.set_param!(env, CPX_PARAM_AGGIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_ALL_MAX, X)
        #CPLEX.set_param!(env, CPX_PARAM_ALL_MIN, X)
        #CPLEX.set_param!(env, CPX_PARAM_APIENCODING, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARCOLNZ, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARCROSSALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_BAREPCOMP, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARGROWTH, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARMAXCOR, X)
        #CPLEX.set_param!(env, CPX_PARAM_BAROBJRNG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARORDER, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARQCPEPCOMP, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARSTARTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BBINTERVAL, X)
        #CPLEX.set_param!(env, CPX_PARAM_BNDSTRENIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_BQPCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_BRDIR, X)
        #CPLEX.set_param!(env, CPX_PARAM_BTTOL, X)
        #CPLEX.set_param!(env, CPX_PARAM_CALCQCPDUALS, X)
        #CPLEX.set_param!(env, CPX_PARAM_CLIQUES, X)
        #CPLEX.set_param!(env, CPX_PARAM_CLOCKTYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_CLONELOG, X)
        #CPLEX.set_param!(env, CPX_PARAM_COEREDIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_COLREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_CONFLICTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_CONFLICTDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_COVERS, X)
        #CPLEX.set_param!(env, CPX_PARAM_CPUMASK, X)
        #CPLEX.set_param!(env, CPX_PARAM_CRAIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTLO, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTPASS, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTSFACTOR, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTUP, X)
        #CPLEX.set_param!(env, CPX_PARAM_DATACHECK, X)
        #CPLEX.set_param!(env, CPX_PARAM_DEPIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_DETTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_DISJCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_DIVETYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_DPRIIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_EACHCUTLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPAGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPINT, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPLIN, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPMRK, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPOPT, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPPER, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPRELAX, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPRHS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FEASOPTMODE, X)
        #CPLEX.set_param!(env, CPX_PARAM_FILEENCODING, X)
        #CPLEX.set_param!(env, CPX_PARAM_FLOWCOVERS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FLOWPATHS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FPHEUR, X)
        #CPLEX.set_param!(env, CPX_PARAM_FRACCAND, X)
        #CPLEX.set_param!(env, CPX_PARAM_FRACCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FRACPASS, X)
        #CPLEX.set_param!(env, CPX_PARAM_GUBCOVERS, X)
        #CPLEX.set_param!(env, CPX_PARAM_HEURFREQ, X)
        #CPLEX.set_param!(env, CPX_PARAM_IMPLBD, X)
        #CPLEX.set_param!(env, CPX_PARAM_INTSOLFILEPREFIX, X)
        #CPLEX.set_param!(env, CPX_PARAM_INTSOLLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_ITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_LANDPCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_LBHEUR, X)
        #CPLEX.set_param!(env, CPX_PARAM_LOCALIMPLBD, X)
        #CPLEX.set_param!(env, CPX_PARAM_LPMETHOD, X)
        #CPLEX.set_param!(env, CPX_PARAM_MCFCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MEMORYEMPHASIS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPCBREDLP, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPEMPHASIS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPINTERVAL, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPKAPPASTATS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPORDIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPORDTYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPSEARCH, X)
        #CPLEX.set_param!(env, CPX_PARAM_MPSLONGNUM, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETEPOPT, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETEPRHS, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETFIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETPPRIIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_NODEFILEIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_NODELIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_NODESEL, X)
        #CPLEX.set_param!(env, CPX_PARAM_NUMERICALEMPHASIS, X)
        #CPLEX.set_param!(env, CPX_PARAM_NZREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_OBJDIF, X)
        #CPLEX.set_param!(env, CPX_PARAM_OBJLLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_OBJULIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_OPTIMALITYTARGET, X)
        #CPLEX.set_param!(env, CPX_PARAM_PARALLELMODE, X)
        #CPLEX.set_param!(env, CPX_PARAM_PERIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_PERLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERDETTIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTEREPAGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTEREPGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERINTSOL, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERNODE, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERTIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHTIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_POPULATELIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_PPRIIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_PREDUAL, X)
        #CPLEX.set_param!(env, CPX_PARAM_PREIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_PREPASS, X)
        #CPLEX.set_param!(env, CPX_PARAM_PROBETIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_QPMAKEPSDIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_QPMETHOD, X)
        #CPLEX.set_param!(env, CPX_PARAM_QPNZREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_QTOLININD, X)
        #CPLEX.set_param!(env, CPX_PARAM_RAMPUPDETTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_RAMPUPDURATION, X)
        #CPLEX.set_param!(env, CPX_PARAM_RAMPUPTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_REINV, X)
        #CPLEX.set_param!(env, CPX_PARAM_RELAXPREIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_RELOBJDIF, X)
        #CPLEX.set_param!(env, CPX_PARAM_REPAIRTRIES, X)
        #CPLEX.set_param!(env, CPX_PARAM_RINSHEUR, X)
        #CPLEX.set_param!(env, CPX_PARAM_ROWREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SCAIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_SCRIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIFTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIFTDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIFTITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIMDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_SINGLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SOLNPOOLAGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_SOLNPOOLGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_SOLUTIONTYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_STARTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_STRONGITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SUBALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_SYMMETRY, X)
        #CPLEX.set_param!(env, CPX_PARAM_THREADS, X)
        #CPLEX.set_param!(env, CPX_PARAM_TILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_TRELIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_TUNINGREPEAT, X)
        #CPLEX.set_param!(env, CPX_PARAM_TUNINGTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_VARSEL, X)
        #CPLEX.set_param!(env, CPX_PARAM_WORKDIR, X)
        #CPLEX.set_param!(env, CPX_PARAM_WORKMEM, X)
        #CPLEX.set_param!(env, CPX_PARAM_WRITELEVEL, X)
        #CPLEX.set_param!(env, CPX_PARAM_ZEROHALFCUTS, X)
    return env
end

function change_objective_coef(lp::CPLEX.Model, objective::Number, value::Number = 1.0)
	num_colums = CPLEX.num_var(lp)

	new_objective = zeros(num_colums)

	new_objective[objective] = value

	CPLEX.set_obj!(lp, new_objective)
end

function change_objective_coef{T <: Number}(lp::CPLEX.Model, objective::Array{T}, value::Array = [])
	if value == []
		value = ones(Float64, length(objective))
	end

	num_colums = CPLEX.num_var(lp)

	new_objective = zeros(num_colums)

	map((x,y) -> new_objective[x] = y, objective,value)

	CPLEX.set_obj!(lp, new_objective)
end

function get_variable_bounds(lp::CPLEX.Model)
	lb = CPLEX.get_varLB(lp)
	ub = CPLEX.get_varUB(lp)

	return lb, ub
end 

function get_variable_bounds(lp::CPLEX.Model, index::Number)
	lb = CPLEX.get_varLB(lp)[index]
	ub = CPLEX.get_varUB(lp)[index]	
	return lb,ub
end

function get_constraint_bounds(lp::CPLEX.Model)
	row_lb = CPLEX.get_constrLB(lp)
	row_ub = CPLEX.get_constrUB(lp)

	return row_lb, row_ub
end

function get_constraint_matrix(lp::CPLEX.Model)
	return CPLEX.get_constr_matrix(lp)
end 

function get_direction(lp::CPLEX.Model)
	direction = CPLEX.get_sense(lp)

	if CPLEX.get_sense(lp) == :Max
		return "maximize"
	else
		return "minimize"
	end
end

function get_objective_value(lp::CPLEX.Model)
	return CPLEX.get_objval(lp)
end

function get_reduced_costs(lp::CPLEX.Model)
	return CPLEX.get_reduced_costs(lp)
end

function get_solution(lp::CPLEX.Model)
	return CPLEX.get_solution(lp)
end

function get_solution_status(lp::CPLEX.Model)
	code = get_solution_status_code(lp)
	if code == 1 #1 CPX_STAT_OPTIMAL
		return "Optimal"
	elseif code == 4 #4 CPX_STAT_INForUNBD
		return "No feasible solution"
	elseif code == 3 #3 CPX_STAT_INFEASIBLE
		return "Infeasible" 
	elseif code == 5 #5 CPX_STAT_OPTIMAL_INFEAS
		return "Feasible" 
	elseif code == 2 #2 CPX_STAT_UNBOUNDED
		return "Unbound" 
	elseif code == 6 #6 CPX_STAT_NUM_BEST
		return "Undefined"
	else 
		return warn(string("Unknown status code: ", code))
	end 

	#7 CPX_STAT_FEASIBLE_RELAXED
	#8 CPX_STAT_OPTIMAL_RELAXED
	#10 CPX_STAT_ABORT_IT_LIM
	#11 CPX_STAT_ABORT_TIME_LIM
	#12 CPX_STAT_ABORT_OBJ_LIM
	#13 CPX_STAT_ABORT_USER
	#20 CPX_STAT_OPTIMAL_FACE_UNBOUNDED
end

function get_solution_status_code(lp::CPLEX.Model)
	return getfield(CPLEX, CPLEX.get_status(lp))
end 

function get_slack(lp::CPLEX.Model)
	return CPLEX.get_constr_duals(lp)
end 

function get_objective(lp::CPLEX.Model)
	return CPLEX.get_obj(lp)
end

function set_direction(lp::CPLEX.Model, direction::AbstractString)
	if lowercase(direction) == "max"
		CPLEX.set_sense!(lp, :Max)
	else
		CPLEX.set_sense!(lp, :Min)
	end
end

function set_col_bounds(lp::CPLEX.Model, index::Number, lb::Number, ub::Number)
	new_lb = CPLEX.get_varLB(lp)
	new_lb[index] = lb
	CPLEX.set_varLB!(lp, new_lb)

	new_ub = CPLEX.get_varUB(lp)	
	new_ub[index] = ub
	CPLEX.set_varUB!(lp, new_ub)
end

function set_col_bounds(lp::CPLEX.Model, index::Number, fixed_value::Number)
	new_lb = CPLEX.get_varLB(lp)
	new_lb[index] = fixed_value
	CPLEX.set_varLB!(lp, new_lb)

	new_ub = CPLEX.get_varUB(lp)	
	new_ub[index] = fixed_value
	CPLEX.set_varUB!(lp, new_ub)
end

function set_objective(lp::CPLEX.Model, objective::Number, value::Number = 1.0)
	new_objective = CPLEX.get_obj(lp)
	new_objective[objective] = value 

	CPLEX.set_obj!(lp, new_objective)
end

function setup_cplex(model::Model, objective::AbstractString = "max")
  # Flux Balance Analysis usingh CPLEX.jl
	nrxns = length(model.c)
	nmets = length(model.b)

	# Construct the CPLEX model
	env = settings_cplex()
	lp = CPLEX.Model(env, "CPLEX")

	if objective == "max"
		set_direction(lp, "max")
	else
		set_direction(lp, "min")
	end 

    CPLEX.add_vars!(lp, model.c, model.lb, model.ub)

    ####
    if length(model.b) == length(model.csense["="])
	   CPLEX.add_constrs!(lp, model.S, fill('=', nmets), model.b)
    else 
        b_lower = deepcopy(model.b)
        b_upper = deepcopy(model.b)
        b_lower[get(model.csense,"<=", [])] = -Inf
        b_upper[get(model.csense,">=", [])] = Inf
        CPLEX.add_rangeconstrs!(lp, model.S, b_lower, b_upper)
    end 
    ####


	return lp
end

function setup_cplex{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
  # Flux Balance Analysis usingh CPLEX.jl

	nmets = length(b)

	# Construct the CPLEX model
	env = CPLEX.Env()
	lp = CPLEX.Model(env, "CPLEX")
	CPLEX.set_param!(env, "CPX_PARAM_LPMETHOD", CPLEX.CPX_ALG_AUTOMATIC) # should?


	if objective == "max"
		set_direction(lp, "max")
	else
		set_direction(lp, "min")
	end 
	
	CPLEX.add_vars!(lp, c, lb, ub)

	CPLEX.add_constrs!(lp, S, fill('=', nmets), b)


	return lp
end

function setup_cplex{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b_lower::Array{T}, b_upper::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
    nmets = length(b_lower)

    # Construct the CPLEX model
    env = CPLEX.Env()
    lp = CPLEX.Model(env, "CPLEX")
    CPLEX.set_param!(env, "CPX_PARAM_LPMETHOD", CPLEX.CPX_ALG_AUTOMATIC) # should?


    if objective == "max"
        set_direction(lp, "max")
    else
        set_direction(lp, "min")
    end 

    CPLEX.add_vars!(lp, c, lb, ub)

    CPLEX.add_rangeconstrs!(lp, S, b_lower, b_upper)
	
    return lp
end

function setup_cplex{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, csense::Array{String}, S::SparseMatrixCSC, objective::AbstractString = "max")
    b_lower = deepcopy(b)
    b_upper = deepcopy(b)
    b_lower[find(csense .== "<=")] = -Inf
    b_upper[find(csense .== ">=")] = Inf

    setup_cplex(c, lb, ub, b_lower, b_upper, S, objective)
end 

function setup_cplex{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, csense::Dict, S::SparseMatrixCSC, objective::AbstractString = "max")
    b_lower = deepcopy(b)
    b_upper = deepcopy(b)
    b_lower[get(csense,"<=", [])] = -Inf
    b_upper[get(csense,">=", [])] = Inf

    setup_cplex(c, lb, ub, b_lower, b_upper, S, objective)

end 

function setup_cplex{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, beq::Array{T}, Seq::SparseMatrixCSC, bineq::Array{T}, Sineq::SparseMatrixCSC, objective::AbstractString = "max")
  # Flux Balance Analysis usingh CPLEX.jl


    # Construct the CPLEX model
    env = CPLEX.Env()
    lp = CPLEX.Model(env, "CPLEX")
    CPLEX.set_param!(env, "CPX_PARAM_LPMETHOD", CPLEX.CPX_ALG_AUTOMATIC) # should?


    if objective == "max"
        set_direction(lp, "max")
    else
        set_direction(lp, "min")
    end 
    
    CPLEX.add_vars!(lp, c, lb, ub)

    CPLEX.add_rangeconstrs!(lp, Sineq, fill(-Inf, length(bineq)), bineq)
    CPLEX.add_constrs!(lp, Seq, fill('=', length(beq)), beq)

    return lp
end

function solve(lp::CPLEX.Model)
	CPLEX.optimize!(lp)
end

function solver_status(lp::CPLEX.Model)
    string(CPLEX.status_symbols[get_solution_status_code(lp)])
end 

