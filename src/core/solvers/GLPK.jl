function Base.show(io::IO, m::GLPK.Prob)

    @printf "%0s \n" "GLPK Model" 
    @printf "   %-15s : %10s \n" "sense " GLPK.get_obj_dir(m) == 2 ? "maximize" : "minimize"
    @printf "   %-15s = %10s \n" "variables" GLPK.get_num_cols(m)
    @printf "   %-15s = %10s \n" "constraints" GLPK.get_num_rows(m)
    @printf "   %-15s = %10s \n" "coefficients" GLPK.get_num_nz(m)
end


# status codes 
# http://artax.karlin.mff.cuni.cz/r-help/library/glpkAPI/html/glpkConstants.html
# To change settings for GLPK, find the parameter you wish 
# to change, uncomment its line, and change the X to the desired
# value.
# Note: You can also, instead of numbers, use fx 
#
# setfield!(P, Symbol("msg_lev"), GLPK.MSG_OFF)
# as equivalent to
# setfield!(P, Symbol("msg_lev"),convert(Int32, 1))
#
function settings_glpk()
    P = SimplexParam()
    setfield!(P, Symbol("msg_lev"),convert(Int32, 1))
    #setfield!(P, Symbol("meth"),convert(Int32, X))
    #setfield!(P, Symbol("pricing"),convert(Int32, X))
    #setfield!(P, Symbol("r_test"),convert(Int32, X))
    #setfield!(P, Symbol("tol_bnd"),convert(Float64, X))
    #setfield!(P, Symbol("tol_dj"),convert(Float64, X))
    #setfield!(P, Symbol("tol_piv"),convert(Float64, X))
    #setfield!(P, Symbol("obj_ll"),convert(Float64, X))
    #setfield!(P, Symbol("obj_ul"),convert(Float64, X))
    #setfield!(P, Symbol("it_lim"),convert(Int32, X))
    #setfield!(P, Symbol("tm_lim"),convert(Int32, X))
    #setfield!(P, Symbol("out_frq"),convert(Int32, X))
    #setfield!(P, Symbol("out_dly"),convert(Int32, X))
    #setfield!(P, Symbol("presolve"),convert(Int32, 1))
    return P
end

function change_objective_coef(lp::GLPK.Prob, objective::Number, value::Number = 1.0)
	num_colums = GLPK.get_num_cols(lp)

	new_objective = zeros(num_colums)

	new_objective[objective] = value

	map((x -> GLPK.set_obj_coef(lp, x, new_objective[x])), collect(1:num_colums))
end

function change_objective_coef{T <: Number}(lp::GLPK.Prob, objective::Array{T}, value::Array=[])
	if value == []
		value = ones(Float64, length(objective))
	end 
	num_colums = GLPK.get_num_cols(lp)

	new_objective = zeros(num_colums)
	new_objective[objective] = value 
	#map((x,y) -> new_objective[x] = y, objective,value)

	map((x -> GLPK.set_obj_coef(lp, x, new_objective[x])), collect(1:num_colums));
end

function get_variable_bounds(lp::GLPK.Prob)
	num_rxns = GLPK.get_num_cols(lp)
	lb = map(x -> GLPK.get_col_lb(lp,x), collect(1:num_rxns))
	lb = vcat(lb...)
	ub = map(x -> GLPK.get_col_ub(lp,x), collect(1:num_rxns))
	ub = vcat(ub...)

	return lb, ub
end

function get_variable_bounds(lp::GLPK.Prob, index::Number)
	lb = GLPK.get_col_lb(lp, index)
	ub = GLPK.get_col_ub(lp, index)
	return lb,ub
end

function get_constraint_bounds(lp::GLPK.Prob)
	num_mets = GLPK.get_num_rows(lp)
	lb = map(x -> GLPK.get_row_lb(lp,x), collect(1:num_mets))
	lb = vcat(lb...)
	ub = map(x -> GLPK.get_row_ub(lp,x), collect(1:num_mets))
	ub = vcat(ub...)

	return lb, ub
end

function get_constraint_matrix(lp::GLPK.Prob)
	num_rxns = GLPK.get_num_cols(lp)
	num_mets = GLPK.get_num_rows(lp)
	num_elements = GLPK.get_num_nz(lp)
	ar = Array(Float64, num_elements)

	ia = map(x -> GLPK.get_mat_col(lp,x)[1], collect(1:num_rxns));
	ia = vcat(ia...)

	ar = map(x -> GLPK.get_mat_col(lp,x)[2], collect(1:num_rxns));
	ar = vcat(ar...)

	ja = map(x -> GLPK.get_mat_row(lp,x)[1], collect(1:num_mets));
	ja = vcat(ja...)
	sort!(ja)

	S = sparse(ia,ja,ar)
end

function get_direction(lp::GLPK.Prob)
	if GLPK.get_obj_dir(lp) == 2
		return "maximize"
	else 
		return "minimize"
	end 
end

function get_objective_value(lp::GLPK.Prob)
	return GLPK.get_obj_val(lp)
end

function get_reduced_costs(lp::GLPK.Prob)
	num_rxns = GLPK.get_num_cols(lp)
	return map((x -> get_col_dual(lp, x)), collect(1:num_rxns))
end

function get_solution(lp::GLPK.Prob)  
	num_colums = GLPK.get_num_cols(lp)
	x = map((x -> get_col_prim(lp, x)), collect(1:num_colums))
end

function get_solution_status(lp::GLPK.Prob)
	code = get_solution_status_code(lp)
	if code == 5
		return "Optimal"
	elseif code == 4
		return "No feasible solution"
	elseif code == 3
		return "Infeasible"
	elseif code == 2
		return "Feasible"
	elseif code == 6
		return "Unbound"
	elseif code == 1
		return "Undefined"
	else 
		return warn(string("Unknown status code: ", code))
	end 
end

function get_solution_status_code(lp::GLPK.Prob)
	return GLPK.get_status(lp)
end 

function get_slack(lp::GLPK.Prob)
	num_mets = GLPK.get_num_rows(lp)
	return map((x -> get_row_dual(lp, x)), collect(1:num_mets))
end 

function get_objective(lp::GLPK.Prob)
	num_rxns = GLPK.get_num_cols(lp)
	return map(x -> GLPK.get_obj_coef(lp,x), collect(1:num_rxns))
end

function set_direction(lp::GLPK.Prob, direction::AbstractString)
	if lowercase(direction) == "max"
		GLPK.set_obj_dir(lp, GLPK.MAX);
	else
		GLPK.set_obj_dir(lp, GLPK.MIN);
	end
end

function set_col_bounds(lp::GLPK.Prob, index::Number, lb::Number, ub::Number)
	if lb == ub
		GLPK.set_col_bnds(lp, index, GLPK.FX, lb, ub)
	else 
		GLPK.set_col_bnds(lp, index, GLPK.DB, lb, ub)
	end
end

function set_col_bounds(lp::GLPK.Prob, index::Number, fixed_value::Number)
	GLPK.set_col_bnds(lp, index, GLPK.FX, fixed_value, fixed_value)
end

function set_objective(lp::GLPK.Prob, objective::Number, value::Number = 1.0)
	GLPK.set_obj_coef(lp, objective, value)
end

function setup_glpk(model::Model, objective::AbstractString = "max")
	lp = GLPK.Prob()
	GLPK.set_obj_dir(lp, GLPK.MAX)

	if contains(lowercase(objective),"min")
		GLPK.set_obj_dir(lp, GLPK.MIN)
	end

	GLPK.add_cols(lp, length(model.c))
	for i = 1:length(model.c)
		GLPK.set_col_bnds(lp, i, GLPK.DB, model.lb[i],model.ub[i])
		if model.lb[i] == model.ub[i]
			GLPK.set_col_bnds(lp, i, GLPK.FX, model.lb[i],model.ub[i])
		end
		GLPK.set_obj_coef(lp,i,model.c[i])
	end

	GLPK.add_rows(lp, length(model.b))

	for i in model.csense["<="]
		GLPK.set_row_bnds(lp,i,GLPK.UB,model.b[i],model.b[i])
	end 
	for i in model.csense[">="]
		GLPK.set_row_bnds(lp,i,GLPK.LO,model.b[i],model.b[i])
	end 
	for i in model.csense["="]
		GLPK.set_row_bnds(lp,i,GLPK.FX,model.b[i],model.b[i])
	end	

	GLPK.load_matrix(lp, model.S)

	return lp
end

function setup_glpk{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T},  S::SparseMatrixCSC, objective::AbstractString = "max")
	setup_glpk(c, lb, ub, b, b, S, objective)
end

function setup_glpk{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b_lower::Array{T}, b_upper::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
	lp = GLPK.Prob()
	
	GLPK.set_prob_name(lp, "sample")
	GLPK.term_out(GLPK.ON)
	
	GLPK.set_obj_dir(lp, GLPK.MAX)

	if contains(lowercase(objective),"min")
		GLPK.set_obj_dir(lp, GLPK.MIN)
	end

	GLPK.add_cols(lp, length(c))
	for i = 1:length(c)
		GLPK.set_col_bnds(lp, i, GLPK.DB, lb[i],ub[i])
		if lb[i] == ub[i]
			GLPK.set_col_bnds(lp, i, GLPK.FX, lb[i],ub[i])
		end
		GLPK.set_obj_coef(lp,i,c[i])
	end

	GLPK.add_rows(lp, length(b_lower))
	for i = 1:length(b_lower)
		if b_lower[i] == b_upper[i]
			GLPK.set_row_bnds(lp,i,GLPK.FX,b_lower[i],b_upper[i])
		else 
			GLPK.set_row_bnds(lp,i,GLPK.DB,b_lower[i],b_upper[i])
		end 
	end	

	GLPK.load_matrix(lp, S)

	return lp
end 

function setup_glpk{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, csense::Array{String}, S::SparseMatrixCSC, objective::AbstractString = "max")
    b_lower = deepcopy(b)
    b_upper = deepcopy(b)
    b_lower[find(csense .== "<=")] = -Inf
    b_upper[find(csense .== ">=")] = Inf
    return setup_glpk(c, lb, ub, b_lower, b_upper, S, objective)
end 

function setup_glpk{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, csense::Dict, S::SparseMatrixCSC, objective::AbstractString = "max")
    b_lower = deepcopy(b)
    b_upper = deepcopy(b)
    b_lower[get(csense,"<=", [])] = -Inf
    b_upper[get(csense,">=", [])] = Inf
    return setup_glpk(c, lb, ub, b_lower, b_upper, S, objective)
end 

function setup_glpk{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, beq::Array{T}, Seq::SparseMatrixCSC, bineq::Array{T}, Sineq::SparseMatrixCSC, objective::AbstractString = "max")
	lp = GLPK.Prob()
	
	GLPK.set_prob_name(lp, "sample")
	GLPK.term_out(GLPK.ON)
	
	GLPK.set_obj_dir(lp, GLPK.MAX)

	if contains(lowercase(objective),"min")
		GLPK.set_obj_dir(lp, GLPK.MIN)
	end

	GLPK.add_cols(lp, length(c))
	for i = 1:length(c)
		GLPK.set_col_bnds(lp, i, GLPK.DB, lb[i],ub[i])
		if lb[i] == ub[i]
			GLPK.set_col_bnds(lp, i, GLPK.FX, lb[i],ub[i])
		end
		GLPK.set_obj_coef(lp,i,c[i])
	end

	GLPK.add_rows(lp, length(beq)+length(bineq))
	for i = 1:length(beq)
		GLPK.set_row_bnds(lp,i,GLPK.FX,beq[i],beq[i])
	end	

	#GLPK.add_rows(lp, length(bineq))
	for i = 1:length(bineq)
		GLPK.set_row_bnds(lp,i+length(beq),GLPK.LO,-bineq[i],0.0)
	end	
	GLPK.load_matrix(lp, vcat(Seq,Sineq))


	return lp
end

function solve(lp::GLPK.Prob)
	return GLPK.simplex(lp, settings_glpk())
end

function solver_status(lp::GLPK.Prob)

	status_code = get_solution_status_code(lp)

	if status_code == 1
		return "GLP_UNDEF" 
	elseif status_code == 2
		return "GLP_FEAS" 
	elseif status_code == 3
		return "GLP_INFEAS" 
	elseif status_code == 4
		return "GLP_NOFEAS" 
	elseif status_code == 5
		return "GLP_OPT" 
	elseif status_code == 6
		return "GLP_UNBND" 
	end
end











