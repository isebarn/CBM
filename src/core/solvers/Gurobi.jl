function Base.show(io::IO, m::Gurobi.Model)

    @printf "%0s \n" "Gurobi Model" 
    @printf "   %-15s : %10s \n" "sense " Gurobi.model_sense(m) == :maximize ? "maximize" : "minimize"
    @printf "   %-15s = %10s \n" "variables" Gurobi.num_vars(m)
    @printf "   %-15s = %10s \n" "constraints" Gurobi.num_constrs(m)
    @printf "   %-15s = %10s \n" "coefficients" Gurobi.num_cnzs(m)
end

function change_objective_coef(lp::Gurobi.Model, objective::Number, value::Number = 1.0)
	obj = zeros(Float64, Gurobi.num_vars(lp))
	obj[objective] = value 
	Gurobi.set_objcoeffs!(lp, obj)
	update_model!(lp)
end

function change_objective_coef{T <: Number}(lp::Gurobi.Model, objective::Array{T}, value::Array=[])
	if value == []
		value = ones(Float64, length(objective))
	end 
	obj = zeros(Float64, Gurobi.num_vars(lp))
	obj[objective] = value 
	Gurobi.set_objcoeffs!(lp, obj)
	update_model!(lp)
end 

function get_variable_bounds(lp::Gurobi.Model)
	lb = Gurobi.get_dblattrarray(lp, "lb", 1, Gurobi.num_vars(lp))
	ub = Gurobi.get_dblattrarray(lp, "ub", 1, Gurobi.num_vars(lp))
	return lb, ub 
end 

function get_variable_bounds(lp::Gurobi.Model, index::Number)
	lb = Gurobi.get_dblattrarray(lp, "lb", index, index)[1]
	ub = Gurobi.get_dblattrarray(lp, "ub", index, index)[1]
	return lb, ub
end

function get_constraint_bounds(lp::Gurobi.Model)
	b = Gurobi.get_dblattrarray(lp, "RHS", 1, Gurobi.num_constrs(lp))
	return b,b 
end 

function get_constraint_matrix(lp::Gurobi.Model)
	return Gurobi.get_constrmatrix(lp)
end 

function get_direction(lp::Gurobi.Model)
	Gurobi.model_sense(lp) == :maximize ? "maximize" : "minimize"
end

function get_objective_value(lp::Gurobi.Model)
	return Gurobi.get_objval(lp)
end

function get_reduced_costs(lp::Gurobi.Model)
	Gurobi.get_dblattrarray(lp, "RC", 1, Gurobi.num_vars(lp))
end 

function get_solution(lp::Gurobi.Model)
	return Gurobi.get_solution(lp)
end

function get_solution_status(lp::Gurobi.Model)
	code = get_solution_status_code(lp)
	if code == 2 
		return "Optimal"
	elseif code == 3
		return "Infeasible"
	elseif code == 4
		return "Infeasible or Unbounded"
	elseif code == 5
		return "Unbound"
	elseif code == 13 #SUBOPTIMAL
		return "Feasible"
	else
		return warn(string("Unknown status code: ", code))
	end 
end

function get_solution_status_code(lp::Gurobi.Model)
	return Gurobi.get_status_code(lp)
end

function get_slack(lp::Gurobi.Model)
	Gurobi.get_dblattrarray(lp, "Pi", 1,Gurobi.num_constrs(lp))
end

function get_objective(lp::Gurobi.Model)
	Gurobi.get_dblattrarray(lp, "Obj", 1, Gurobi.num_vars(lp))
end 

function set_direction(lp::Gurobi.Model, direction::AbstractString)
	lowercase(direction) == "max" ? Gurobi.set_sense!(lp, :maximize) : Gurobi.set_sense!(lp, :minimize)
	update_model!(lp)
end

function set_col_bounds(lp::Gurobi.Model, index::Number, lb::Number, ub::Number)
	lower, upper = get_variable_bounds(lp)
	lower[index] = lb
	upper[index] = ub
	Gurobi.set_dblattrarray!(lp, "LB", 1, num_vars(lp), lower)
	Gurobi.set_dblattrarray!(lp, "UB", 1, num_vars(lp), upper)
	update_model!(lp)
end 

function set_col_bounds(lp::Gurobi.Model, index::Number, fixed_value::Number)
	set_col_bounds(lp, index, fixed_value, fixed_value)
end 

function set_objective(lp::Gurobi.Model, objective::Number, value::Number = 1.0)
	obj = Gurobi.get_dblattrarray(lp, "Obj", 1, num_vars(lp))
	obj[objective] = value
	Gurobi.set_dblattrarray!(lp, "Obj", 1, num_vars(lp), obj)
	update_model!(lp)
end

function setup_gurobi(model::Model, objective::AbstractString = "max")
	env = Gurobi.Env()  
	setparam!(env,"OutputFlag",0)
	setparam!(env,"Method",0)

	lp = Gurobi.Model(env, "god", :maximize)

	add_cvars!(lp, model.c, model.lb, model.ub) 
	update_model!(lp)

	add_constrs!(lp, model.S, fill('=', length(model.b)), model.b)

	update_model!(lp)

	return lp	
end

#function setup_gurobi{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
#function setup_gurobi{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, beq::Array{T}, Seq::SparseMatrixCSC, bineq::Array{T}, Sineq::SparseMatrixCSC, objective::AbstractString = "max")

function solve(lp::Gurobi.Model)
	return Gurobi.optimize(lp)
end

