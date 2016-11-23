function Base.show(io::IO, m::Clp.ClpCInterface.ClpModel)

    @printf "%0s \n" "Clp Model" 
    @printf "   %-15s : %10s \n" "sense " Clp.ClpCInterface.optimization_direction(m) == -1 ? "maximize" : "minimize"
    @printf "   %-15s = %10s \n" "variables" Clp.ClpCInterface.get_num_cols(m)
    @printf "   %-15s = %10s \n" "constraints" Clp.ClpCInterface.get_num_rows(m)
    @printf "   %-15s = %10s \n" "coefficients" Clp.ClpCInterface.get_num_elements(m)
end

function settings_clp()
    lp = Clp.ClpCInterface.ClpModel()

    Clp.ClpCInterface.set_log_level(lp,0)
    return lp
end

function change_objective_coef(lp::ClpModel, objective::Number, value::Number = 1.0)
    obj = zeros(Float64, number_cols(lp))
    obj[objective] = value
    chg_obj_coefficients(lp, obj)
end

function change_objective_coef{T <: Number}(lp::ClpModel, objective::Array{T}, value::Array = [])
    if value == []
        value = ones(Float64, length(objective))
    end

    obj = zeros(Float64, number_cols(lp))
    obj[objective] = value
    chg_obj_coefficients(lp, obj)
end

function get_variable_bounds(lp::ClpModel)
    return get_col_lower(lp), get_col_upper(lp)
end

function get_variable_bounds(lp::ClpModel, index::Number)
    lb,ub = get_variable_bounds(lp)
    return lb[index], ub[index]
end

function get_constraint_bounds(lp::ClpModel)
    return get_row_lower(lp), get_row_upper(lp)
end

function get_constraint_matrix(lp::ClpModel)
    ClpCInterface.get_constraint_matrix(lp)
end

function get_direction(lp::ClpModel)
   optimization_direction(lp) == -1 ? "maximize" : "minimize"
end

function get_objective_value(lp::ClpModel)
    objective_value(lp)
end

function get_reduced_costs(lp::ClpModel)
    get_reduced_cost(lp)
end

function get_solution(lp::ClpModel)
    get_col_solution(lp)
end

function get_solution_status(lp::ClpModel)
    code = get_solution_status_code(lp)
    if code == 0
        return "Optimal"
    else 
        return "Not optimal"
    end
end 

function get_solution_status_code(lp::ClpModel)
    status(lp)
end

function get_slack(lp::ClpModel)
    dual_column_solution(lp)
end

function get_objective(lp::ClpModel)
    get_obj_coefficients(lp)
end

function set_direction(lp::ClpModel, direction::AbstractString)
    set_optimization_direction(lp, lowercase(direction) == "max" ? -1.0 : 1.0)
end

function set_col_bounds(lp::ClpModel, index::Number, lb::Number, ub::Number)
    lower, upper = get_variable_bounds(lp)
    lower[index] = lb 
    upper[index] = ub
    chg_column_lower(lp, lower)
    chg_column_upper(lp, upper)
end

function set_col_bounds(lp::ClpModel, index::Number, fixed_value::Number)
    set_col_bounds(lp, index, fixed_value, fixed_value)
end

function set_objective(lp::ClpModel, objective::Number, value::Number = 1.0)
    obj = Clp.ClpCInterface.objective(lp)
    obj[objective] = value
    chg_obj_coefficients(lp, obj)
end

function setup_clp(model::Model, objective::AbstractString = "max")
    lp = settings_clp()

    load_problem(lp,model.S, model.lb, model.ub, model.c, model.b,model.b)
    set_direction(lp, objective)
    
    return lp
end 

#function setup_clp{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
#end

#function setup_clp{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, beq::Array{T}, Seq::SparseMatrixCSC, bineq::Array{T}, Sineq::SparseMatrixCSC, objective::AbstractString = "max")
#end

function solve(lp::ClpModel)
    primal(lp,0)
end

function solver_status(lp::ClpModel)
    get_solution_status(lp)
end 


