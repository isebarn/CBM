














@doc """
	convert_to_reversible!(model)

Converts the model by combining duplicate reactions that previously flowed in
 opposite directions

""" convert_to_reversible!

@doc """
	convert_to_reversible(model)

Converts a copy of the model by combining duplicate reactions that previously flowed in
 opposite directions

""" convert_to_reversible
















@doc """ 
    change_objective_coef(lp, objective, value)

Change the objective value for the entire lp model. ``objective`` may be 
a single index or an array of indices 

### Example
    change_objective_coef(lp, [2,5], [1, 0.5])
""" change_objective_coef

@doc """ 
    lb, ub = get_variable_bounds(lp, [index])

Return the column bounds for the lp model. Optionally, choose an index to return the bounds 
of a single column 

### Example 
    lb, ub = get_variable_bounds(lp, 13)

""" get_variable_bounds

@doc """ 
    b_lb, b_ub = get_constraint_bounds(lp)

Return the row bounds for the lp model
""" get_constraint_bounds

@doc """ 
    S = get_constraint_matrix(lp)

Return a sparse matrix filled with the coefficients of the reactions
""" get_constraint_matrix

@doc """ 
    direction = get_direction(lp)

Return either "maximize" or "minimize"
""" get_direction

@doc """ 
    obj_val = get_objective_value(lp)

Return the objective value of the lp model. Note that the model must have been solved
""" get_objective_value

@doc """ 
    red_cost = get_reduced_costs(lp)

Return the reduced cost vector of the lp model. Note that the model must have been solved
""" get_reduced_costs

@doc """ 
    sol_vector = get_solution(lp)

Return the solution vector of the lp model. Note that the model must have been solved
""" get_solution

@doc """ 
    sol_stat = get_solution_status(lp)

Return the solution status of the lp model. Returns one of the following 

* Optimal
* No feasible solution
* Infeasible
* Feasible
* Unbound
* Undefined

""" get_solution_status

@doc """ 
    stat_code = get_solution_status_code(lp)

Return the solution status code for a solved lp model. The solution status codes 
differ between solvers.
""" get_solution_status_code

@doc """ 
    slack = get_slack(lp)

Return the slack vector of the lp model. Note that the model must have been solved

""" get_slack

@doc """ 
    objective = get_objective(lp)

Return the objective vector of the lp model. This corresponds to model.c in the Model
""" get_objective

@doc """ 
    set_direction(lp, direction)

Set the direction of optimization. Direction may be "max" or "min"
""" set_direction

@doc """ 
    set_col_bounds(lp, index, lb, ub)
    set_col_bounds(lp, index, fixed_value)

Set the bounds for a column specified by **index** to either a lower- and upper- bound or to 
a **fixed value**

### Example 
    set_col_bounds(lp, 11, 0, 10)

To set the bounds for column (reaction) 11 to 0 and 10 
""" set_col_bounds

@doc """ 
    set_objective(lp, objective, value)

Change the objective value for a single column (reaction). This does not affect other reactions
""" set_objective

@doc """ 
    solve(lp)

Solve the lp model
""" solve
