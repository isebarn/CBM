

@doc """
    fva(model,[optPercentage = 100])

Determine the flux variability of the model

### Optional arguements:
    optPercentage: 0 - 1 (fraction)
    
Used to fix the biomass production rate anywhere between 0 - 1


#### Examples
Run flux variability with biomass lower bound at 90% of maximum

    min, max, minX, maxX = fva(model, 0.9)

Return a min flux vector, max flux vector, minimum primal flux vector and maximum 
primal flux vector
""" fva









@doc """
    add_metabolite!(model, metabolite)

Add a metabolite to the model
""" add_metabolite!

@doc """
    add_reaction(model,  reaction_name, reaction, [lb == -1000, ub = 1000])

Add a reaction to a copy of the model

### Optional arguements:
    lb \& ub
* add lower/upper bounds, lower bound must be less/equal to upper bound

Return the new model with the added reaction.

*Note: Original model remains unchanged...see add_reaction! to permanently add a reaction*
""" add_reaction


@doc """
    add_reaction!(model, reaction_name, reaction, [lb == -1000, ub = 1000])

Permanently add a reaction to the model

### Optional arguements:
    lb \& ub
* add lower/upper bounds, lower bound must be less/equal to upper bound
""" add_reaction!



@doc """
    change_objective(model, reaction::String[, objective = 1.0])

Copy the model and change its objective

### Optional arguements
    objective: Sets the objective coefficient. Default is 1.0


Returns the model duplicate with a new objective




    change_objective(model, reaction_index::Number[, objective = 1.0])

Copy the model and change its objective

### Optional arguements
    objective: Sets the objective coefficient. Default is 1.0


Returns the model duplicate with a new objective
""" change_objective

@doc """
    change_objective!(model, reaction::String[, objective = 1.0])

Change the model's objective
### Optional arguements
    objective: Sets the objective coefficient. Default is 1.0



    change_objective!(model, reaction_index::Number[, objective = 1.0])

Change the model's objective
### Optional arguements
    objective: Sets the objective coefficient. Default is 1.0
""" change_objective!



@doc """
    change_reaction_bounds!(model, reaction, lb, ub)

Change the lower \& upper bounds of a reaction

""" change_reaction_bounds!

@doc """
    change_reaction_bounds(model, reaction, lb, ub)

Return a duplicate of the model with the lower/upper bounds of a reaction changed



    change_reaction_bounds(model, reaction, value[, bound])

    change_reaction_bounds(model, reaction_index, value[, bound])

Return a duplicate of the model after changing a bound

### Optional arguements:
    bound:
* 'l' to change the lower bound
* 'u' to change the upper bound
* 'b' to fix lower and upper to the same value (default)
""" change_reaction_bounds


@doc """
	convert_to_irreversible(model)

Finds all reactions that are bi-directional and splits them into one-way reactions 

""" convert_to_irreversible

@doc """
	convert_to_irreversible!(model)

Finds all reactions in a copy of the model that are bi-directional 
 and splits them into one-way reactions 

""" convert_to_irreversible!

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
    reduce_model(model)
 
 
Reduces the model by removing:
* Reactions whose flow is fixed at 0.0
* Metabolites that never appear after cutoff reactions are removed
* Genes that have no effect after cutoff reactions are removed



""" reduce_model!

@doc """
    reduce_model!(model)
 
 
Returns a reduced copy of the model
* Reactions whose flow is fixed at 0.0
* Metabolites that never appear after cutoff reactions are removed
* Genes that have no effect after cutoff reactions are removed
 


""" reduce_model

@doc """
    remove_gene(model, gene::String)
Return a duplicate of the model with a gene removed

    remove_gene(model, gene::Array{String, N})
Return a duplicate of the model with an array of genes removed
""" remove_gene

@doc """
    remove_reaction(model, reaction)
\n
    remove_reaction(model, reaction_index)

Return a duplicate of the model with a reaction removed
""" remove_reaction

@doc """
    remove_reaction!(model, reaction)
\n
    remove_reaction!(model, reaction_index)
Remove a reaction from the model
""" remove_reaction!

@doc """
    remove_reactions(model, reaction)
\n
    remove_reactions(model, reaction_indices)

Return a duplicate of the model with a list of reactions removed
""" remove_reactions

@doc """
    remove_reactions!(model, reaction)
\n
    remove_reactions!(model, reaction_indices)
Remove a list of reactions from the model
""" remove_reactions!


@doc """
    find_blocked_reactions(model, [tolerance = 1e-10])

Return a list of reactions that exhibit zero flux
""" find_blocked_reactions

@doc """
    find_deadend_metabolites(model, lbfilter = true)

Locate all metabolites that are only either produced or consumed but 
not both

### Optional arguements:
    lbfilter: true / false
Filter the result for metabolites who participate in reactions that have `positive` lower bounds
""" find_deadend_metabolites

@doc """
    find_essential_genes(model, min_growth)

Find those genes which, if knocked out, would result in a biomass-production below that 
specified by **min_growth**

### Examples
    julia> find_essential_genes(model,0.01)
    2-element Array{String,1}:
     "b2416"
     "b2415"

Note that **min_growth** represents the actual growth, not a percentage. For min growth of 
0.1%, do

    julia> find_essential_genes(model,0.001 * fba(model).obj)
    2-element Array{String,1}:
     "b2416"
     "b2415"
""" find_essential_genes


@doc """
    find_exchange_reactions(model[, output = 'index')

Find all exhange reactions in the model

### Optional arguements: 
    output:
* 'index': return an array of indices of the exchange reactions
* 'name' : return an array of the names of the exchange reactions

""" find_exchange_reactions

@doc """
    find_reactions_from_genes(model, gene [, output = 'index'])

Find the reactions in which a specific ***gene*** participates

### Optional arguements:
    output:
* 'index': return an array of indices of the reactions
* 'name' : return an array of the names of the reactions

""" find_reactions_from_gene

@doc """
    find_reactions_from_metabolites(model, metabolite [, output = 'index'])

Find the reactions in which a specific ***metabolite*** participates

### Optional arguements:
    output:
* 'index': return an array of indices of the reactions
* 'name' : return an array of the names of the reactions

""" find_reactions_from_metabolite

@doc """
    gene_deletion(model, num_genes)

Perform fba for all gene combinations knock-outs from 1:n 

### Examples
    knockouts = gene_deletion(model,2)

Return a GeneDeletion object which has the following fields:  
*r_g - Reaction/Gene which displays the reactions disabled by gene combinations
*r_f - Reaction/Flow which displays the flow after disabling reactions 
*g_f - Gene/Flow which displays the flow for each gene combination knockout 
""" gene_deletion



@doc """
    knockout_genes(model, genes)

Perform fba after a knocking out selected genes genes    

Return a GeneKnockout type with fields ``growth``, ``flux``, ``disabled`` and ``affected``

### Example

    julia> knockout_genes(model, [model.genes[4], model.genes[16], model.genes[99]])

    Type: GeneKnockout
                            Fields 
                            growth :                      0.374230
                              flux :              Array{Float64,1}
                          disabled :                       [12,71]
                          affected ::
                                                   gene:      reactions:
                                              gene_016 :            [12]
                                              gene_004 :             [3]
                                              gene_099 :            [71]

So knocking out genes 4, 16 and 99, would **touch** reactions 12, 3 and 71, respectively, but only reactions 12 and 71 would be disabled                                              

""" knockout_genes

@doc """
    knockout_reactions(model, reactions)

Attempt to find the simplest combination of genes that need to be knocked out to disable the selected reactions 
""" knockout_reactions

@doc """
    print_reaction_formula(model,reaction::String)
    print_reaction_formula(model,index::Number)

Print out the formula for a reaction
""" print_reaction_formula


@doc """
    solution = robustness_analysis(model, reactions; objective[, pts = [], direction = "max"])

Perform a robustness analysis for any number of fixed reactions, specified by indices or reaction names 

**Example**
To perform a robustness analysis on reaction 13 againts reactions 5,8 and 11, where the point resolution
for reactions 5,8 and 11 is 4,2 and 10, respectively

    robustness_analysis(model, [5,8,11], 13, [4,2,10], "max")
    Robustness Analysis 
                          result :   (2,3,4) Array 
                          ranges : 
                                  reaction :          range: 
                                         5 :     (-0.0,20.0) 
                                         8 :      (0.0,20.0) 
                                        11 :    (8.39,175.0) 



The method returns a **Robustness** object with fields 

 * **solution.result**

    a multidimensional matrix which contains the flux value of the objective for each point
 
 * **solution.reactions**

    names of the reactions

 * **solution.ranges**

    the flux values of the reactions at specific points 

        solution.ranges[1]
        4-element Array{Float64,1}:
         -1.3884e-28
          6.66667   
         13.3333    
         20.0 

    So reaction 1 (reaction 5, ACONTb) is fixed at 6.6666 at point 2 

To see the value of the objective, reaction 13, when reaction 5 
is fixed at point 1, reaction 8 at point 2 and reaction 11 at 
point 10

    solution[1,2,10]
    6.365820260756199e-16

To see every point where the objective is between 0.5 and 0.9
    
    range(solution, 0.5, 0.9)
    4-element Array{Tuple{Int64,Int64,Int64},1}:
     (1,2,3)
     (2,2,3)
     (3,2,3)
     (4,2,3)

To get the entire vector for reaction 5 when reaction 8 is fixed at 
point 2 and reaction 11 is fixed at point 6

    solution[:,2,6]
    4-element Array{Float64,1}:
     0.425313 
     0.314254 
     0.203194 
     0.0921351  

Return a matrix with the objective reactions optimal value at each  
value of the control reaction

**Plotting**

It is reccomended to use PyPlot to plot in Julia.

using PyPlot, to plot a 2D image of the effect of reaction 5
while reaction 8 is fixed at point 2 and reaction 11 at point 6

    PyPlot.plot(a[:,2,6])

to plot a 3D image in PyPlot for reactions 5 and 8 while 
reaction 11 is fixed at point 6

    PyPlot.surf(a[:,:,6])
""" robustness_analysis

@doc """
    synthetic_lethal_genes(model, cutoff = 0.1, num_runs = 100)

Check the model for essential genes, conditionally essential genes and non-essential genes

Returns a ``SLG`` type, which contains the fields ``ess``, ``cond_ess`` and ``non_ess``,::

* ``cutoff`` represents the minimum biomass flux as a fraction of the wild-type flux. 

* ``num_runs`` indicates how many times the algorithm runs, higher number gives better results, but takes longer.

**Example**

To find essentiality with biomass fixed at ``0.1`` ::

    julia> slg = synthetic_lethal_genes(model, 0.1, 200)
    run 200 of 200 Number of conditionally essential genes found: 115

    Type: SLG
         Field                  Contains    Size  
           ess                Essential:       7  
      cond_ess    Conditional Essential:     115  
       non_ess            Non-Essential:      15 
""" synthetic_lethal_genes



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
