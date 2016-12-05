# returns the following:
# rules, new_genes, gene_lookup, effects, rxn_genes
##################
# rules:
# rules read fx "x[1] && x[12]" where x[n] is a gene
##################
# new_genes 
# new_genes are simply x[1] ... x[num_genes]
##################
# gene_lookup
# dictionary with all original gene names AND new_genes as keys and vals 
# so the dict contains both 
# b0001 => x[1] 
# and
# x[1] => b0001
##################
# effects
# dict saying which genes affect which reactions
# 1 => [1,2,3] if gene number 1 affects reactions 1,2 and 3
##################
# rxn_genes 
# dict which includes the genes that appear in every reaction 
# 55 => [3,71] if the reaction rule for reaction 55 has genes 3 and 71 
function prep_genes(model)
    # abstract genes 
    genes = deepcopy(model.genes)

    # RETURN generate vector of x[1] .. x[num_genes]
    num_genes = length(genes) 
    num_rxns = length(model.c) 

    new_genes = Array(String, num_genes)
    for (i,v) in enumerate(genes)
        new_genes[i] = string("x[",i, "]")
    end

    # RETURN dictionary to easily switch between
    gene_lookup = Dict(zip([new_genes; genes], [genes; new_genes]))

    # RETURN abstract rules 
    rules = deepcopy(model.rxn_rules)

    for (i,rule) in enumerate(rules)
        for element in split(rule)
            if haskey(gene_lookup, element)
                rules[i] = replace(rules[i], element, gene_lookup[element])
            end 
        end 
    end 

    # RETURN the effect every gene makes, key => val is gene number and reactions touched
    effects = []
    for gene in new_genes 
        push!(effects, [])

        for (i,rule) in enumerate(rules)
            if contains(rule, gene)
                push!(effects[end], i)
            end 
        end 
    end 
    effects = convert(Array{Array{Int64,1},1}, effects)
    effects = Dict(zip(1:num_genes, effects))

    # RETURN the genes present in every reaction 
    rxn_genes = []
    for i in 1:num_rxns
        push!(rxn_genes, [])
    end 

    for (key,val) in effects
        for v in val
            push!(rxn_genes[v], key)
        end 
    end 
    rxn_genes = convert(Array{Array{Int64,1},1}, rxn_genes)
    rxn_genes = Dict(zip(1:num_rxns, rxn_genes))

    return rules, new_genes, gene_lookup, effects, rxn_genes
end 

# WILL BE DEPRECATED, use prep_genes()
# Use this to convert gene names in model.genes & model.rxn_rules 
# into a 'safe' format gene_xxxx
# This is due to some genes being called b1 and b12, so b1 
# partially matches b12 in rules and results in errors
# returns new_gene_names, new_rules, new_old_genes, saved_rules  
# new_gene_names is a list with reactions named gene_xx1, gene_xx2 ... gene_NNN
# new_rules is a list with all rules, with genes converted to those matching
# the genes in 'new gene names'
# new_old_genes is a dict with a mapping for new_gene_names => original_gene_names
# saved_rules is just the old reaction rules
function abstract_gene_names(genes, rules)
    num_genes = length(genes)
    num_digits = ndigits(num_genes)
    num_rules = length(rules)

    new_gene_names = Array(String, num_genes)
    new_rules = deepcopy(rules)


    for (i,v) in enumerate(genes)
        new_gene_names[i] = string("gene_","0"^(num_digits - ndigits(i)), i)
    end

    saved_names = Dict(zip(genes, new_gene_names))

    # Make sure no parenthesis touches a gene name
    for (i,v) in enumerate(new_rules)
        new_rules[i] = replace(new_rules[i], ")", " ) ")
        new_rules[i] = replace(new_rules[i], "(", " ( ")
    end 

    for (i,rule) in enumerate(new_rules)
        rule = split(rule)

        for (j,v) in enumerate(rule )
            if haskey(saved_names, v)
                rule[j] = saved_names[v]
            end
        end 

        rule = join(rule, " ")
        new_rules[i] = rule 
    end 

    saved_rules = Dict(zip(1:num_rules, rules))


    return new_gene_names, new_rules, invert_dict(saved_names), saved_rules  
end 

function abstract_gene_names(model)
    return abstract_gene_names(model.genes, model.rxn_rules)
end 

# Returns the names of all genes that are in a reaction rule
# Note that rule needs to be converted to abstract
function genes_in_reaction(reaction_rule::String)
    rule = split(reaction_rule)
    indices = find(map(x -> startswith(x, "gene_"), rule))
    genes = rule[indices]

    return unique(genes)
end

# fetch the genes that appear in a reaction rule
function genes_in_reaction_rule(reaction)
    tmp = deepcopy(reaction)
    tmp = replace(tmp, '(', " ")
    tmp = replace(tmp, ')', " ")
    tmp = replace(tmp, '|', " ")
    tmp = replace(tmp, '&', " ")
    tmp = replace(tmp, "false", " ")
    tmp = replace(tmp, r"#=\b\S+?\b=#", "")
    tmp = split(tmp)
end

# set a gene to 'false' in a rule. Note use along with abstract_gene_names
function disable_genes_in_rules{T<:AbstractString}(gene::AbstractString, rules::Array{T})
    for (i,v) in enumerate(rules)
        rules[i] = replace(rules[i], gene, "false")
    end 
    return
end 

function disable_genes_in_rules{T<:AbstractString}(genes::Array, rules::Array{T})
    map(x -> disable_genes_in_rules(x, rules), genes)
    return
end

# set all genes to 'true'. Note use along with abstract_gene_names
function enable_genes_in_rules(rules)
    for (idx, rule) in enumerate(rules), val in split(rule)
        if startswith(val, "gene_")
            rules[idx] = replace(rules[idx], val, "true")
        end 
    end 
    return 
end 

# wrapper for disable_genes_in_rules and enable_genes_in_rules
function convert_rule_to_boolean{T <: AbstractString}(disabled_genes, rules::Array{T})
    disable_genes_in_rules(disabled_genes, rules)
    enable_genes_in_rules(rules)
end

# used to fetch old gene names from a list of new gene names
function convert_from_abstract{T<:String}(dict, gene_list::Array{T})
    for (i,v) in enumerate(gene_list)
        gene_list[i] = dict[v]
    end 
    return 
end 

# find the location of a gene in model.genes from its string name
function find_gene(model::Model, gene::String)
    if gene_exists(model, gene)
        return find(gene .== model.genes)
    else
        warn("Gene does not exist")
        return
    end
end

# check if gene exists in model.genes
function gene_exists(model::Model, gene::String)
    any(gene .== model.genes)
end

# search all reaction gene rules for a gene
# Returns a list of indices corresponding to the reactions which a gene/list of genes
# appears in
function find_gene_effects(rules::Array{String}, gene::String)
    find(x -> contains(x, gene), rules)
end 

function find_gene_effects(rules::Array{String}, genes::Array{String})
    map(x -> find_gene_effects(rules, x), genes)
end 

function find_gene_effects(model::Model, gene::String)
    genes, rules, new_old, _ = abstract_gene_names(model)
    old_new = invert_dict(new_old)
    find(x -> contains(x, old_new[gene]), rules)
end

function find_gene_effects{T <: String}(model::Model, genes::Array{T})
    return map(x -> find_gene_effects(model, x), genes)
end

function find_gene_effects(model::Model, gene_indices::Array{Number})
    if (any(gene_indices .<= 0)) | (any(gene_indices .> length(model.genes)))
        warn("One or more invalid gene indices")
    end 

    return find_gene_effects(model, model.genes[gene_indices])
end

function find_gene_effects(model::Model, gene_index::Number)
    return find_gene_effects(model, [gene_index])
end



# change this to whichever solver you want to use 
# Note: only edit the global variable value
if !isdefined(:settings_default_lp)
    global settings_default_lp = "glpk"
end

# go through a list of blocked reactions and 
# solve the lp, and restore
# sorts the blocked reactions lexiographically
# to minimize solver calls
function analyze_list_of_blocked_reaction_lp(lp, rxns)
    progress_meter = ProgressMeter.Progress(length(rxns), 0.5, "\n\u1b[1FFVA:", 50)
    combos = sort(rxns, lt=lexless)
    num_combos = length(combos)

    wt = answer_lp(lp)

    lb, ub = get_variable_bounds(lp)

    flow = Array(Float64, num_combos)

    ####
    for d in combos[1] 
        set_col_bounds(lp, d, 0.0, 0.0)
    end 

    flow[1] = answer_lp(lp)

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
        
        flow[i] = answer_lp(lp)
        next!(progress_meter)
    end

    for e in disable 
        set_col_bounds(lp, e, lb[e], ub[e])
    end     

    flow /= wt
    flow = round(flow, 16)

    result = Dict(zip(combos, flow))

    return result
end 

# Solve lp and return:
# the objective value if optimial solution is found
# 0.0 otherwise
function answer_lp(lp)
    solve(lp)
    if get_solution_status(lp) == "Optimal"
        return get_objective_value(lp)
    else 
        return 0.0
    end 
end

# return a FBAsolution. lp must be solved
function construct_lp_solution(lp)
    f = get_objective_value(lp)
    x = get_solution(lp)
    y = get_slack(lp)
    w = get_reduced_costs(lp)
    success = get_solution_status(lp) 
    info = SolverInfo(settings_default_lp, solver_status(lp))

    return  FBAsolution(f,x,y,w,success, info);
end

function createlpmodel()
  # Create a LP model structure and set solver specific parameters
  # THINK: Check memory allocation/deallocation
  # THINK: Support user-defined parameters
  global cobra_lp_solver
  if cobra_lp_solver == "cplex"
    # Use Primal Simplex
    m = JuMP.Model(solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_LPMETHOD=1))
  elseif cobra_lp_solver == "gurobi"
    m = JuMP.Model(solver=GurobiSolver(OutputFlag=0, Method=0))
  elseif cobra_lp_solver == "glpk"
    m = JuMP.Model(solver=GLPKSolverLP()) # THINK: Fails!!!!
  else
    error("Solver ", cobra_lp_solver, "not implemented yet.")
  end
  m
end

function createmilpmodel()
  # Create a MILP model structure and set solver specific parameters
  # THINK: Support user-defined parameters
  global cobra_mip_solver
  if cobra_mip_solver == "cplex"
    # CPX_PARAM_SCRIND=0
    m = JuMP.Model(solver=CplexSolver())
  elseif cobra_mip_solver == "gurobi"
    #OutputFlag=0
    m = JuMP.Model(solver=GurobiSolver())
  elseif cobra_mip_solver == "glpk"
    m = JuMP.Model(solver=GLPKSolverMIP()) # THINK: Fails!!!!
  else
    error("Solver ", cobra_lp_solver, "not implemented yet.")
  end
  m
end

function createqpmodel()
  global cobra_qp_solver
  if cobra_qp_solver == "cplex"
    # Use Primal Simplex for QPs
    m = JuMP.Model(solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_QPMETHOD=1))
  elseif cobra_qp_solver == "gurobi"
    assert("Not tested yet!")
    m = JuMP.Model(solver=GurobiSolver(OutputFlag=0, Method=0))
  elseif cobra_qp_solver == "glpk"
    error("GLPK cannot solve quadratic problems.")
  else
    error("Solver ", cobra_qp_solver, "not implemented yet.")
  end
  m
end

# change the default solver 
function set_solver(solver)
  Solver = ""
  if solver == "glpk"
    Solver = "GLPK"
  elseif solver == "cplex"
    Solver = "CPLEX"
  elseif solver == "gurobi"
    Solver = "Gurobi"
  elseif solver == "clp"
    Solver = "Clp"
  else 
    warn("choose one of the following solvers:")
    warn("glpk")
    warn("cplex")
    warn("gurobi")
    warn("clp")
  return
  end 


  if !any(readdir(Pkg.dir()) .== Solver)
    warn(Solver * " is not installed")
    warn("Please run: Pkg.add(" * Solver * ") and restart Julia")
    return 
  end 

  global settings_default_lp = solver
end

# setup an lp model, GLPK, CPLEX, Gurobi or Clp depending on that is set as the default solver
function setup_lp(model::Model, objective::String = "max", solver::String = "")
	if solver == ""
		solver = settings_default_lp
	end 
	if solver == "glpk"
		return setup_glpk(model, objective)
	elseif solver == "gurobi"
		return setup_gurobi(model, objective)
	elseif solver == "clp"
		return setup_clp(model, objective)
	elseif solver == "cplex"
		# Should be an elif, else should trigger an error (unknown solver) 
		return setup_cplex(model, objective)
	else 
		warn("Unknown solver. Check settings.jl for the return value of settings_lp_solver()")
	end
end

# THINK: This function should not be structured differently from the above one
function setup_lp{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, S::SparseMatrixCSC, objective::String = "max", solver::String = "glpk")
	if solver == "CPLEX"
		return setup_cplex(c, lb, ub, b, S, objective)
	elseif solver == "Gurobi"
		return setup_gurobi(c, lb, ub, b, S, objective)
	else
		return setup_glpk(c, lb, ub, b, S, objective)
	end
end

# find the minimum and maximum values for a chosen objective 
function min_max_lp(lp, objective::Number = 0, value::Number = 1.0)

	if (objective != 0) & isinteger(objective)
		change_objective_coef(lp, objective, value)
	end

	set_direction(lp, "min")
	solve(lp)
	min = get_objective_value(lp)

	set_direction(lp, "max")
	solve(lp)
	max = get_objective_value(lp)

	return min,max
end


# Delete a column from a sparse matrix 
function delete_column(struct::AbstractSparseMatrix, n::Number)
	reaction_indices = collect(nzrange(struct, n))
	sparse_size = length(struct.nzval)
	reactions = get_sparse_columns(struct)
	
	new_j = deleteat!(deepcopy(struct.rowval), reaction_indices)
	new_i = deleteat!(reactions, reaction_indices)
	new_a = deleteat!(deepcopy(struct.nzval), reaction_indices)

	for i in reaction_indices[1]:length(new_i)
		new_i[i] -= 1
	end

	j = convert(Array{Int64}, new_j)
	i = convert(Array{Int64}, new_i)
	a = convert(Array{Float64}, new_a)

	return sparse(j,i,a)
end

# add a column to a sparse matrix 
function add_column(struct::AbstractSparseMatrix)
	struct = hcat(struct,sparse(zeros(struct.m)))
end

# add a row to a sparse matrix 
function add_row(struct::AbstractSparseMatrix)
	struct = vcat(struct,sparse(zeros(struct.n)'))
end

function get_sparse_columns(struct::AbstractSparseMatrix)
	reactions = []

	for i in 1:struct.n
		push!(reactions, i*ones(Int64, length(nzrange(struct, i))))
	end
	return vcat(reactions...)
end

function statistics(model::Model, verbose::Bool = false)
	# THINK: Include comment/name field
	rxns = length(model.rxns)
	mets = length(model.mets)
	genes = length(model.genes)
	if verbose
		@printf("Reactions: %d\nMetabolites: %d\nGenes: %d\n", rxns, mets, genes)
	end
	rxns, mets, genes
end

# Fill up empty model fields with netural data 
# to keep all reaction/metabolite/gene fields the 
# same length
function fix_empty_fields(model)
    nrxns = length(model.rxns)
    nmets = length(model.mets)
    ngenes = length(model.genes)


    if length(model.rxn_name) == 0
        model.rxn_name = fill("", nrxns)
    end

    if length(model.rxn_rules) == 0
        model.rxn_rules = fill("", nrxns)
    end

    if length(model.rxn_subsystem) == 0
        model.rxn_subsystem = fill("", nrxns)
    end

    if length(model.met_formula) == 0
        model.met_formula = fill("", nmets)
    end

    if length(model.met_name) == 0
        model.met_name = fill("", nmets)
    end

    if length(model.gene_name) == 0
        model.gene_name = fill("", ngenes)
    end   
end

# wrapper for all fixing functions 
function fix_model(model)
    fix_empty_fields(model)
    fix_gene_rules(model)
    fix_genes(model)
    fix_extra_dicts(model)
end 

# converts the arrays in rxn_extra/met_extra/gene_extra 
# dicitonaries to Array{Any}
function fix_extra_dicts(model)

    for (key,value) in model.rxn_extra
        model.rxn_extra[key] = convert(Array{Any}, value)
    end 

    for (key,value) in model.met_extra
        model.met_extra[key] = convert(Array{Any}, value)
    end 

    for (key,value) in model.gene_extra
        model.gene_extra[key] = convert(Array{Any}, value)
    end 
end 

# models may have odd string errors, such as ill-detectable 
# whitespaces, and/or mushed agains parentheses or genes 
function fix_gene_rules(model)
    model.rxn_rules = convert(Array{String}, model.rxn_rules)
    model.genes = convert(Array{String}, model.genes)

    
    # take care of whitespaces
    model.rxn_rules = split.(model.rxn_rules)
    model.rxn_rules =  map(z -> map(x -> map(y -> isspace(y) ? ' ' : y, x), z), model.rxn_rules)
    model.rxn_rules = map(x -> rstrip.(x), model.rxn_rules)
    model.rxn_rules = map(x -> lstrip.(x), model.rxn_rules)
    model.rxn_rules = map(x -> join(x, " "), model.rxn_rules)
    
    
    # detect where "and" or "or" are mushed to genes
    # some gene gene100 and gene10 will disrupt
    # each other so I sort them so the longest 
    # genes are checked first
    # i.e if a rule reads "gene100and gene12"
    # and genes gene100 and gene10 exist, what happens
    # that if gene10 is checked I will do 
    # replace(rule, "gene10", " gene10 ")
    # which results in the rule to read
    # "gene10 0and gene12"
    num_genes = length(model.genes)
    gene_name_lengths = sort(unique(length.(model.genes)), rev=true)
    gene_idx_by_name_length = map(y -> find(map(x -> length(x) == y, model.genes)), gene_name_lengths)
    gene_idx_by_name_length = vcat(gene_idx_by_name_length...)
    length_sorted_genes = model.genes[gene_idx_by_name_length]
    genes = length_sorted_genes
    num_digits = ndigits(num_genes)

    new_gene_names = Array(String, num_genes)

    for (i,v) in enumerate(genes)
        new_gene_names[i] = string("gene_","0"^(num_digits - ndigits(i)), i)
    end

    saved_names = Dict(zip(genes, new_gene_names))
    rules = deepcopy(model.rxn_rules)
    for (i,gene) in enumerate(genes) 
        for (j,v) in enumerate(rules)
            if contains(v, gene)
                rules[j] = replace(rules[j], gene, saved_names[gene])
            end 
        end 
    end 

    rules = map(x -> replace(x, "(", " ( "), rules)
    rules = map(x -> replace(x, ")", " ) "), rules)

    rules = map(x -> replace(x, "and", " and "), rules)
    rules = map(x -> replace(x, "&&", " && "), rules)
    rules = map(x -> replace(x, "and", "&&"), rules)

    rules = map(x -> replace(x, "or", " or "), rules)
    rules = map(x -> replace(x, "||", " || "), rules)
    rules = map(x -> replace(x, "or", "||"), rules)

    saved_names = invert_dict(saved_names)

    for (key,val) in saved_names
        for (i, rule) in enumerate(rules)
            rules[i] = replace(rules[i], key, val)
        end 
    end 

    model.rxn_rules = rules

    return
end 

# remove genes that dont appear in any rules 
function fix_genes(model)
    # remove genes that dont appear in any rules 
    never_appear = []

    for (i,gene) in enumerate(model.genes)
        if isempty(find(x -> contains(x, gene), model.rxn_rules))
            push!(never_appear, i)
        end 
    end 

    deleteat!(model.genes, never_appear)
    deleteat!(model.gene_name, never_appear)


    for (key,val) in model.gene_extra 
        deleteat!(model.gene_extra[key], never_appear)
    end 
end 
function sendto(p::Int; args...)
    for (nm, val) in args        
        @spawnat(p, eval(Main, Expr(:(=), nm, val)))
    end
end

function sendto(ps::Vector{Int}; args...)
    for p in ps
        sendto(p; args...)
    end
end

function disable_cores()
    rmprocs(procs()[2:end])
end 

function initialize_cores()
    addprocs(Sys.CPU_CORES-1)
    extra_cores = procs()[2:end]

    @sync for id in extra_cores
        @async remotecall_fetch(include, id, Pkg.dir() * "/CBM/src/CBM.jl")
    end 
end 

function initcores()
    addprocs(Sys.CPU_CORES-1)
    extra_cores = procs()

    include(Pkg.dir() * "/CBM/src/core/multicore.jl")
end 


function clone(model)
	return deepcopy(model)
end

function find_metabolite(model::Model, metabolite::String)
	if metabolite_exists(model, metabolite)
		return findin(model.mets,[metabolite])[1]
	else
		return
	end
end

function metabolite_exists(model::Model, metabolite::String)
	any(model.mets .== metabolite)
end

function find_reaction{T<:String}(model::Model, reaction::T)
	index_rxns = find(x -> lowercase(x) == lowercase(reaction), model.rxns)
	index_rxn_name = find(x -> lowercase(x) == lowercase(reaction), model.rxn_name)

	return vcat([index_rxns; index_rxn_name]...)
end

function find_reaction{T<:String}(model::Model, reactions::Array{T})
	indices = map(x -> find_reaction(model,x), reactions)
	return indices
end

function reaction_exists(model::Model, reaction::String)
	any(model.rxns .== reaction)
end 

function find_text(text)
	for folder in folders, file in readdir(folder)
		file = pwd() * "/" * folder * "/" * file
		f = open(file)
		try
			for (line_number, line) in enumerate(eachline(f))
				if contains(line, text)
					place = collect(search(line, text))
					@printf "%50s\x1b[31m\x1b[1m ::%4s:: \x1b[0m%s\x1b[32m\x1b[1m%s\x1b[0m%s" file line_number line[1:place[1]-1] line[place] line[place[end]+1:end]
				end
			end
		catch
			println("*** Read of "*file*" failed")
		end
	end
end

function invert_dict(dict, warning::Bool = false)
	vals = collect(values(dict))
	dict_length = length(unique(vals))

	if dict_length < length(dict)
		if warning
			warn("Keys/Vals are not one-to-one")
		end 

		linked_list = Array[]

		for i in vals 
			push!(linked_list,[])
		end 

		new_dict = Dict(zip(vals, linked_list))

		for (key,val) in dict 
			push!(new_dict[val],key)
		end
	else
		key = collect(keys(dict))

		counter = 0
		for (k,v) in dict 
			counter += 1
			vals[counter] = v
			key[counter] = k
		end
		new_dict = Dict(zip(vals, key))
	end 

	return new_dict
end


"""
Wrapper for `load_json()` and `load_matlab()`

    model = load_model(filename)
"""
function load_model(filename::String; fix = true)
    all_files = readdir(dirname(filename))

    if endswith(filename, ".json")
        return load_json(filename)
    elseif endswith(filename, ".mat")
        return load_matlab(filename)
    elseif any(all_files .== basename(filename) * ".json")
        return load_model(filename * ".json")
    elseif any(all_files .== basename(filename) * ".mat")
        return load_model(filename * ".mat")
    else 
        warn("Model not found")
    end 
end 

"""
Returns a dictionary with variables saved in a `.mat` file 

    open_mat_file(filename)

**Example**

    open_mat_file("/home/user/variables.mat")
     Dict{String,Any} with 2 entries:
     "var_a" => "[1,2,3,4]"
     "var_b" => "hello"
"""
function open_mat_file(filename::String)
    filename = endswith(filename, ".mat") ? filename : filename * ".mat"

    if any(readdir(dirname(filename)) .== basename(filename))
        mat_file = matread(filename)
    else 
        warn("File not found, does it end with '.mat' ? ")
    end 
end 

"""
Use to save variables to `.mat` format 

    save_mat_file(filename; args...)

**Example**

If you have the variables `a` and `b` defined as 

    a = [1,2,3,4]
    b = "hello"

You can save them by calling 

    save_mat_file("/home/user/variables.mat", var_a = a, var_b = b)

and they can be accessed within `Julia` by calling 

    open_mat_file("/home/user/variables.mat")
     Dict{String,Any} with 2 entries:
     "var_a" => "[1,2,3,4]"
     "var_b" => "hello"
"""
function save_mat_file(filename::String; args...)
    filename = endswith(filename, ".mat") ? filename : filename * ".mat"
    varnames = [] 
    varvals = [] 

    for (nm, val) in args 
        push!(varnames, string(nm))       
        push!(varvals, string(val))       
    end


    MAT.matwrite(filename, Dict(zip(varnames, varvals)))    
end 
