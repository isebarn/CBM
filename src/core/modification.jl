"""
    add_metabolite!(model, metabolite)

Add a metabolite to the model
"""
function add_metabolite!(model::Model, metabolite::String; b::Number = 0, csense::String = "=")
	add_metabolite!(model, metabolite, b, csense)
end 

function add_metabolite!(model::Model, metabolite::String, b::Number = 0, csense::String = "=")
	if !metabolite_exists(model, metabolite)
		push!(model.mets, metabolite)
		push!(model.met_name, metabolite)
		push!(model.met_formula, "")
		push!(model.b, b)
		model.S = add_row(model.S)
		push!(model.csense[csense], length(model.mets))

		for (key,value) in model.met_extra
			push!(model.met_extra[key], nothing)
		end 

	else
		warn(metabolite, " already exists in the model")
	end
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


"""
    add_reaction(model,  reaction_name, reaction, [lb == -1000, ub = 1000])

Add a reaction to a copy of the model

### Optional arguements:
    lb \& ub
* add lower/upper bounds, lower bound must be less/equal to upper bound

Return the new model with the added reaction.

*Note: Original model remains unchanged...see add_reaction! to permanently add a reaction*
"""
function add_reaction(model::Model, rxn_name::String, reaction::String; lb::Number = -1000.0, ub::Number = 1000.0)
	return add_reaction(model, rxn_name, reaction, lb, ub)
end 

function add_reaction(model::Model, rxn_name::String, reaction::String, lb::Number = -1000.0, ub::Number = 1000.0)
	tmp = deepcopy(model)
	tmp = add_reaction!(tmp, rxn_name, reaction, lb, ub)
end

"""
    add_reaction!(model, reaction_name, reaction, [lb == -1000, ub = 1000])

Permanently add a reaction to the model

### Optional arguements:
    lb \& ub
* add lower/upper bounds, lower bound must be less/equal to upper bound
"""
function add_reaction!(model::Model, rxn_name::String, reaction::String; lb::Number = -1000.0, ub::Number = 1000.0)
	return add_reaction!(model, rxn_name, reaction, lb, ub)
end 

function add_reaction!(model::Model, rxn_name::String, reaction::String, lb::Number = -1000.0, ub::Number = 1000.0)
    if contains(reaction, "<=>")
        middle = "<=>"
    elseif contains(reaction, "<=")
        middle = "<="
    elseif contains(reaction, "=>")
        middle = "=>"
    else 
        warn("Reaction must contain <=>, <= or =>")
        return 
    end 

    # use this for faster lookup and its less "messy" to use than
    # to use find()
    metdict = Dict(zip(model.mets, collect(1:length(model.mets))))

    model.S = add_column(model.S)
    col = model.S.n 

    reaction = replace(reaction, " ", "")
    sides = split(reaction, middle);
    left = sides[1];
    right = sides[2];

    if length(left) != 0
        left_mets = split.(split(left, " + "))
        left_mets = map(y -> filter(x -> x != "", y), left_mets)
        left_mets = map(x -> convert(Array{Any}, x), left_mets)
        map(x -> length(x) == 1 ? unshift!(x, "1.0") : x, left_mets)
        for x in left_mets 
            met = convert(String, x[2])
            if !haskey(metdict, met)
                add_metabolite!(model, met)
                metdict[met] = length(model.mets)
            end 
	    end 
        left_mets = map(x -> [-parse(x[1]), metdict[x[2]]], left_mets)

        for i in left_mets 
            model.S[Int(i[2]), col] = i[1]
        end 
    end 

    if length(right) != 0
        right_mets = split.(split(right, " + "))
        right_mets = map(y -> filter(x -> x != "", y), right_mets)
        right_mets = map(x -> convert(Array{Any}, x), right_mets)
	    map(x -> length(x) == 1 ? unshift!(x, "1.0") : x, right_mets)
        for x in right_mets 
            met = convert(String, x[2])
	        if !haskey(metdict, met)
	        	add_metabolite!(model, met)
	        	metdict[met] = length(model.mets)
	        end 
	    end         
        right_mets = map(x -> [parse(x[1]), metdict[x[2]]], right_mets)

        for i in right_mets 
            model.S[Int(i[2]), col] = i[1]
        end 
    end

    if rxn_name == ""
    	rxn_name = string("reaction_", length(model.c)) 
    end 
    
    push!(model.rxns, rxn_name)
    push!(model.lb, lb)
    push!(model.ub, ub)
    push!(model.c, 0.0)
    push!(model.rxn_name, rxn_name)
    push!(model.rxn_rules, "")
    push!(model.rxn_subsystem, "")

    for (key,val) in model.rxn_extra
    	push!(model.rxn_extra[key], "")
    end 

    return model
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------

# single reaction by index. returns new model
"""
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
"""
function change_objective(model::Model, reaction_index::Number, objective::Number = 1.0 )
	tmp = deepcopy(model)
	change_objective!(tmp, [reaction_index], [objective])
	return tmp
end

# single reaction by index. returns new model. named arguements
function change_objective(model::Model, reaction_index::Number; objective::Number = 1.0 )
	tmp = deepcopy(model)
	change_objective!(tmp, [reaction_index], [objective])
	return tmp
end

# single reaction by index. modifies the model
"""
    change_objective!(model, reaction::String[, objective = 1.0])

Change the model's objective
### Optional arguements
    objective: Sets the objective coefficient. Default is 1.0



    change_objective!(model, reaction_index::Number[, objective = 1.0])

Change the model's objective
### Optional arguements
    objective: Sets the objective coefficient. Default is 1.0
"""
function change_objective!(model::Model, reaction_index::Number, objective::Number = 1.0 )
	change_objective!(model, [reaction_index], [objective])
end
# single reaction by index. modifies the model. named arguements
function change_objective!(model::Model, reaction_index::Number; objective::Number = 1.0 )
	change_objective!(model, [reaction_index], [objective])
end

# multiple reactions by indices. returns new model. handles errors
function change_objective{T,V <: Number}(model::Model, reaction_indices::Array{T}, objectives::Array{V} = Number[] )

	if length(objectives) == 0
		objectives = ones(Float64, length(reaction_indices))
	end 

	if length(reaction_indices) != length(objectives)
		warn("You must supply as many reactions as objectives")
		return 
	end 

	if any(x -> (x <= 0) | (x > length(model.c)) , reaction_indices)
		warn("Invalid reaction indices")
		return 
	end 

	if length(unique(reaction_indices)) != length(reaction_indices)
		warn("Reaction indices must be unique")
		return 
	end 

	if !isinteger(reaction_indices)
		warn("Not all indices are integers")
		return 
	end

	# convert reaction_indices in case the user supplied values
	# like [1.0 6.0] instead of [1 6]
	reaction_indices = convert(Array{Int64}, reaction_indices)

	tmp = deepcopy(model)
	change_objective!(tmp, reaction_indices, objectives)
	return tmp
end

# multiple reactions by indices. returns new model. named arguements
function change_objective{T,V <: Number}(model::Model, reaction_indices::Array{T}; objectives::Array{V} = Number[] )
	tmp = change_objective(model, reaction_indices, objectives)
	return tmp
end

# multiple reactions by indices !. handles errors
# This is the 'handler' function
function change_objective!{T,V <: Number}(model::Model, reaction_indices::Array{T}, objectives::Array{V} = Number[] )
	if objectives == []
		objectives = ones(Float64, length(reaction_indices))
	end 

	if length(reaction_indices) != length(objectives)
		warn("You must supply as many reactions as objectives")
		return 
	end 

	if any(x -> (x <= 0) | (x > length(model.c)) , reaction_indices)
		warn("Invalid reaction indices")
		return 
	end 

	if length(unique(reaction_indices)) != length(reaction_indices)
		warn("Reaction indices must be unique")
		return 
	end 

	if !isinteger(reaction_indices)
		warn("Not all indices are integers")
		return 
	end

	model.c = zeros(Float64, length(model.c))
	model.c[reaction_indices] = objectives
end

# multiple reactions by indices. modifies the model. named arguements
function change_objective!{T,V <: Number}(model::Model, reaction_indices::Array{T}; objectives::Array{V} = Number[] )
	change_objective!(model, reaction_indices, objectives)
end

# single reaction by name. returns new model. 
function change_objective(model::Model, reaction_name::String, objective::Number = 1.0)
	reaction_index = find_reaction(model, reaction_name)[1]
	tmp = deepcopy(model)
	change_objective!(tmp, [reaction_index], [objective])
	return tmp
end

# single reaction by name. returns new model. named arguements
function change_objective(model::Model, reaction_name::String; objective::Number = 1.0)
	reaction_index = find_reaction(model, reaction_name)[1]
	tmp = deepcopy(model)
	change_objective!(tmp, [reaction_index], [objective])
	return tmp	
end
# single reaction by name. modifies the model. 
function change_objective!(model::Model, reaction_name::String, objective::Number = 1.0)
	reaction_index = find_reaction(model, reaction_name)[1]
	change_objective!(model, [reaction_index], [objective])
end
# single reaction by name. modifiels the model. named arguements
function change_objective!(model::Model, reaction_name::String; objective::Number = 1.0)
	reaction_index = find_reaction(model, reaction_name)[1]
	change_objective!(model, [reaction_index], [objective])
end

# multiple reactions by names. returns new model 
function change_objective(model::Model, reaction_names::Array{String}, objectives::Array{Number} = Number[])
	reaction_indices = find_reaction(model, reaction_names)
	tmp = deepcopy(model)
	change_objective(tmp, reaction_indices, objectives)
	return tmp
end
# multiple reactions by names. returns new model. named arguements
function change_objective(model::Model, reaction_names::Array{String}; objectives::Array{Number} = Number[])
	reaction_indices = find_reaction(model, reaction_names)
	tmp = deepcopy(model)
	change_objective(tmp, reaction_indices, objectives)
	return tmp
end

# multiple reactions by names. modifies the model
function change_objective!(model::Model, reaction_names::Array{String}, objectives::Array{Number} = Number[])
	reaction_indices = find_reaction(model, reaction_names)
	change_objective!(model, reaction_index, objective)
end
# multiple reactions by names. modifies the model. named arguements
function change_objective!(model::Model, reaction_names::Array{String}; objectives::Array{Number} = Number[])
	reaction_indices = find_reaction(model, reaction_names)
	change_objective!(model, reaction_names, objective)
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


# THINK: "both" has different meaning in Matlab Cobra
# (sets lb and ub to the same value)

function change_reaction_bounds(model::Model, reaction::String, value::Number, bound::String = "both")
	tmp = deepcopy(model)
	change_reaction_bounds!(tmp, reaction, value, bound)
	return tmp
end

"""
    change_reaction_bounds!(model, reaction, lb, ub)

Change the lower \& upper bounds of a reaction

"""
function change_reaction_bounds!(model::Model, reaction::String, value::Number, bound::String = "both")

	location = find_reaction(model, reaction)
	if length(location) == 0
		warn("The reaction was not found in the model")
	end 
	change_reaction_bounds!(model, location[1], value, value, bound)
end

###

function change_reaction_bounds(model::Model, reaction::String, lb::Number, ub::Number, bound::String = "both")

	tmp = deepcopy(model)
	change_reaction_bounds!(tmp, reaction, lb, ub, bound)

	return tmp
end

function change_reaction_bounds!(model::Model, reaction::String, lb::Number, ub::Number, bound::String = "both")

	location = find_reaction(model, reaction)
	if length(location) == 0
		warn("The reaction was not found in the model")
	end 

	change_reaction_bounds!(model, location[1], lb, ub, bound)
end

###

function change_reaction_bounds(model::Model, reaction_index::Number, value::Number, bound::String = "both")
	tmp = deepcopy(model)
	change_reaction_bounds!(tmp, reaction_index, value, bound)
	return tmp	
end 

function change_reaction_bounds!(model::Model, reaction_index::Number, value::Number, bound::String = "both")
	if (reaction_index < 0) | (reaction_index > length(model.c))
		warn("Invalid reaction indices")
		return 
	end 

	if !isinteger(reaction_index)
		warn("Not all indices are integers")
		return 
	end
	
	change_reaction_bounds!(model, reaction_index, value, value, bound)
end 

###
"""
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
"""
function change_reaction_bounds(model::Model, reaction_index::Number, lb::Number, ub::Number, bound::String = "both")
	if (reaction_index < 0) | (reaction_index > length(model.c))
		warn("Invalid reaction indices")
		return 
	end 

	if !isinteger(reaction_index)
		warn("Not all indices are integers")
		return 
	end

	tmp = deepcopy(model)
	change_reaction_bounds!(tmp, reaction_index, lb, ub, bound)
	return tmp
end

function change_reaction_bounds!(model::Model, reaction_index::Number, lb::Number, ub::Number, bound::String = "both")	if length(reaction_index[1]) == 0
		warn("The reaction doesnt exist in the model")
		return 
	end

	if lb > ub
		warn("Lower bound cannot be set greater than upper bound")
		return
	end

	if contains(lowercase(bound), "l")
		model.lb[reaction_index] = convert(Float64, lb)
	elseif contains(lowercase(bound), "u")
		model.ub[reaction_index] = convert(Float64, ub)
	else
		model.lb[reaction_index] = convert(Float64, lb)
		model.ub[reaction_index] = convert(Float64, ub)
	end
end

#-----------------------------------------------------
#-----------------------------------------------------


"""
    convert_to_irreversible(model)

Finds all reactions that are bi-directional and splits them into one-way reactions 

"""
function convert_to_irreversible(model::Model)
	tmp = deepcopy(model)
	return convert_to_irreversible!(tmp)
end

"""
    convert_to_irreversible!(model)

Finds all reactions in a copy of the model that are bi-directional 
 and splits them into one-way reactions 

"""
function convert_to_irreversible!(model::Model)
    reversible_reactions = find_reversible(model)
    new_reactions = collect(1:length(reversible_reactions)) + model.S.n
    for i in reversible_reactions
        model.S = add_column(model.S)
        indexes = nzrange(model.S, i)
        mets = model.S.rowval[indexes]
        coeffs = model.S.nzval[indexes]
        model.S[mets, model.S.n] = coeffs 
    end 

    model.lb = vcat(model.lb, zeros(length(new_reactions)))
    model.ub = vcat(model.ub, deepcopy(model.ub[reversible_reactions]))
    model.ub[reversible_reactions] = 0.0

    model.c = vcat(model.c, model.c[reversible_reactions])

    model.rxns = vcat(model.rxns, deepcopy(model.rxns[reversible_reactions]))
    model.rxns[reversible_reactions] .*= "_left"
    model.rxns[new_reactions] .*= "_right"

    model.rxn_name = vcat(model.rxn_name, deepcopy(model.rxn_name[reversible_reactions]))
    model.rxn_name[reversible_reactions] .*= "_left"
    model.rxn_name[new_reactions] .*= "_right"

    model.rxn_rules = vcat(model.rxn_rules, model.rxn_rules[reversible_reactions])

    model.rxn_subsystem = vcat(model.rxn_subsystem, model.rxn_subsystem[reversible_reactions])

    for (key,val) in model.rxn_extra
        vals = vcat()
        model.rxn_extra[key] = vcat(val, deepcopy(val[reversible_reactions]))
    end 

    return model
end 

function find_reversible(model::Model)
	reversible_reactions = find(x -> (model.lb[x] < 0) & (model.ub[x] > 0), collect(1:length(model.c)))
end 

# -------------------------------------------------------------------
# -------------------------------------------------------------------

function convert_to_reversible!(model::Model)
    duplicates = []
    for i in 1:model.S.n
        mets = model.S.rowval[nzrange(model.S,i)]
        coeffs = model.S.nzval[nzrange(model.S,i)]

        for j in i+1:model.S.n 
            if model.S.rowval[nzrange(model.S,j)] == mets 
                if model.S.nzval[nzrange(model.S,j)] == coeffs
                    push!(duplicates, j)
                    model.ub[i] = deepcopy(model.ub[j])
                end 
            end 
        end 
    end 

    new_indices = setdiff(1:model.S.n, duplicates)

    S = sparse(zeros(model.S.m, length(new_indices)))
    for i in new_indices 
        mets = model.S.rowval[nzrange(model.S,i)]
        coeffs = model.S.nzval[nzrange(model.S,i)]
        S[mets, i] = coeffs
    end 
    model.S = S

    
    deleteat!(model.rxns, duplicates)
    deleteat!(model.lb, duplicates)
    deleteat!(model.ub, duplicates)
    deleteat!(model.c, duplicates)
    deleteat!(model.rxn_name, duplicates)
    deleteat!(model.rxn_rules, duplicates)
    deleteat!(model.rxn_subsystem, duplicates)

    for (key,val) in model.rxn_extra
        model.rxn_extra[key] = val[new_indices]
    end 

    return model
end

function convert_to_reversible(model::Model)
	tmp = deepcopy(model)
	convert_to_reversible(tmp)
end 

# -------------------------------------------------------------------
# -------------------------------------------------------------------


"""
    reduce_model!(model)
 
 
Returns a reduced copy of the model
* Reactions whose flow is fixed at 0.0
* Metabolites that never appear after cutoff reactions are removed
* Genes that have no effect after cutoff reactions are removed
 


"""
function reduce_model(model::Model)
    tmp = deepcopy(model)
    reduce_model!(tmp) 
end

"""
    reduce_model(model)
 
 Reduces the model by removing:
* Reactions whose flow is fixed at 0.0
* Metabolites that never appear after cutoff reactions are removed
* Genes that have no effect after cutoff reactions are removed

"""
function reduce_model!(model::Model)
    reference = fba(model, "max").x

    min_flux, max_flux = fva(model, 0)

    useless = find(abs(min_flux - max_flux) .< 1e-3)
    useless = useless[find(abs(reference[useless]) .< 1e-3)]

    new_indices = setdiff(1:model.S.n, useless)

    S = sparse(zeros(model.S.m, length(new_indices)))
    for (i,v) in enumerate(new_indices) 
        mets = model.S.rowval[nzrange(model.S,v)]
        coeffs = model.S.nzval[nzrange(model.S,v)]
        S[mets, i] = coeffs
    end 
    model.S = S

    
    deleteat!(model.rxns, useless)
    deleteat!(model.lb, useless)
    deleteat!(model.ub, useless)
    deleteat!(model.c, useless)
    deleteat!(model.rxn_name, useless)
    deleteat!(model.rxn_rules, useless)
    deleteat!(model.rxn_subsystem, useless)

    for (key,val) in model.rxn_extra
        model.rxn_extra[key] = val[new_indices]
    end 

    return model
end 

# -------------------------------------------------------------------
# -------------------------------------------------------------------

"""
    remove_gene(model, gene::String)
Return a duplicate of the model with a gene removed

    remove_gene(model, gene::Array{String, N})
Return a duplicate of the model with an array of genes removed
""" 
function remove_gene(model::Model, gene_index::Number)
	
	tmp = deepcopy(model)
	remove_gene!(tmp, gene_index)
	
	return tmp
end

function remove_gene!(model::Model, gene_index::Number)
	
	if (gene_index <= 0) | (gene_index > length(model.genes))
		warn("Invalid gene index")
		return 
	end 
	gene_names = convert(Array{String}, [model.genes[gene_index]])
	remove_gene!(model, gene_names)
end

function remove_gene(model::Model, gene_name::String)
	tmp = deepcopy(model)
	remove_gene!(tmp, gene_name)
	return tmp
end

function remove_gene!(model::Model, gene_name::String)
	if !in(gene_name, model.genes)
		warn("Gene does not exist")
		return
	end 
	gene_names = convert(Array{String}, [gene_name])
	remove_gene!(model, gene_names)
end

function remove_gene{T <: Number}(model::Model, gene_indices::Array{T})
	tmp = deepcopy(model)
	remove_gene!(tmp, gene_indices)
	return tmp
end

function remove_gene!{T <: Number}(model::Model, gene_indices::Array{T})
	gene_indices = unique(gene_indices)
	if !isempty(find(x -> (x <= 0) | (x > length(model.genes)), gene_indices))
		warn("One or more indices are invalid")
		return
	end 

	gene_names = model.genes[gene_indices]
	gene_names = convert(Array{String}, gene_names)

	remove_gene!(model, gene_names)
end


function remove_gene{T <: String}(model::Model, gene_names::Array{T})
	tmp = deepcopy(model)
	remove_gene!(tmp, gene_names)
	return tmp
end

function remove_gene!{T <: String}(model::Model, gene_names::Array{T})
	gene_names = unique(gene_names)
	if any(x -> !in(x,model.genes), gene_names)
		warn("One or more gene does not appear in the model")
		return 
	end 
	
	genes_old = convert(Array{String}, gene_names)
	geneknockout = knockout_genes(model,genes_old)
	killed = geneknockout.disabled
	affected = geneknockout.affected

	indices = unique(sort(vcat(collect(values(affected))...)))
	safe_genes, safe_rules, new_old, _ = abstract_gene_names(model)
	old_new = invert_dict(new_old)
	genes_new = map(x -> old_new[x], genes_old)

	safe_rules = safe_rules[indices]

	for (i,v) in enumerate(safe_rules), gene in genes_new
	    safe_rules[i] = replace(safe_rules[i], gene, " false ")
	end

	for (i,v) in enumerate(safe_rules), gene in genes_in_reaction(v)
	    safe_rules[i] = replace(safe_rules[i], gene, new_old[gene])
	end
	model.rxn_rules[indices] = safe_rules
	model.lb[killed] = 0.0
	model.ub[killed] = 0.0
end

#------------------------------------------
#------------------------------------------


"""
    remove_reaction(model, reaction)
\n
    remove_reaction(model, reaction_index)

Return a duplicate of the model with a reaction removed
"""
function remove_reaction(model::Model, reaction::String)
    tmp = deepcopy(model)
    remove_reaction!(tmp, reaction)
    return tmp
end

"""
    remove_reactions!(model, reaction)
\n
    remove_reactions!(model, reaction_indices)
Remove a list of reactions from the model
"""
function remove_reaction!(model::Model, reaction::String)
    index = find_reaction(model, reaction)

    if length(index) == 0
        warn("Reaction not found in model")
        return 
    end 

    remove_reaction!(model,index)
end

function remove_reaction(model::Model, index::Number)
    tmp = deepcopy(model)
    remove_reaction!(tmp, index)
    return tmp
end 

function remove_reaction!(model::Model, index::Number)
    if !isinteger(index)
        warn("Reaction index must be integer")
        return
    end 

    if (index <= 0) | (index > length(model.rxns))
        warn("Invalid reaction index")
        return
    end 

    index = convert(Int64, index)
    remove_reaction!(model, [index])
end 

function remove_reaction{T <: String}(model::Model, reactions::Array{T})
    tmp = deepcopy(model)
	remove_reaction!(tmp,reactions)
    return tmp
end

function remove_reaction!{T <: String}(model::Model, reactions::Array{T})
	index = find_reaction(model, reactions)

    if any(x -> isempty(x), index)
        warn("One or more reactions not found in the model")
        return
    end 

	index = sort(unique(vcat(index...)))
	remove_reaction!(model, index)
end 

function remove_reaction{T <: Number}(model::Model, indices::Array{T})
    tmp = deepcopy(model)
	remove_reaction!(tmp, indices)
    return tmp
end

function remove_reaction!{T <: Number}(model::Model, indices::Array{T})
    indices = sort(unique(vcat(indices...)))
	deleteat!(model.rxns,indices)
    deleteat!(model.lb,indices)
    deleteat!(model.ub,indices)
    deleteat!(model.c,indices)
	deleteat!(model.rxn_rules,indices)
	deleteat!(model.rxn_name,indices)
	deleteat!(model.rxn_subsystem,indices)

    for (key,value) in model.rxn_extra 
        deleteat!(model.rxn_extra[key], indices)
    end 

	sort!(indices, rev=true)
	for i in indices
        model.S = delete_column(model.S, i)
		model.rxn_gene_mat = delete_column(model.rxn_gene_mat, i)
    end
end

#------------------------------------------------------
#------------------------------------------------------




