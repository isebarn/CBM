function find_blocked_reactions(model::Model; tolerance::Number = 1e-10)
	return find_blocked_reactions(model, tolerance)
end 

function find_blocked_reactions(model::Model, tolerance::Number = 1e-10)

	minF,maxF = fva(model,0);
	index = Int64[];

	# Find all min/max fluxes that approach 0
	for i = 1:length(minF)
		if (abs(minF[i]) < tolerance) && (abs(maxF[i]) < tolerance)
			index = [index; model.rxns[i]];
		end
	end

	return index;
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function find_deadend_metabolites(model::Model; lbfilter::Bool = true)
	return find_deadend_metabolites(model, lbfilter)
end

function find_deadend_metabolites(model::Model, lbfilter::Bool = true)
	dead_end = zeros(model.S.m)
	reactions = get_sparse_columns(model.S)

	# determine pos/neg numbers
	negative = (x -> model.S.nzval[x] < 0.0)
	positive = (x -> model.S.nzval[x] > 0.0)

	# reactions with a POSITIVE lower bound
	lb_filter = (x -> model.lb[reactions[x]] >= 0.0)

	for i = 1:length(model.b)

		# Need to find for each met, the reactions it participates in
		# Then filter out every met thats either always a producer or always
		# a consumer
		# Optionally I filter out those that have a positive lower bound >= 0.0
		met_s_location = findin(model.S.rowval, i)

		if lbfilter 
			met_s_location = filter(lb_filter, met_s_location)
		end

		if length(met_s_location) < 2
			continue
		end
		
		# arrays of neg/pos occurrences
		consume = map(negative, met_s_location) 
		produce = map(positive, met_s_location)

		# any(a) returns true if any element of a is true
		if any(consume)
			dead_end[i] -= 1
		end
		if any(produce)
			dead_end[i] += 1
		end	
	end

	return find(dead_end)
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function find_essential_genes(model::Model, min_growth::Number)
    # find all single gene deletions 
    ess_genes = gene_deletion(model,1).g_f

    # find those single gene deletions that bring growth below "min_growth"
    for (key,val) in ess_genes
        if val > min_growth
            delete!(ess_genes,key)
        end 
    end 
    ess_genes = collect(keys(ess_genes))
    ess_genes = vcat(ess_genes...)
    return ess_genes
end 

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function find_exchange_reactions(model::Model; output::String = "index")
	return find_exchange_reactions(model, output)
end 

function find_exchange_reactions(model::Model, output::String = "index")
	# Need to find all reactions that only appear ONCE in the reaction indices
	# as these are those that only map to A SINGLE metabolite
	count_single_occurences = (x -> length(nzrange(model.S,x)) == 1)

	index = map(count_single_occurences, collect(1:length(model.c)))
	index = find(index)
	if contains(lowercase(output),"name")
		return model.rxns[index]
	end

	return index
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function find_reactions_from_gene(model::Model, gene::String; output::String = "index")
	return find_reactions_from_gene(model, gene, output)
end

function find_reactions_from_gene(model::Model, gene::Int64; output::String = "index")
	return find_reactions_from_gene(model, gene, output)
end  

function find_reactions_from_gene(model::Model, gene::String, output::String = "index")
	index = find_gene_effects(model, gene)

	if in(lowercase(output),["n", "name"])
		return model.rxns[index]
	end 

	return index;
end

function find_reactions_from_gene(model::Model, gene::Int64, output::String = "index")
	gene = model.genes[gene]
	return find_reactions_from_gene(model, gene, output)
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function find_reactions_from_metabolite(model::Model, metabolite::String; output::String = "index")
	return find_reactions_from_metabolite(model, metabolite, output)
end 

function find_reactions_from_metabolite(model::Model, metabolite::String, output::String = "index")
	index = []

	if metabolite_exists(model,metabolite)
		push!(index,findin(model.S.rowval, find_metabolite(model,metabolite)))
	end

	index = vcat(index...)
	
	reactions = get_sparse_columns(model.S)
	
	index = reactions[index]

	if lowercase(output) == "formula" || lowercase(output[1]) == 'f'
		for i in index
			print_reaction_formula(model,model.rxns[i])
		end

	elseif lowercase(output) == "name" || lowercase(output[1]) == 'n'
		return model.rxns[index]
		
	else
		return index
	end
end

function find_reactions_from_metabolite(model::Model, metabolite::Int64; output::String = "index")
	return find_reactions_from_metabolite(model, metabolite, output)
end

function find_reactions_from_metabolite(model::Model, metabolite::Int64, output::String = "index")
	metabolite = model.mets[metabolite]
	return find_reactions_from_metabolite(model, metabolite, output)
end


# -------------------------------------------------------------------
# -------------------------------------------------------------------


function gene_deletion(model::Model, n::Number)
    gene_list, rule_list, saved_names, _ = abstract_gene_names(model)


    lp = setup_lp(model)
    wt = fba(model).f

    # create a dictionary containing the reactions that are affected 
    # by each geme, so that geneX => [rxnX1,...rxnXn]
    rxn_affected_by_gene = find_gene_effects(rule_list, gene_list)
    rxn_affected_by_gene = Dict(zip(gene_list,rxn_affected_by_gene))

    # empty dicts, yet to be filled 
    r_g_dict = Dict()
    r_f_dict = Dict()
    truth_table = Dict()


    for i in 1:n 

        # every gene combination
        combos = combinations(gene_list,i)
        num_combos = length(combos)
        
        # 1: Finding knockout effects
        # every gene combination affects some combination 
        # of reactions, rxns_affected contains the indices of 
        # these reactions 
        rxns_affected = Array(Array{Any, 1}, num_combos)
        for (index,combo) in enumerate(combos)
            rxns = map(x -> rxn_affected_by_gene[x], combo)
            rxns = sort(unique(vcat(rxns...)))
            
            rxns_affected[index] = rxns 
        end 
        
        # 28% of looptime    
        # 2: Convert rules to true/false statements
        # convert rules to "safe mode", where every gene is 
        # simply named gene_xxx, for string matching safety
        converted_rules = Array(Array{Any,1}, num_combos)
        for (index,combo) in enumerate(combos)
            rules = rule_list[rxns_affected[index]]
            convert_rule_to_boolean(combo, rules)

            converted_rules[index] = rules
        end 
        conv_rule_list = deepcopy(converted_rules)
        converted_rules = vcat(converted_rules...)

        
        # 19% of looptime 
        # 3: Evaluating every true/false statement
        # les út allar convertaðar reglur og eval(parse)
        new_truth_table = Dict(zip(converted_rules, fill(true, length(converted_rules))))
        rules_to_test = setdiff(collect(keys(new_truth_table)), collect(keys(truth_table)))
        truth_table = merge(new_truth_table, truth_table)
        for key in rules_to_test
            try 
            truth_table[key] = eval(parse(key))
        catch 
            println(key)
        end 
        end
        
        # 31% of looptime
        # 4: Finding disabled reactions and connecting to gene knockouts
        # map reactions to gene combinations and insert into a dictionary
        dead_reactions = Array(Array{Any,1}, num_combos)
        for (index, combo) in enumerate(combos)
            rules = conv_rule_list[index]

            dead = find(map(x -> !truth_table[x], rules))
            dead_reactions[index] = rxns_affected[index][dead]
        end 
        
        # 5: not sure
        # create a dictionary that comtains the reactions that 
        # get disabled for every particular gene combination
        unique_rxns = unique(dead_reactions)
        combo_arrays = Array(Array{Any, 1}, length(unique_rxns))
        (y -> map(x -> y[x] = [], 1:length(y)))(combo_arrays)
        rxn_combo_dict = Dict(zip(unique_rxns, combo_arrays))
        for (index,combo) in enumerate(combos)
            push!(rxn_combo_dict[dead_reactions[index]], combo)
        end 

        # 10% of total
        # 6: all the solver stuff
        # calculate fba for every combination of disabled reactions
        rxn_flow_dict = analyze_list_of_blocked_reaction_lp(lp, unique_rxns)
        
        # 7: Unioning dictionaries
        # find the union of the rxn_combo_dict and former loops
        intersects = intersect(collect(keys(rxn_combo_dict)), collect(keys(r_g_dict)))
        for clash in intersects
            r_g_dict[clash] = vcat(r_g_dict[clash], rxn_combo_dict[clash])
        end 
        r_g_dict = merge(rxn_combo_dict, r_g_dict)
        
        # 8: unioning more dictionaries
        # find the union of rxns_flow from other loops and this loop
        intersects = intersect(collect(keys(rxn_flow_dict)), collect(r_f_dict))
        for clash in intersects
            r_f_dict[clash] = vcat(r_f_dict[clash], rxn_flow_dict[clash])
        end 
        r_f_dict = merge(rxn_flow_dict, r_f_dict)

    end 

        
    for (key,val) in r_g_dict, (i,v) in enumerate(val), (j,k) in enumerate(v)
        v[j] = saved_names[k]
    end

    c_f_dict = Dict()
    for (key,val) in r_g_dict, gene in val
        c_f_dict[convert(Array{UTF8String},gene)] = r_f_dict[key]
    end

    return GeneDeletion(r_g_dict, r_f_dict, c_f_dict)
end 

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function knockout_genes(model::Model, gene::String)
    return knockout_genes(model, [gene])
end 

function knockout_genes(model::Model, gene_index::Number)
    if (gene_index <= 0) | (gene_index > length(model.genes))
        warn("Gene index is invalid")
    end 

    if !isinteger(gene_index)
        warn("Gene index must be integer")
    end 

    return knockout_genes(model, [gene_index])
end 

function knockout_genes{T <: Number}(model::Model, gene_indices::Array{T})
    gene_indices = unique(gene_indices)
    if !isempty(find(x -> (x <= 0) | (x > length(model.genes)), gene_indices))
        warn("One or more indices are invalid")
        return
    end   

    return knockout_genes(model, model.genes[gene_indices])  
end 

function knockout_genes{T <: String}(model::Model, genes::Array{T})
    if any(x -> !in(x,model.genes), genes)
        warn("One or more gene does not appear in the model")
        return 
    end 

    # find wild type growth
    wt = fba(model).f

    lp = setup_lp(model)
    
    # create a safe collection of rules
    gene_list, rule_list, saved_names, _ = abstract_gene_names(model)
    inverse_saved_names = invert_dict(saved_names)
    
    safe_genes = deepcopy(genes)
    convert_from_abstract(inverse_saved_names, safe_genes)
    rxns_affected = find_gene_effects(rule_list, safe_genes)
    gene_effect = Dict(zip(safe_genes, rxns_affected))
    rxns_affected = sort(unique(vcat(rxns_affected...)))

    rules_to_check = rule_list[rxns_affected]
    genes_to_kill = map(x -> inverse_saved_names[x], genes)
    convert_rule_to_boolean(genes_to_kill, rules_to_check)

    rxns_to_kill = []
    for (index, value) in enumerate(rules_to_check)
        push!(rxns_to_kill, !eval(parse(value)))
    end 


    killed_reactions = rxns_affected[find(rxns_to_kill)]

    map(x -> set_col_bounds(lp, x, 0.0, 0.0), killed_reactions)
    solve(lp)

    return GeneKnockout(get_objective_value(lp), get_solution(lp), killed_reactions, gene_effect)
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function knockout_reactions(model::Model, reaction::String)
    rxn_idx = find(model.rxns .== reaction)[1]

    if isempty(rxn_idx)
        warn("Reaction " * reaction * " not found")
        return 
    end 

    return knockout_reactions(model, [rxn_idx])
end

function knockout_reactions{T <: String}(model::Model, reactions::Array{T})
    rxn_idx = find(model.rxns .== reactions)

    if length(rxn_idx) != length(reactions)
        warn("Not all reactions found in model")
    end 

    return knockout_reactions(model, rxn_idx) 
end 

function knockout_reactions(model::Model, reactions::Number)
    return knockout_reactions(model, [reactions])
end  

function knockout_reactions{T <: Number}(model::Model, reactions::Array{T})
    gene_list, rule_list, saved_names, _ = abstract_gene_names(model)
    gene_list = Dict(zip(gene_list, collect(1:length(gene_list))))

    # get the gene reaction rules for the reactions I want to disable
    knockout_rules = rule_list[reactions]

    # check if any of the rules is "", i.e empty i.e not "knockoutable"
    sanity = find(map(x -> !contains(x, "gene_"), knockout_rules))
    if length(sanity) > 0
        for index in sanity
            warn("Reaction " * string(reactions[index]), " cannot be knocked out")
        end 
        deleteat!(knockout_rules, sanity)
    end

    # find the genes that appear in the reactions, a1,a2,a3...b1,b2,b3....n1,n2... 
    # put that list of genes into an array of arrays, RG = [[a1,a2,a3...], [b1,b2,b3...],......[n1,n2...]]
    r_g_array = Array[]
    for (index, rule) in enumerate(knockout_rules)
        push!(r_g_array, genes_in_reaction(rule))
    end 

    # locate the smallest gene knockout every reaction needs 
    Combo_Array = Array[]
    for (index, rule) in enumerate(knockout_rules)#, n in 1:length(r_g_array[index])
        push!(Combo_Array, [])
        knockout_found = false 
        
        n = 0
        while (n <= length(r_g_array[index])) & (!knockout_found)
            n += 1
            for combo in combinations(r_g_array[index], n)
                splt = split(rule)
                for gene in combo, (i,v) in enumerate(splt)
                    splt[i] = replace(splt[i], gene, "false")
                end 

                for (i,v) in enumerate(splt)
                    if startswith(v,"gene_")
                        splt[i] = "true"
                    end 
                end 

                if !eval(parse(join(splt, " ")))
                    knockout_found = true 
                    push!(Combo_Array[index],combo)
                end 

            end 
        end 
    end

    Reaction_Matrices = Array[]
    for (i,v) in enumerate(Combo_Array)
        push!(Reaction_Matrices, zeros(Int64, length(v[1]), length(v)))
        for (j,k) in enumerate(v)
            for(o,p) in enumerate(k)
                Reaction_Matrices[i][o,j] = gene_list[p]
            end 
        end
    end
    return Reaction_Matrices
    # With help from Tasos Papastylianou
    # http://stackoverflow.com/a/39087691/4551168
    # Find the smalest combinations of columns between matrices 
    # where every matrix represents reactions and columns represent the 
    # smallest set of genes that can disable the reaction 
    #########################################################################
    toJagged(x) = [x[:,i] for i in 1:size(x,2)];
    JaggedMatrices = [toJagged(x) for x in Reaction_Matrices];

    Combined = [unique(i) for i in JaggedMatrices[1]];
    for n in 2:length(JaggedMatrices)
        Combined = [unique([i;j]) for i in Combined, j in JaggedMatrices[n]];
    end

    # Find the minimum
    Lengths         = [length(s) for s in Combined];
    Minima          = findin(Lengths, min(Lengths...));
    SubscriptsArray = ind2sub(size(Lengths), Minima);
    #########################################################################

    genes = Combined[SubscriptsArray[1]]

    gene_list = invert_dict(gene_list)
    
    Genes = Array(AbstractString, length(genes[1]), length(genes))

    for (i,v) in enumerate(genes), (j,gene) in enumerate(v)
        Genes[j,i] = saved_names[gene_list[gene]]
    end 

    return Genes
end 



function print_medium(model)
	# List medium contents
	idx = find_exchange_reactions(model)
	for i in idx
	  if model.lb[i] < 0
	    @printf("%s\t%1.4f", model.rxns[i], model.lb[i])
	  end
	end
end

function print_reaction_formula(model::Model,reaction::String)
	id = find_reaction(model, reaction)
	id = vcat(id...)

	for i in id
		print_reaction_formula(model, i)
	end

end

function print_reaction_formula(model::Model, index::Number)

		sparse_indices = collect(nzrange(model.S,index))
		coefficients = model.S.nzval[sparse_indices]
		metabolites = model.S.rowval[sparse_indices]
		key = deepcopy(model.mets[index])

		right = ""
		left = ""

		for (i,v) in enumerate(coefficients)

			if v > 0
				right = string(right,v," ", model.mets[metabolites[i]], " + ")
			end
			if v < 0
				#left = string(left,metabolites[i]," + ")
				left = string(left,abs(v)," ", model.mets[metabolites[i]], " + ")
			end
		end

		right = right[1:end-2]
		left = left[1:end-2]
		println(left, " -> ", right)
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function robustness_analysis(model::Model, reactions::Array{String}; objective::String = "", pts = [], direction::String = "max")
    return robustness_analysis(model, reacitons, objective, pts, direction)
end 

function robustness_analysis(model::Model, rxns::Array{String}, objective::String = "", pts = [], direction::String = "max")
    rxns_idx = map(x -> find(x .== model.rxns), rxns)
    rxns_idx = filter(x -> length(x) == 1, rxns_idx)
    rxns_idx = vcat(rxns_idx...)

    if length(rxns_idx) != length(rxns)
        warn("One or more reactions not found in model.rxns")
        return 
    end 

    if objective == ""
        objective = model.rxns[find(model.c)[1]]
    end 

    obj_idx = find(model.rxns .== objective)

    if isempty(obj_idx)
        warn("Objective reaction not found in model.rxns")
    end 

    obj_idx = obj_idx[1]

    if pts == []
        pts = fill(20, length(rxns_idx))
    end 

    if length(pts) != length(rxns)
        warn("Please supply as many points as you supply reactions")
        return 
    end 

    return robustness_analysis(model, rxns_idx, obj_idx, pts, direction)
end 

function robustness_analysis{T <: Number}(model::Model, rxns::Array{T}; objective::Number = 0, pts = [], direction::String = "max")
    robustness_analysis(model, rxns; objective, pts, direction)
end 

function robustness_analysis{T <: Number}(model::Model, rxns::Array{T}, objective::Number = 0, pts = [], direction::String = "max")
	if objective == 0
		objective = find(model.c)[1]
	end
    if pts == []
        pts = fill(20, length(rxns))
    end 


    # just to relieve the code of endless length(model.c)
    num_rxns = length(rxns)

    # create a collection of indices for every single "fixation" of the reactions 
    # so that if I mean to fix reaction 5 in 4 different points and reaction 8 in 2 points 
    # combos will initially be [1,1], [1,2], [1,3], [1,4], ... [4,3], [4,4]
    # but idx will filter out those indices that I wont want, such as [1,4], point 4 is 
    # out of range for reaction 8, as I only want 1 and 2 for 8
    # so in the end, I will have 
    # [1,1],  [1,2],  [2,1], ... , [4,1], [4,2]
    combos = collect(combinations(repmat(1:maximum(pts), 1, num_rxns), num_rxns))
    combos = unique(combos)
    combos = sort(combos, lt=lexless)

    idx = map(y -> map(x -> combos[x][y] > pts[y], 1:length(combos)), 1:length(pts))
    idx = map(x -> find(x), idx)
    idx = union(idx...)
    idx = sort(idx, rev=true)
    map(x -> splice!(combos, x), idx)
    num_combos = length(combos)

    linspaces = map(x -> Array(Float64, x), pts)

    results = Array(Float64, num_combos)

    lp = setup_lp(model, direction) 

    for (i,v) in enumerate(rxns)
        minF, maxF = min_max_lp(lp, v)
        linspaces[i] = collect(linspace(minF, maxF, pts[i]))
    end

    change_objective_coef(lp, objective)

    last_combo = zeros(Int64, length(rxns))
    for (i, combo) in enumerate(combos) 
        combos_that_change = find(last_combo .!= combo)
        map((x,y) -> set_col_bounds(lp, rxns[x], linspaces[x][y]), combos_that_change, combo)

        results[i] = answer_lp(lp)
    end 

    # format the solution
    result_matrix = reshape(results, tuple(combos[end]...))

    return Robustness(result_matrix, rxns, linspaces)
end 

# -------------------------------------------------------------------
# -------------------------------------------------------------------



function synthetic_lethal_genes(model, cutoff = 10, nruns = 100)
    if cutoff < 1
        warn("Cutoff set as " * string(cutoff))
        warn("Cutoff is to be defined as a percentage not as a fraction of 100")
    end 

    # just a counter 
    cond_count = 0

    # calculate the wt_growth and gene_flow first.
    # see description of prep_genes() for information. The purpose is to 
    # format the data for quicker and easier lookup
    wt_growth = fba(model).f * (cutoff/100)
    gene_flow = gene_deletion(model,1).g_f
    rules, new_genes, gene_lookup, effects, rxn_genes = prep_genes(model)

    ngenes = length(model.genes)
    nrxns = length(model.rxns)

    # quick_gene is a dict with keys 1 => model.genes[1] ....
    # for faster lookup when replacing them in strings
    quick_gene = Dict(zip(1:ngenes, new_genes))
    quick_essential = Dict()
    quick_cond = Dict()
    
    # will be filled with evaluated booleans 
    # rxn_rules are converted to true/false statements 
    # that need to be evaluated, and there is tremendous benefit 
    # in only evaluating once instdead of repeadidly evaluating fx
    # the expression "true | false" 
    truth_table = Dict()

    # will be filled with shuffled up indices of genes
    perms = zeros(Int64, ngenes, nruns)

    # fill up the essential  array of essential genes
    # and sort them. In addition, create a dictionary quickgene 
    # for quick lookup
    for (key,val) in gene_flow
        if val < wt_growth
            quick_essential[find(key .== model.genes)[1]] = val
        end 
    end 

    # fill up the permutation matrix "perm" of gene indices and 
    # turn all "essentials" inside to 0 so that we dont 
    # repeatedly check genes that we know will kill
    for i = 1:nruns
        perms[:,i] = randperm(ngenes)
    end 
    perms = map(x -> haskey(quick_essential, x) ? 0 : x, perms)

    lp = setup_lp(model)
    lb, ub = get_variable_bounds(lp)
    for i in 1:nruns
        print(string("\n\u1b[1Frun ", i, " of ", nruns, " Number of conditionally essential genes found: ", cond_count))
        
        map(x -> set_col_bounds(lp, x, lb[x], ub[x]), 1:nrxns)


        tmp_rules = deepcopy(rules)

        # must do this every round, restard the lp and get its bounds 

        disabled = []

        for j in 1:ngenes 
            gene_idx = perms[j,i]
            if gene_idx == 0
                continue 
            end

            affected_idx = effects[gene_idx]
            affected_rules = deepcopy(tmp_rules[affected_idx])

            for (k,rule) in enumerate(affected_rules)
                affected_rules[k] = replace(affected_rules[k], quick_gene[gene_idx], " false ")
                affected_rules[k] = replace(affected_rules[k], r"x\[[0-9]*\]", " true ")
            end 


            for expression in affected_rules
                if !haskey(truth_table, expression)
                    truth_table[expression] = eval(parse(expression))
                end 
            end 

            # find the indices of the reactions that a gene 
            # disables, but only the indices of those reactions 
            # that haven yet been disabled
            dead = []
            for (k,r) in enumerate(affected_rules)
                if !truth_table[r] 
                    push!(dead, affected_idx[k])
                end 
            end 
            dead = setdiff(dead, disabled)

            map(x -> set_col_bounds(lp, x, 0.0, 0.0), dead)

            if answer_lp(lp) < wt_growth
                if !haskey(quick_cond, gene_idx)
                    quick_cond[gene_idx] = 0
                    cond_count += 1
                end 
                map(x -> set_col_bounds(lp, x, lb[x], ub[x]), dead)
            else 
                for (k,rule) in enumerate(tmp_rules)
                    tmp_rules[k] = replace(tmp_rules[k], quick_gene[gene_idx], " false ")
                end 
            end 
        end 
    end
    println()

    cond_ess = sort(collect(keys(quick_cond)))
    essential = sort(collect(keys(quick_essential)))

    non_essential = setdiff(collect(1:ngenes), cond_ess)
    non_essential = setdiff(non_essential, essential)
    
    return essential, cond_ess, non_essential
end 



