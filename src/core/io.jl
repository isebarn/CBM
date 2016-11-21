function load_json(file)
    try
        file = (endswith(file,".json")) ? JSON.parsefile(file) : JSON.parsefile(file*".json")
    catch y
        println(y)
        return
    end 

    metabolites = get(file, "metabolites", Dict());
    description = get(file, "id", file);
    reactions = get(file, "reactions", Dict());
    metabolites = get(file, "metabolites", Dict());
    genes = get(file, "genes", Dict());
    constraints = get(file, "constraints", "");

    met_main_fields = ["id", "name", "formula"];
    getmets(x) = map(y -> get(y, x, ""), metabolites)
    met_formula = getmets("formula")
    met_name = getmets("name")
    mets = getmets("id")

    rxn_main_fields = ["name", "id", "rxn_subsystem", "gene_reaction_rule"];
    getrxns(x) = map(y -> get(y, x, ""), reactions)
    rxns = getrxns("id")
    rxn_name = getrxns("name")
    rxn_subsystem = getrxns("subsystem")
    rxn_rules = getrxns("gene_reaction_rule")

    gene_main_fields = ["name", "id"];
    getgenes(x) = map(y -> get(y, x, ""), genes)
    gene = getgenes("id")
    gene_name = getgenes("name")

    constraint_main_fields = ["b", "csense"];
    if constraints != ""
        b = constraints["b"]
        csense = constraints["csense"]
    else 
        b = zeros(Float64, length(mets))
        csense = Dict(zip(["=", "<=", ">="], [collect(1:length(b)), [], []]))
    end 

    stochiometric_fields = ["lower_bound", "upper_bound", "objective_coefficient", "metabolites"]
    lb = getrxns("lower_bound")
    ub = getrxns("upper_bound")
    c = map(x -> get(x, "objective_coefficient", 0.0), reactions)

    ###########
    # Creating the Stochiometric Sparse Matrix 
    # for faster lookup, instead of searching an array for a met, search a dict
    tmpdict = Dict(zip(mets, collect(1:length(metabolites))))
    S = zeros(length(metabolites), length(reactions))

    for (i,reaction) in enumerate(reactions)
        for (key, val) in reaction["metabolites"]
            S[tmpdict[key], i] = val
        end 
    end 

    S = sparse(S)
    ###########

    ###########
    # Creating the Gene Reaction Matrix
    # for faster lookup, instead of searching an array for a gene, search a dict
    tmpdict = Dict(zip(gene, collect(1:length(gene))))
    
    # need to fix the rules for parsing efficiency
    rxn_rules = map(x -> replace(rxn_rules[x], " ( ", " ( "), 1:length(rxn_rules));
    rxn_rules = map(x -> replace(rxn_rules[x], " ) ", " ) "), 1:length(rxn_rules));
    rxn_rules = map(x -> replace(rxn_rules[x], " or ", " || "), 1:length(rxn_rules));
    rxn_rules = map(x -> replace(rxn_rules[x], " and ", " && "), 1:length(rxn_rules));

    # create tmprules to get rid of everything but gene ids
    tmprules = deepcopy(rxn_rules);
    tmprules = map(x -> replace(tmprules[x], "(", ""), 1:length(tmprules));
    tmprules = map(x -> replace(tmprules[x], ")", ""), 1:length(tmprules));
    tmprules = map(x -> replace(tmprules[x], "|", ""), 1:length(tmprules));
    tmprules = map(x -> replace(tmprules[x], "&", ""), 1:length(tmprules));
    tmprules = split.(tmprules)

    rxn_gene_mat = zeros(length(reactions), length(gene))

    for (i,v) in enumerate(tmprules)
        for k in v
            rxn_gene_mat[i, tmpdict[k]] = 1.0
        end 
    end 

    rxn_gene_mat = sparse(rxn_gene_mat)
    ######

    ############### All extra fields 
    met_extra_keys = map(x -> collect(keys(x)), metabolites);
    met_extra_keys = map(x -> setdiff(x, met_main_fields), met_extra_keys);
    met_extra_keys = unique(vcat(met_extra_keys...));
    met_extra_values = map(x -> getmets(x), met_extra_keys);
    met_extra = Dict(zip(met_extra_keys, met_extra_values));

    rxn_extra_keys = map(x -> collect(keys(x)), reactions);
    rxn_extra_keys = map(x -> setdiff(x, rxn_main_fields), rxn_extra_keys);
    rxn_extra_keys = map(x -> setdiff(x, stochiometric_fields), rxn_extra_keys);
    rxn_extra_keys = unique(vcat(rxn_extra_keys...));
    rxn_extra_values = map(x -> getrxns(x), rxn_extra_keys);
    rxn_extra = Dict(zip(rxn_extra_keys, rxn_extra_values));

    gene_extra_keys = map(x -> collect(keys(x)), genes);
    gene_extra_keys = map(x -> setdiff(x, gene_main_fields), gene_extra_keys);
    gene_extra_keys = unique(vcat(gene_extra_keys...));
    gene_extra_values = map(x -> getgenes(x), gene_extra_keys);
    gene_extra = Dict(zip(gene_extra_keys, gene_extra_values));
    ###################

    floatconvert(x) = convert(Array{Float64}, x)
    stringconvert(x) = convert(Array{String}, x)


    lb = floatconvert(lb)
    ub = floatconvert(ub)
    c = floatconvert(c)
    b = floatconvert(b)

    rxns = stringconvert(rxns)
    mets = stringconvert(mets)
    gene = stringconvert(gene)
    rxn_name = stringconvert(rxn_name)
    rxn_rules = stringconvert(rxn_rules)
    rxn_subsystem = stringconvert(rxn_subsystem)
    met_formula = stringconvert(met_formula)
    met_name = stringconvert(met_name)
    gene_name = stringconvert(gene_name)

    return Model(rxns,mets,gene,S,lb,ub,c,b,csense,rxn_gene_mat,rxn_name,rxn_rules,rxn_subsystem,rxn_extra,met_formula,met_name,met_extra,gene_name,gene_extra,description)
    ##    rxns
    ##    mets
    ##    genes
    ##    S
    ##    lb
    ##    ub
    ##    c
    ##    b
    ##    csense
    ##    rxn_gene_mat
    ##    rxn_name
    ##    rxn_rules
    ##    rxn_subsystem
    ##    rxn_extra
    ##    met_formula
    ##    met_name
    ##    met_extra
    ##    gene_name
    ##    gene_extra
    ##    description
end 

function load_matlab(filename; model_key = "", x...)
    mat_file = matread(filename)
    # extract filename from full path, ie a/b/model.mat => model.mat
    model_name = filename[search(filename, r"\w*.mat")]

    # This is necessary because there exist inconsistend model.mat files
    if haskey(mat_file, "model")
        model_dict = mat_file["model"]
    elseif haskey(mat_file, model_name)
        model_dict = mat_file[model_name]
    elseif length(collect(keys(mat_file))) == 1
        model_dict = mat_file[collect(keys(mat_file))[1]]
    elseif model_key != ""
        model_dict = mat_file[model_key]
    else 
        warn("Cannot determine which variable in " * model_name * ".mat is the model")
        warn("The keys stored in the file" * model_name * ".mat are:")
        counter = 0
        for (key,val) in mat_file 
            counter += 1
            warn(string(counter) * " : " * key)
        end 
        warn("Please call load_matlab(" * filename * ", X), where X is one of the keys")
        return 
    end 

    # need to check if the model_dict has every necessary key 
    # note csense and description are optional so I dont check if they are in the model_dict 
    model_keys = collect(keys(model_dict))
    extra_args = x
    extra_keys = map(y -> y[1], extra_args)
    extra_vals = map(y -> y[2], extra_args)
    

    fields = ["rxns","mets","genes","S","lb","ub","c","b","csense","rxnGeneMat","rxnNames","grRules","subSystems","metFormulas","metNames","description"]
    
    if !isempty(setdiff(fields, vcat(model_keys, ["csense", "description"])))
        missing = setdiff(fields, model_keys)
        warn("Not all necessary keys are present in the file")
        warn("Missing keys are: ")
        counter = 0
        for i in missing 
            counter += 1
            warn(string(counter) * " : " * i)
        end 

        println()

        warn("Unused keys in the matlab file are:")

        counter = 0

        for i in setdiff(vcat(model_keys, ["csense", "description"]), fields)
            counter += 1
            warn(string(counter) * " : " * i)
        end 

        println()

        warn("You can run load_matlab(filename, " * (model_key != "" ? model_key : " ") * "X = \"Y\")")
        warn("Where X is the missing field, and Y is key of that field in the matlab file")
    end 


    # create the model and start filling in 
    model = Model([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])

    
    # intersect the fields array with the model_keys array and insert those values
    model_fields = ["rxns","mets","genes","S","lb","ub","c","b","csense","rxn_gene_mat","rxn_name","rxn_rules","rxn_subsystem","met_formula","met_name","description"]
    for (i,v) in enumerate(model_fields)
        setfield!(model, Symbol(v), get(model_dict, fields[i], [])[1:end])
    end

    # use the varargs fields, if there are any supplied 
    for i in 1:length(extra_keys)
        if in(extra_keys[i], fields)
            idx = find(map(y -> extra_keys[i] == y, fields))[1]
            setfield!(model, Symbol(model_fields[i]), get(model_dict, extra_vals[i], [])[1:end])
        end 
    end 

    # get the rxn_extra, met_extra and gene_extra fields, if there are any in the model 
    extra_fields = setdiff(model_keys, fields)
    rxn_extra = Dict()
    met_extra = Dict()
    gene_extra = Dict()

    for i in extra_fields 
        extra = model_dict[i][1:end]

        if length(extra) == length(model.rxns)
            rxn_extra[i] = extra
        elseif length(extra) == length(model.mets)
            met_extra[i] = extra 
        elseif length(extra) == length(model.genes)
            gene_extra[i] = extra 
        end 
    end 

    model.rxn_extra = rxn_extra
    model.met_extra = met_extra
    model.gene_extra = gene_extra

    # check the csense field...tricky business
    if length(model.csense) == 0 
        model.csense = Dict(zip(["=", "<=", ">="], [collect(1:length(model.b)), [], []]))
    elseif (typeof(model.csense) <: String)
        if length(model.csense) == length(model.rxns)
            L = []
            E = []
            U = []
            for i in 1:length(model.csense)
                if model.csense[i] == 'L'
                    push!(L, i)
                elseif model.csense[i] == 'U'
                    push!(U, i)
                else 
                    push!(E, i)
                end 
            end 
            model.csense = Dict(zip(["<=", "=", ">="], [L, E, U]))
        else 
            model.csense = Dict(zip(["=", "<=", ">="], [collect(1:length(model.b)), [], []]))
        end 
    else 
        model.csense = Dict(zip(["=", "<=", ">="], [collect(1:length(model.b)), [], []]))
    end 


    # need to reshape this because of the [1:end] call everywhere
    # which flattens out the sparse matrix
    model.S = sparse(reshape(model.S, length(model.mets), length(model.rxns)))

    return model
end 


function load_matlab_variable(file)
	return matread(file)
end

function load_matlab_array(file, variable_name)
	array = load_matlab_variable(file)[variable_name]
	array = array'
	array = array[:]
	
	try
		array = convert(Array{Int64}, array)
	catch 
		array = convert(Array{Float64}, array)
	end 
end

function write_matlab_variable(destination, variable_name, variable)
	file = matopen(destination, "w")
	try 
		write(file, variable_name, variable)
	catch y
		warn("Error writing variable")
		warn(y)
	end 
	close(file)
end 


function convert_to_uft8!(array)
	return convert(Array{UTF8String}, array)
end
function load_table(filename, delimiter = ',')
    rxn_table = readdlm(filename * "_rxns.csv", delimiter);
    met_table = readdlm(filename * "_mets.csv", delimiter);
    gene_table = readdlm(filename * "_genes.csv", delimiter);
    constraint_table = readdlm(filename * "_constraints.csv", delimiter);

    # Sometimes all elements in a column are type SubString{String}, so I locate these columns 
    # and convert their elements to String
    #
    # Locate the columns
    rxn_substring_cols = find(map(x -> typeof(x) <: SubString, rxn_table[1,:]))
    met_substring_cols = find(map(x -> typeof(x) <: SubString, met_table[1,:]))
    gene_substring_cols = find(map(x -> typeof(x) <: SubString, gene_table[1,:]))
    constraint_table_substring_cols = find(map(x -> typeof(x) <: SubString, constraint_table[1,:]))
    #
    # Convert the columns elements 
    map(x -> map(y -> rxn_table[y,x] = convert(String, rxn_table[y,x]), 1:size(rxn_table)[1]), rxn_substring_cols);
    map(x -> map(y -> met_table[y,x] = convert(String, met_table[y,x]), 1:size(met_table)[1]), met_substring_cols);
    map(x -> map(y -> gene_table[y,x] = convert(String, gene_table[y,x]), 1:size(gene_table)[1]), gene_substring_cols);
    map(x -> map(y -> constraint_table[y,x] = convert(String, constraint_table[y,x]), 1:size(constraint_table)[1]), constraint_table_substring_cols);


    model = Model([],[],[],sparse(zeros(size(met_table)[1],0)),[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])

    # Metabolites 
    ##########################
    model.mets = met_table[:,1]
    model.met_formula = met_table[:,2]
    model.met_name = met_table[:,3]
    ##########################

    # Reactions 
    ##########################
    model.rxns = rxn_table[:,1]
    model.lb = rxn_table[:,2]
    model.ub = rxn_table[:,3]
    model.c = rxn_table[:,4]
    model.rxn_name = rxn_table[:,5]
    model.rxn_rules = rxn_table[:,6]
    model.rxn_subsystem = rxn_table[:,7]
    metdict = Dict(zip(model.mets, collect(1:length(model.mets))))
    for (i,v) in enumerate(rxn_table[:,8])
        model.S = addrxn(model.S, v, metdict)
    end 
    ##########################


    # Genes 
    ##########################
    model.genes = gene_table[:,1]
    model.gene_name = gene_table[:,2]
    ##########################


    # Constraints 
    ##########################
    model.b = constraint_table[:,1]

    model.csense = Dict()
    model.csense["<="] = find(x -> x == "<=", constraint_table[:,2])
    model.csense[">="] = find(x -> x == "=>", constraint_table[:,2])
    model.csense["="] = find(x -> x == "=", constraint_table[:,2])
    ##########################

    return model
end 

function addrxn(S, rxn, metdict)
    if contains(rxn, "<=>")
        middle = "<=>"
    elseif contains(rxn, "<=")
        middle = "<="
    elseif contains(rxn, "=>")
        middle = "=>"
    else 
        warn("Reaction must contain '<=>', '<=' or '=>'")
        return 
    end 

    S = add_column(S)
    col = S.n 

    sides = split(rxn, middle);
    left = sides[1];


    if length(left) != 0
        left_mets = split.(split(left, " + "))
        left_mets = map(y -> filter(x -> x != "", y), left_mets)
        left_mets = map(x -> convert(Array{Any}, x), left_mets)
        map(x -> length(x) == 1 ? unshift!(x, "1.0") : x, left_mets)
        left_mets = map(x -> [-parse(x[1]), metdict[x[2]]], left_mets)

        for i in left_mets 
            S[Int(i[2]), col] = i[1]
        end 
    end 

    if length(sides) == 1
        return S 
    end 

    right = sides[2];
    if length(right) != 0
        right_mets = split.(split(right, " + "))
        right_mets = map(y -> filter(x -> x != "", y), right_mets)
        right_mets = map(x -> convert(Array{Any}, x), right_mets)
        map(x -> length(x) == 1 ? unshift!(x, "1.0") : x, right_mets)
        right_mets = map(x -> [parse(x[1]), metdict[x[2]]], right_mets)



        for i in right_mets 
            S[Int(i[2]), col] = i[1]
        end 
    end 

    return S
end 
function export_json(model, filename)
    model_dict = Dict()

    ####### Metabolites
    model_dict["metabolites"] = []
    for (i,v) in enumerate(model.mets)
        metabolite = Dict()
        metabolite["id"] = model.mets[i]
        metabolite["formula"] = model.met_formula[i]
        metabolite["name"] = model.met_name[i]

        for (key,val) in model.met_extra
            metabolite[key] = val[i]
        end

        push!(model_dict["metabolites"], metabolite)
    end 

    ####### Description
    model_dict["id"] = model.description

    ####### Reaction
    model_dict["reactions"] = []
    for (i,v) in enumerate(model.rxns)
        reaction = Dict()
        reaction["id"] = model.rxns[i]
        reaction["name"] = model.rxn_name[i]
        reaction["subsystem"] = model.rxn_subsystem[i]
        reaction["gene_reaction_rule"] = model.rxn_rules[i]
        reaction["lower_bound"] = model.lb[i]
        reaction["upper_bound"] = model.ub[i]
        reaction["objective_coefficient"] = model.c[i]

        # get metabolites and coefficients from model.S
        idx = model.S.rowval[nzrange(model.S,i)] ;
        rxn_mets = model.mets[idx];
        met_coeffs = model.S.nzval[nzrange(model.S,i)];
        reaction["metabolites"] = Dict(zip(rxn_mets, met_coeffs));

        for (key,val) in model.rxn_extra 
            reaction[key] = val[i]
        end 

        push!(model_dict["reactions"], reaction)
    end 

    ####### Genes
    model_dict["genes"] = []
    for (i,v) in enumerate(model.genes)
        gene = Dict()
        gene["id"] = model.genes[i]
        gene["name"] = model.gene_name[i]

        for (key,val) in model.gene_extra 
            gene[key] = val[i]
        end

        push!(model_dict["genes"], gene)
    end 

    ####### Constraints
    model_dict["constraints"] = Dict(zip(["b", "csense"], [model.b, model.csense]))


    ####### Make JSON string and write it to file
    jsonstring = json(model_dict)

    file = open(filename, "w")
    write(file, jsonstring)
    close(file)
end 
function export_matlab(model, filename)
    matlab_model = Dict()

    matlab_model["mets"] = model.mets 
    matlab_model["metNames"] = model.met_name 
    matlab_model["metFormulas"] = model.met_formula 
    matlab_model["rxns"] = model.rxns
    matlab_model["rxnNames"] = model.rxn_name 
    matlab_model["subSystems"] = model.rxn_subsystem 
    matlab_model["lb"] = model.lb 
    matlab_model["ub"] = model.ub 
    matlab_model["c"] = model.c 
    matlab_model["b"] = model.b 
    matlab_model["S"] = model.S 
    matlab_model["rxnGeneMat"] = model.rxn_gene_mat 
    matlab_model["grRules"] = model.rxn_rules 
    matlab_model["genes"] = model.genes 
    matlab_model["gene_name"] = model.gene_name

    # csense conversion from dictionary to array of L, E and U 
    ###########################################################
    csense_length = length(model.csense["<="]) + length(model.csense["="]) + length(model.csense[">="]);
    csense = Array(Char, csense_length);
    csense[model.csense["<="]] = 'L';
    csense[model.csense["="]] = 'E';
    csense[model.csense[">="]] = 'U';
    csense = join(csense)
    matlab_model["csense"] = csense
    ###########################################################

    for (key,val) in model.rxn_extra 
        matlab_model[key] = val 
    end 

    for (key,val) in model.met_extra 
        matlab_model[key] = val 
    end 

    for (key,val) in model.gene_extra 
        matlab_model[key] = val 
    end 

    exported_model = Dict(zip(["model"], [matlab_model]))

    MAT.matwrite(filename, exported_model)
end 

function gene_rule_generator(model)
	tmp = deepcopy(model.rules)

	for (i,v) in enumerate(tmp) 
		genes = genes_in_reaction(v)

		for gene in genes 
			g = string("x(", find_gene(model, gene), ")")
			tmp[i] = replace(tmp[i], gene, g)
		end 
	end 
	return tmp
end 

function sparse_gene_generator(model)
	sparse = spzeros(length(model.c), length(model.genes))

	for (i,v) in enumerate(model.genes)
		rxns = find_reactions_from_gene(model, v)

		for r in rxns 
			sparse[r,i] = 1
		end 
	end 

	return sparse
end 

function export_table(model, filename; delimiter = ',', overwrite = false )
    export_table(model, filename, delimiter, overwrite)
end 

function export_table(model, filename, delimiter = ',', overwrite = false )
    println(filename)
    if !overwrite 
        exists = false 
        if isfile(filename * "_rxns.csv")
            warn("File " * filename * "_rxns.csv exists in this location")
            exists = true 
        end
        if isfile(filename * "_mets.csv")
            warn("File " * filename * "_mets.csv exists in this location")
            exists = true 
        end
        if isfile(filename * "_genes.csv")
            warn("File " * filename * "_genes.csv exists in this location")
            exists = true 
        end
        if isfile(filename * "_constraints.csv")
            warn("File " * filename * "_constraints.csv exists in this location")
            exists = true 
        end

        if exists 
            warn("To overwrite call exp_table(model, " * filename * ", overwrite = true)")
            return 
        end 
    end 


    ################# Reactions 
    ##########################################
    rxn_data = deepcopy(model.rxns);
    rxn_data = hcat(rxn_data, model.lb);
    rxn_data = hcat(rxn_data, model.ub);
    rxn_data = hcat(rxn_data, model.c);
    rxn_data = hcat(rxn_data, model.rxn_name);
    rxn_data = hcat(rxn_data, model.rxn_rules);
    rxn_data = hcat(rxn_data, model.rxn_subsystem);

    ##### Create equations from model.S
    ###########################################################
    equations = [] 
    for i in 1:length(model.c)
        rxn_mets = model.mets[rowvals(model.S)[nzrange(model.S,i)]]
        coeffs = model.S.nzval[nzrange(model.S,i)]

        left = find(x -> x < 0, coeffs)
        if !isempty(left)
            leftside = map(x -> (coeffs[x] == -1 ? "" : string(abs(coeffs[x]))), left)
            leftside .*= " " 
            leftside .*= rxn_mets[left]
            leftside = join(leftside, " + ")
        else 
            leftside = ""
        end 

        right = find(x -> x > 0, coeffs)
        if !isempty(right)
            rightside = map(x -> (coeffs[x] == 1 ? "" : string(coeffs[x])), right)
            rightside .*= " " 
            rightside .*= rxn_mets[right]
            rightside = join(rightside, " + ")
        else 
            rightside = ""
        end 

        bounds = sign(model.lb[i]) + sign(model.ub[i])

        if bounds == 0.0 
            direction = " <=> "
        elseif bounds == -1.0
            direction = " <= "
        else 
            direction = " => "
        end 

        equation = leftside * direction * rightside
        equation = rstrip(equation)
        equation = lstrip(equation)

        push!(equations, equation)
    end 
    rxn_data = hcat(rxn_data, equations);
    ###########################################################




    ################### Metabolites
    ###########################################

    met_data = deepcopy(model.mets);
    met_data = hcat(met_data, model.met_formula);
    met_data = hcat(met_data, model.met_name);




    ################### Genes
    ###########################################

    gene_data = deepcopy(model.genes);
    gene_data = hcat(gene_data, model.gene_name);




    ################### Constraints
    ###########################################

    constraint_data = deepcopy(model.b)

    # csense conversion from dictionary to array of L, E and U 
    ###########################################################
    csense_length = length(model.csense["<="]) + length(model.csense["="]) + length(model.csense[">="]);
    csense = Array(String, csense_length);
    csense[model.csense["<="]] = "<=";
    csense[model.csense["="]] = "=";
    csense[model.csense[">="]] = "=>";
    constraint_data = hcat(constraint_data, csense)
    ###########################################################
    
    writedlm(filename * "_rxns.csv",rxn_data , delimiter)
    writedlm(filename * "_mets.csv",met_data , delimiter)
    writedlm(filename * "_genes.csv",gene_data , delimiter)
    writedlm(filename * "_constraints.csv",constraint_data , delimiter)
    #write_delimiter_file(rxn_data, "_rxns", filename, delimiter)
    #write_delimiter_file(met_data, "_mets", filename, delimiter)
    #write_delimiter_file(gene_data, "_genes", filename, delimiter)
    #write_delimiter_file(constraint_data, "_constraints", filename, delimiter)
end 
