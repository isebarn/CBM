function test_add_reaction(model)
    value1 = 1.31592;
	value2 = 21.2434;
    tmp = deepcopy(model)

    # testing the first one
    test_model_1 = add_reaction(tmp,"test1", "pyr_e <=> pyr_c");
    @test_approx_eq value1 round(fba(test_model_1).obj,5)
    @test test_model_1 != tmp
    @test test_model_1 != model 

    @test_approx_eq length(tmp.c) length(model.c) 
    @test_approx_eq length(tmp.lb) length(model.lb) 
    @test_approx_eq length(tmp.ub) length(model.ub) 
    @test_approx_eq length(tmp.rxns) length(model.rxns) 
    @test_approx_eq length(tmp.mets) length(model.mets) 

    @test_approx_eq length(test_model_1.c) length(model.c) + 1
    @test_approx_eq length(test_model_1.lb) length(model.lb) + 1
    @test_approx_eq length(test_model_1.ub) length(model.ub) + 1
    @test_approx_eq length(test_model_1.rxns) length(model.rxns) + 1
    @test_approx_eq length(test_model_1.mets) length(model.mets) 
    @test_approx_eq test_model_1.S.m tmp.S.m 
    @test test_model_1.S != tmp.S
    @test_approx_eq test_model_1.S.n tmp.S.n + 1
    test_model_2 = add_reaction(tmp, "test1", "pyr_e <=> test_c");
    test_model_2 = add_reaction!(test_model_2, "test2", "test_e <=> test_c");
    test_model_2 = add_reaction!(test_model_2, "test3", "test_e <=> ");
    @test_approx_eq_eps(value2, fba(test_model_2).obj, 1e-4)
    @test test_model_2 != tmp

    @test_approx_eq length(tmp.c) length(model.c) 
    @test_approx_eq length(tmp.lb) length(model.lb) 
    @test_approx_eq length(tmp.ub) length(model.ub) 
    @test_approx_eq length(tmp.rxns) length(model.rxns) 
    @test_approx_eq length(tmp.mets) length(model.mets) 
    @test test_model_2 != model 
    @test_approx_eq length(test_model_2.c) length(model.c) + 3
    @test_approx_eq length(test_model_2.lb) length(model.lb) + 3
    @test_approx_eq length(test_model_2.ub) length(model.ub) + 3
    @test_approx_eq length(test_model_2.rxns) length(model.rxns) + 3
    @test_approx_eq length(test_model_2.mets) length(model.mets) + 2
    @test_approx_eq test_model_2.S.m tmp.S.m + 2
    @test test_model_2.S != tmp.S 
    @test_approx_eq test_model_2.S.n tmp.S.n + 3
end

function test_change_objective(model)
	value = 20.0
	safe_model = deepcopy(model)

	safe_model.c = zeros(Float64, length(safe_model.c))
	safe_model.c[13] = 1
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, 13, 1.0)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, 13; objective = 1.0)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, 13; objective = 1)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, model.rxns[13]; objective = 1)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, model.rxns[13], 1.0)).obj
	@test_approx_eq safe_model.c model.c

	safe_model.c = zeros(Float64, length(safe_model.c))
	safe_model.c[18] = 1
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, 18, 1.0)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, 18; objective = 1.0)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, 18; objective = 1)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, model.rxns[18]; objective = 1)).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, model.rxns[18], 1.0)).obj
	@test_approx_eq safe_model.c change_objective(model, model.rxns[18], 1).c
	@test_approx_eq sum(safe_model.c) sum(change_objective(model, model.rxns[18], 1).c)

	safe_model.c = zeros(Float64, length(safe_model.c))
	safe_model.c[5] = 1
	safe_model.c[18] = 1
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, [5,18])).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, [5,18],[1,1])).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, [5,18],[1.0,1.0])).obj
	@test_approx_eq sum(safe_model.c) sum(change_objective(model, [5,18]).c)
	@test_approx_eq sum(safe_model.c) sum(change_objective(model, [5,18], [1,1]).c)
	@test_approx_eq sum(safe_model.c) sum(change_objective(model, [5,18], [1.0,1.0]).c)


	safe_model.c = zeros(Float64, length(safe_model.c))
	safe_model.c[5] = 0.5
	safe_model.c[18] = 1
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, [5,18], [0.5, 1.0])).obj
	@test_approx_eq fba(safe_model).obj fba(change_objective(model, [5.0,18.0]; objectives = [0.5, 1.0])).obj
	@test_approx_eq sum(safe_model.c) sum(change_objective(model, [5,18], [0.5, 1.0]).c)
	@test_approx_eq sum(safe_model.c) sum(change_objective(model, [5,18], [0.5, 1.0]).c)	
end


function test_change_reaction_bounds(model)
	play_model = deepcopy(model)
	safe_model = deepcopy(model)
	safe_model.lb[13] = 0.5
	safe_model.ub[13] = 0.5
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, model.rxns[13], 0.5, "both")).obj
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, "both").lb[13] 0.5
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, "both").ub[13] 0.5
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, "l").ub[13] play_model.ub[13]
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, "l").lb[13] 0.5
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, "u").ub[13] 0.5
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, "u").lb[13] play_model.lb[13]
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, model.rxns[13], 0.5, "u")).obj
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, model.rxns[13], 0.5)).obj
	safe_model.ub[13] = 0.9
	@test abs(fba(safe_model).obj - fba(change_reaction_bounds(model, model.rxns[13], 0.5, "l")).obj) < 1e-6

	safe_model = deepcopy(model)
	safe_model.lb[13] = 0.5
	safe_model.ub[13] = 0.7
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, model.rxns[13], 0.5, 0.7,  "both")).obj
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, model.rxns[13], 0.5, 0.7)).obj
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, model.rxns[13], 0.5, 0.7, "u")).obj
	safe_model.ub[13] = 0.9
	@test abs(fba(safe_model).obj - fba(change_reaction_bounds(model, model.rxns[13], 0.5, 0.9, "l")).obj) < 1e-6
	
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, 0.7,  "both").lb[13] 0.5
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, 0.7,  "both").ub[13] 0.7
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, 0.7,  "l").lb[13] 0.5
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, 0.7,  "l").ub[13] play_model.ub[13]
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, 0.7,  "u").lb[13] play_model.lb[13]
	@test_approx_eq change_reaction_bounds(model, model.rxns[13], 0.5, 0.7,  "u").ub[13] 0.7

	safe_model = deepcopy(model)
	safe_model.lb[13] = 0.5
	safe_model.ub[13] = 0.7
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, 13, 0.7,  "both")).obj
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, 13, 0.7,  "u")).obj
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, 13, 0.7)).obj
	safe_model.lb[13] = 0.7
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, 13, 0.7, "both")).obj
	safe_model.ub[13] = 0.9
	@test abs(fba(safe_model).obj - fba(change_reaction_bounds(model, 13, 0.9, "u")).obj) < 1e-6
	@test abs(fba(safe_model).obj - fba(change_reaction_bounds(model, 13, 0.1, "l")).obj) < 1e-6
	@test_approx_eq change_reaction_bounds(model, 13, 0.7,  "both").lb[13] 0.7
	@test_approx_eq change_reaction_bounds(model, 13, 0.7,  "both").ub[13] 0.7
	@test_approx_eq change_reaction_bounds(model, 13, 0.7,  "l").lb[13] 0.7
	@test_approx_eq change_reaction_bounds(model, 13, 0.7,  "l").ub[13] play_model.ub[13]
	@test_approx_eq change_reaction_bounds(model, 13, 0.7,  "u").lb[13] play_model.lb[13]
	@test_approx_eq change_reaction_bounds(model, 13, 0.7,  "u").ub[13] 0.7

	safe_model = deepcopy(model)
	safe_model.lb[13] = 0.5
	safe_model.ub[13] = 0.7
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, 13, 0.5, 0.7,  "both")).obj
	@test_approx_eq fba(safe_model).obj fba(change_reaction_bounds(model, 13, 0.5, 0.7,  "u")).obj
	safe_model.ub[13] = 0.9
	@test abs(fba(safe_model).obj - fba(change_reaction_bounds(model, 13, 0.5, 0.7,  "l")).obj) < 1e-6
	@test_approx_eq change_reaction_bounds(model, 13, 0.5, 0.7,  "both").lb[13] 0.5
	@test_approx_eq change_reaction_bounds(model, 13, 0.5, 0.7,  "both").ub[13] 0.7
	@test_approx_eq change_reaction_bounds(model, 13, 0.5, 0.7,  "l").lb[13] 0.5
	@test_approx_eq change_reaction_bounds(model, 13, 0.5, 0.7,  "l").ub[13] play_model.ub[13]
	@test_approx_eq change_reaction_bounds(model, 13, 0.5, 0.7,  "u").lb[13] play_model.lb[13]
	@test_approx_eq change_reaction_bounds(model, 13, 0.5, 0.7,  "u").ub[13] 0.7
end
function test_fba(model)
	value = 0.87392

    fbaA = fba(model).obj
	@test_approx_eq_eps value fbaA 1e-5
	#@test value == round(fba(model,"max","gurobi").obj,5)
	@test_approx_eq_eps value fba(model,"max").obj 1e-5
end
function test_find_blocked_reactions(model::Model)
	index = [26,27,29,34,45,47,52,63];
	values = find_blocked_reactions(model);
	@test values == model.rxns[index];
end
function test_find_deadend_metabolites(model)
	ans = [30;32;37;49;60]

	@test ans == find_deadend_metabolites(model)
end
function test_find_exchange_reactions(model::Model)
	index = [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39];
	values = find_exchange_reactions(model);

	@test index == values
end

function test_find_reactions_from_gene(model::Model)
	index = [2,14,58,69,70]
	values = find_reactions_from_gene(model, "s0001")

	@test index == values
end
function test_find_reactions_from_metabolite(model::Model)
	index = [4,5,11,12,13,15,16,18,41,46,52,53,54,58,62,76,79,81]
	value = find_reactions_from_metabolite(model, "h2o_c")

	@test index == value
end
function test_fva(model::Model)

    # this is a necessary fix 
    # the ecoli_core.mat model 
    # which the /CBM/test/test_vars/fva_test.mat
    # values test values are calculated
    # from has these two values swapped 
    # relative to the e_coli_core.json model 
    # so we must do that aswell
    # in addition, the ecoli_core.mat has 
    # the upper limit for ATPM at 1000 while 
    # e_coli_core.json has it at 8.39
    tmp = deepcopy(model)
    tmp.S.rowval[111] = 42
    tmp.S.rowval[112] = 44
    tmp.ub[11] = 1000

    fva_test = open_mat_file(Pkg.dir() * "/CBM/test/test_vars/fva_test.mat")

    minF0 = fva_test["minF0"]
    maxF0 = fva_test["maxF0"]
    
    minF90 = fva_test["minF90"]
    maxF90 = fva_test["maxF90"]
    
    minF100 = fva_test["minF100"]
    maxF100 = fva_test["maxF100"]

    minF, maxF = fva(tmp, 0)
    @test_approx_eq minF minF0
    @test_approx_eq maxF maxF0

    minF, maxF = fva(tmp, 0.9)
    @test_approx_eq_eps minF minF90 1e-3
    @test_approx_eq_eps maxF maxF90 1e-3

    minF, maxF = fva(tmp, 1)
    @test_approx_eq_eps minF minF100 1e-3
    @test_approx_eq_eps maxF maxF100 1e-3
end 
function test_gene_deletion(model)
    single_gene = open_mat_file(Pkg.dir() * "/CBM/test/test_vars/sngl_gene.mat")

    matsol = single_gene["sol"]
    genes = single_gene["genelist"]


    sol = gene_deletion(model,1).g_f 

    @test sum(matsol) - sum(collect(values(sol))) < 1e-6

    solution_order = [] 

    for gene in genes 
        push!(solution_order, sol[[gene]])
    end 

    @test !any(abs(matsol - solution_order) .> 1e-3)

    double_matrix = zeros(Float64, length(genes), length(genes))


    sol = gene_deletion(model,2).g_f

    for i in 1:length(genes)
        for j in i+1:length(genes)
            if haskey(sol, [genes[i], genes[j]])
                double_matrix[i,j] = sol[[genes[i], genes[j]]]
            elseif haskey(sol, [genes[j], genes[i]])
                double_matrix[i,j] = sol[[genes[j], genes[i]]]
            else 
                warn("something wrong")
                break 
            end 
        end 
    end 

    double_matrix += double_matrix' 

    for i in 1:length(genes)
        double_matrix[i,i] = sol[[genes[i]]]
    end 

    double_gene = open_mat_file(Pkg.dir() * "/CBM/test/test_vars/dbl_gene.mat")["ss"]

    @test !any(abs(double_gene - double_matrix) .> 1e-3)
end 
function test_remove_reaction(model)
    eps = 1e-2

    test = deepcopy(model)
    tmp = deepcopy(model)

    tmp.lb[12] = 0
    tmp.ub[12] = 0
    @test_approx_eq_eps(fba(tmp).obj,fba(remove_reaction(test, model.rxns[12])).obj, eps)
    remove_reaction!(test, model.rxns[12])
    @test_approx_eq_eps(fba(tmp).obj, fba(test).obj, eps)
    @test_approx_eq length(model.rxns) length(test.rxns) + 1
    @test_approx_eq length(model.lb) length(test.lb) + 1
    @test_approx_eq length(model.ub) length(test.ub) + 1
    @test_approx_eq length(model.c) length(test.c) + 1
    @test_approx_eq length(model.rxn_rules) length(test.rxn_rules) + 1
    @test_approx_eq length(model.rxn_name) length(test.rxn_name) + 1
    @test_approx_eq length(model.rxn_subsystem) length(test.rxn_subsystem) + 1    


    tmp.lb[8] = 0
    tmp.ub[8] = 0

    @test_approx_eq_eps(fba(tmp).obj,fba(remove_reaction(test, model.rxns[8])).obj, eps)

    remove_reaction!(test, model.rxns[8])
    @test_approx_eq length(model.rxns) length(test.rxns) + 2
    @test_approx_eq length(model.lb) length(test.lb) + 2
    @test_approx_eq length(model.ub) length(test.ub) + 2
    @test_approx_eq length(model.c) length(test.c) + 2
    @test_approx_eq length(model.rxn_rules) length(test.rxn_rules) + 2
    @test_approx_eq length(model.rxn_name) length(test.rxn_name) + 2
    @test_approx_eq length(model.rxn_subsystem) length(test.rxn_subsystem) + 2


    test = deepcopy(model)

    tmp = deepcopy(model)  
    tmp.lb[12] = 0
    tmp.ub[12] = 0

    @test_approx_eq_eps(fba(tmp).obj,fba(remove_reaction(test, 12)).obj, eps)
    
    remove_reaction!(test, 12)
    @test_approx_eq length(model.rxns) length(test.rxns) + 1
    @test_approx_eq length(model.lb) length(test.lb) + 1
    @test_approx_eq length(model.ub) length(test.ub) + 1
    @test_approx_eq length(model.c) length(test.c) + 1
    @test_approx_eq length(model.rxn_rules) length(test.rxn_rules) + 1
    @test_approx_eq length(model.rxn_name) length(test.rxn_name) + 1
    @test_approx_eq length(model.rxn_subsystem) length(test.rxn_subsystem) + 1
    
    tmp.lb[8] = 0
    tmp.ub[8] = 0
    @test_approx_eq_eps(fba(tmp).obj,fba(remove_reaction(test, 8)).obj, eps)
    remove_reaction!(test, 8)
    @test_approx_eq_eps(fba(tmp).obj, fba(test).obj, eps)
    @test_approx_eq length(model.rxns) length(test.rxns) + 2
    @test_approx_eq length(model.lb) length(test.lb) + 2
    @test_approx_eq length(model.ub) length(test.ub) + 2
    @test_approx_eq length(model.c) length(test.c) + 2
    @test_approx_eq length(model.rxn_rules) length(test.rxn_rules) + 2
    @test_approx_eq length(model.rxn_name) length(test.rxn_name) + 2
    @test_approx_eq length(model.rxn_subsystem) length(test.rxn_subsystem) + 2    

    test = deepcopy(model)
    tmp = deepcopy(model)
    tmp.lb[12] = 0
    tmp.ub[12] = 0
    tmp.lb[8] = 0
    tmp.ub[8] = 0
    @test_approx_eq_eps(fba(tmp).obj,fba(remove_reaction(test, [test.rxns[8], test.rxns[12]])).obj, eps)
    remove_reaction!(test, [test.rxns[8], test.rxns[12]])
    @test_approx_eq_eps(fba(tmp).obj, fba(test).obj, eps)
    @test_approx_eq length(model.rxns) length(test.rxns) + 2
    @test_approx_eq length(model.lb) length(test.lb) + 2
    @test_approx_eq length(model.ub) length(test.ub) + 2
    @test_approx_eq length(model.c) length(test.c) + 2
    @test_approx_eq length(model.rxn_rules) length(test.rxn_rules) + 2
    @test_approx_eq length(model.rxn_name) length(test.rxn_name) + 2
    @test_approx_eq length(model.rxn_subsystem) length(test.rxn_subsystem) + 2    

    test = deepcopy(model)
    tmp = deepcopy(model)
    tmp.lb[12] = 0
    tmp.ub[12] = 0
    tmp.lb[8] = 0
    tmp.ub[8] = 0
    @test_approx_eq_eps(fba(tmp).obj,fba(remove_reaction(test, [8,12])).obj, eps)
    remove_reaction!(test, [8,12])
    @test_approx_eq_eps(fba(tmp).obj, fba(test).obj, eps)
    @test_approx_eq length(model.rxns) length(test.rxns) + 2
    @test_approx_eq length(model.lb) length(test.lb) + 2
    @test_approx_eq length(model.ub) length(test.ub) + 2
    @test_approx_eq length(model.c) length(test.c) + 2
    @test_approx_eq length(model.rxn_rules) length(test.rxn_rules) + 2
    @test_approx_eq length(model.rxn_name) length(test.rxn_name) + 2
    @test_approx_eq length(model.rxn_subsystem) length(test.rxn_subsystem) + 2    
end
function test_lp(lp, model)
    tmp = clone(model)
    tol = 1e-6
    @test_approx_eq_eps answer_lp(lp) fba(tmp).obj tol

    change_objective_coef(lp, 11)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 11)).obj tol
    change_objective_coef(lp, 28)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 28)).obj tol
    change_objective_coef(lp, 7)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 7)).obj tol
    change_objective_coef(lp, 13)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 13)).obj tol

    change_objective_coef(lp, 11, 0.5)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 11, 0.5)).obj tol
    change_objective_coef(lp, 28, 0.5)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 28, 0.5)).obj tol
    change_objective_coef(lp, 7, 0.5)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 7, 0.5)).obj tol
    change_objective_coef(lp, 13, 0.5)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 13, 0.5)).obj tol

    change_objective_coef(lp, [5,7])
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, [5,7])).obj tol
    @test get_objective(lp) == change_objective(tmp, [5,7]).c
    
    change_objective_coef(lp, [5,7],[1,1])
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, [5,7],[1,1])).obj tol
    @test get_objective(lp) == change_objective(tmp, [5,7],[1,1]).c
    
    change_objective_coef(lp, [5,7],[0.5,0.7])
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, [5,7],[0.5,0.7])).obj tol
    @test get_objective(lp) == change_objective(tmp, [5,7],[0.5,0.7]).c
    
    change_objective_coef(lp, [13,16],[0.5,0.7])
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, [13,16],[0.5,0.7])).obj tol
    @test get_objective(lp) == change_objective(tmp, [13,16],[0.5,0.7]).c
    
    change_objective_coef(lp, [7,5],[0.5,0.7])
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, [7,5],[0.5,0.7])).obj tol
    @test get_objective(lp) == change_objective(tmp, [7,5],[0.5,0.7]).c
    
    change_objective_coef(lp, 13)
    @test_approx_eq_eps answer_lp(lp) fba(change_objective(tmp, 13)).obj tol
    @test get_objective(lp) == change_objective(tmp, 13).c


    @test get_variable_bounds(lp) == (tmp.lb, tmp.ub)
    
    set_col_bounds(lp, 5, -500,500)
    @test get_variable_bounds(lp)[1] == (change_reaction_bounds(tmp, 5, -500, 500).lb)
    @test get_variable_bounds(lp)[2] == (change_reaction_bounds(tmp, 5, -500, 500).ub)

    set_col_bounds(lp, 5, 0,0)
    @test get_variable_bounds(lp)[1] == (change_reaction_bounds(tmp, 5, 0, 0).lb)
    @test get_variable_bounds(lp)[2] == (change_reaction_bounds(tmp, 5, 0, 0).ub)

    set_col_bounds(lp, 5, 0)
    @test get_variable_bounds(lp)[1] == (change_reaction_bounds(tmp, 5, 0).lb)
    @test get_variable_bounds(lp)[2] == (change_reaction_bounds(tmp, 5, 0).ub)

    set_col_bounds(lp, 5, tmp.lb[5], tmp.ub[5])
    @test_approx_eq_eps answer_lp(lp) fba(tmp).obj tol

    set_col_bounds(lp, 13, 0.5,0.5)
    @test_approx_eq_eps answer_lp(lp) fba(change_reaction_bounds(tmp, 13, 0.5, 0.5)).obj tol
    
    set_col_bounds(lp, 13, tmp.lb[13], tmp.ub[13])
    
    set_col_bounds(lp, 13, 0.5)
    @test_approx_eq_eps answer_lp(lp) fba(change_reaction_bounds(tmp, 13, 0.5)).obj tol
    
    set_col_bounds(lp, 13, tmp.lb[13], tmp.ub[13])
    @test_approx_eq_eps answer_lp(lp) fba(tmp).obj tol
end

function test_robustness_analysis(model)
    matsol = open_mat_file(Pkg.dir() * "/CBM/test/test_vars/robustness.mat")
    single = matsol["single"]
    double = matsol["double"]

    sol = robustness_analysis(model, [8]).result 
    @test_approx_eq sol single

    sol = robustness_analysis(model, [5,8]).result'
    @test_approx_eq sol double
end 
function test_fast_cc()
    testfile = open_mat_file(Pkg.dir() * "/CBM/test/test_vars/fastcctest.mat")
    A = testfile["testA"]

    tmp = load_matlab(Pkg.dir() * "/CBM/test/test_vars/consistRecon1.mat")

    @test A == fast_cc(tmp, 1e-3)
end
function test_fast_core()
    testfile = open_mat_file(Pkg.dir() * "/CBM/test/test_vars/fastcoretest.mat")
    C = testfile["testC"]
    A = testfile["testA"]

    tmp = load_matlab(Pkg.dir() * "/CBM/test/test_vars/consistRecon1.mat")

    @test A == fast_core(C, tmp, 1e-3)
end
function test_solvers(model)
    supported_solvers = ["GLPK", "CPLEX", "Gurobi", "Clp"]
    available_solvers = map(x -> any(x .== readdir(Pkg.dir())), supported_solvers)

    if available_solvers[1]
        test_lp(setup_glpk(model), model)
    end 
    if available_solvers[2]
        test_lp(setup_cplex(model), model)
    end 
    if available_solvers[3]
        test_lp(setup_gurobi(model), model)
    end 
    if available_solvers[4]
        test_lp(setup_clp(model), model)
    end 
end
