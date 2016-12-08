using CBM
using Base.Test


function test_all()
    model = load_json(Pkg.dir() * "/CBM/Models/e_coli_core.json")

    println("test_add_reaction")
    test_add_reaction(model)
    println("test_change_objective")
    test_change_objective(model)
    println("test_change_reaction_bounds")
    test_change_reaction_bounds(model)
    println("test_fba")
    test_fba(model)
    println("test_find_blocked_reactions")
    test_find_blocked_reactions(model)
    println("test_find_deadend_metabolites")
    test_find_deadend_metabolites(model)
    println("test_find_exchange_reactions")
    test_find_exchange_reactions(model)
    println("test_find_reactions_from_gene")
    test_find_reactions_from_gene(model)
    println("test_find_reactions_from_metabolite")
    test_find_reactions_from_metabolite(model)
    println("test_fva")
    test_fva(model)
    #println("test_gene_deletion")
    #test_gene_deletion(model)
    println("test_remove_reaction")
    test_remove_reaction(model)
    #println("test_robustness_analysis")
    #test_robustness_analysis(model)
    #println("test_solvers")
    #test_solvers(model)
end 

test_all()