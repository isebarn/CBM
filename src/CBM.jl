module CBM
using Base.Test
solvers = map(x -> in(x, readdir(Pkg.dir())) ? x : "" , ["GLPK", "CPLEX", "Gurobi", "Clp"])
err_rd, err_wt = redirect_stderr()


using JSON
using ProgressMeter
using JuMP
using MAT
using HDF5
using Combinatorics

include("core/structure.jl")

if in("GLPK", solvers)
    using GLPK
    include("core/solvers/GLPK.jl")
end 

if in("CPLEX", solvers)
    using CPLEX
    include("core/solvers/CPLEX.jl")
end 

if in("Gurobi", solvers)
    using Gurobi
    include("core/solvers/Gurobi.jl")
end 

if in("Clp", solvers)
    using Clp
    using Clp.ClpCInterface
    include("core/solvers/Clp.jl")
end 

include("core/solvers/help.jl")

    export FBAsol 
    export MOMAsol 
    export Model 
    export SLG
    export RobustnessAnalysis 
    export FBAsolution 
    export SolverInfo
    export GeneKnockout

include("core/tools.jl")
    export abstract_gene_names
    export genes_in_reaction
    export disable_genes_in_rules
    export enable_genes_in_rules
    export convert_rule_to_boolean
    export convert_from_abstract
    export find_gene
    export gene_exists
    export find_gene_effects
    export genes_in_reaction_rule
    export setup_lp
    export setup_global_lp
    export min_max_lp
    export analyze_list_of_blocked_reaction_lp
    export answer_lp
    export construct_lp_solution
    export set_solver
    export delete_column
    export add_column
    export add_row
    export get_sparse_columns
    export clone
    export find_metabolite
    export metabolite_exists
    export find_reaction
    export reaction_exists
    export fix_empty_fields
    export fix_model
    export fix_extra_dicts
    export fix_gene_rules
    export fix_genes
    export initialize_cores
    export disable_cores
    export load_model
    export open_mat_file
    export save_mat_file

include("core/io.jl")
    export load_json
    export load_matlab
    export load_table
    export export_json
    export export_matlab
    export export_table

include("core/fba.jl")
    export fba
    export fva
    export fba_norm2
    export fba_norm1
    export fba_norm0
    export fva 
    export moma
    
    if in("GLPK", solvers)
        export settings_glpk
        export setup_glpk
    end 
    
    if in("CPLEX", solvers)
        export settings_cplex
        export setup_cplex
    end 
    
    if in("Gurobi", solvers)
        export setup_gurobi
    end 
    
    if in("Clp", solvers)
        export settings_clp
        export setup_clp
    end 

    export change_objective_coef
    export get_variable_bounds
    export get_constraint_bounds
    export get_constraint_matrix
    export get_direction
    export get_objective_value
    export get_reduced_costs
    export get_solution
    export get_solution_status
    export get_solution_status_code
    export get_slack
    export get_objective
    export set_direction
    export set_col_bounds
    export set_objective
    export solve
    export solver_status

include("core/simulations.jl")
    export find_blocked_reactions
    export synthetic_lethal_genes
    export robustness_analysis
    export find_deadend_metabolites
    export find_essential_genes
    export find_exchange_reactions
    export find_exchange_reactions
    export find_reactions_from_gene
    export find_reactions_from_metabolite
    export gene_deletion
    export knockout_genes
    export knockout_reactions
    export print_medium
    export print_reaction_formula
    export reaction_info

include("core/reconstruction.jl")
    export fast_cc
    export fast_core

# must bring PCBM into the scope of CBM
include("core/multicore.jl")
using PCBM

include("core/modification.jl")
    export add_metabolite!
    export add_metabolite
    export add_reaction!
    export add_reaction
    export change_objective!
    export change_objective
    export change_reaction_bounds!
    export change_reaction_bounds
    export convert_to_irreversible!
    export convert_to_irreversible
    export find_reversible
    export convert_to_reversible!
    export convert_to_reversible
    export reduce_model!
    export reduce_model
    export remove_gene!
    export remove_gene
    export remove_reaction!
    export remove_reaction

include("core/tests.jl")
    export test_add_reaction
    export test_change_objective
    export test_change_reaction_bounds
    export test_fba
    export test_find_blocked_reactions
    export test_find_deadend_metabolites
    export test_find_exchange_reactions
    export test_find_reactions_from_gene
    export test_find_reactions_from_metabolite
    export test_fva
    export test_gene_deletion
    export test_remove_reaction
    export test_lp
    export test_robustness_analysis
    export test_fast_cc
    export test_fast_core
    export test_solvers

# print errrors that arent method refefinition
outerr = String(readavailable(err_rd))
outerr = split(outerr, "\n")

for i in outerr
    if contains(i, "WARNING: Method definition")
        continue
    end  
    println(i)
end 

redirect_stderr(STDOUT)

end 
