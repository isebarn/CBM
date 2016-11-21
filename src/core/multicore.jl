using GLPK
using CPLEX
# 1. BASIC 
# 2. FVA 
# 3. GLPK

################### BASIC ####################
##############################################
#
# Basic functions 
# setup_lp()
# solve()
# get_objective_value()
# get_objective_value()
# change this to whichever solver you want to use 
# Note: only edit the global variable value
if !isdefined(:settings_default_lp)
    global settings_default_lp = "glpk"
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

function setup_lp{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, S::SparseMatrixCSC, objective::String = "max")
    if settings_default_lp == "glpk"
        setup_glpk(c,lb,ub,b,S,objective)
    elseif settings_default_lp == "cplex"
        setup_cplex(c,lb,ub,b,S,objective)
    else 
        warn("no solver chosen")
        return 
    end 
end

function solve()
    solve(lp)
end

function get_objective_value()
    get_objective_value(lp)
end

function get_solution()
    get_solution(lp)
end

function set_objective(objective::Number, value::Number = 1.0)
    set_objective(lp, objective, value)
end

function set_direction(direction::String)
    set_direction(lp, direction)
end

function set_col_bounds(index::Number, lb::Number, ub::Number)
    set_col_bounds(lp, index, lb, ub)
end

function set_col_bounds(index::Number, fixed_value::Number)
    set_col_bounds(lp, index, fixed_value)
end
###############################################

#################### FVA ######################
###############################################
# Necessary for FVA in parallel:
#   - solve_objective_obj(objective::Array{T})

#      this version is designed for FVA, returns the objective value
#   - solve
#   - get_objective_value
#   - set_objective
#   - solve_objective_obj_val(objective::Array{T})

#      this version is designed for FVA, returns the objective value and solution
#   - solve
#   - get_objective_value
#   - get_solution
#   - set_objective
function solve_objective_obj{T <: Number}(objective::Array{T})
    result = Array(Float64, length(objective))
    for (i,v) in enumerate(objective)
        set_objective(v)
        solve()
        result[i] = get_objective_value()
        set_objective(v, 0.0)
    end 
    return result
end 

function solve_objective_obj_sol{T <: Number}(objective::Array{T})
    result = Array(Float64, length(objective))
    solution = Array(Float64, GLPK.get_num_cols(lp), length(objective))
    for (i,v) in enumerate(objective)
        set_objective(v)
        solve()
        result[i] = get_objective_value()
        solution[:,i] = get_solution()
        set_objective(v, 0.0)
    end 
    return (result, solution)
end 
###############################################

############ GLPK implementations ##############
################################################
function get_objective_value(lp::GLPK.Prob)
    return GLPK.get_obj_val(lp)
end

function get_solution(lp::GLPK.Prob)  
    num_colums = GLPK.get_num_cols(lp)
    x = map((x -> get_col_prim(lp, x)), collect(1:num_colums))
end

function settings_glpk()
    P = SimplexParam()
    setfield!(P, Symbol("msg_lev"),convert(Int32, 1))
    #setfield!(P, Symbol("meth"),convert(Int32, X))
    #setfield!(P, Symbol("pricing"),convert(Int32, X))
    #setfield!(P, Symbol("r_test"),convert(Int32, X))
    #setfield!(P, Symbol("tol_bnd"),convert(Float64, X))
    #setfield!(P, Symbol("tol_dj"),convert(Float64, X))
    #setfield!(P, Symbol("tol_piv"),convert(Float64, X))
    #setfield!(P, Symbol("obj_ll"),convert(Float64, X))
    #setfield!(P, Symbol("obj_ul"),convert(Float64, X))
    #setfield!(P, Symbol("it_lim"),convert(Int32, X))
    #setfield!(P, Symbol("tm_lim"),convert(Int32, X))
    #setfield!(P, Symbol("out_frq"),convert(Int32, X))
    #setfield!(P, Symbol("out_dly"),convert(Int32, X))
    #setfield!(P, Symbol("presolve"),convert(Int32, X))
    return P
end

function set_col_bounds(lp::GLPK.Prob, index::Number, lb::Number, ub::Number)
    if lb == ub
        GLPK.set_col_bnds(lp, index, GLPK.FX, lb, ub)
    else 
        GLPK.set_col_bnds(lp, index, GLPK.DB, lb, ub)
    end
end

function set_col_bounds(lp::GLPK.Prob, index::Number, fixed_value::Number)
    GLPK.set_col_bnds(lp, index, GLPK.FX, fixed_value, fixed_value)
end

function set_direction(lp:: GLPK.Prob, direction::String)
    if lowercase(direction) == "max"
        GLPK.set_obj_dir(lp, GLPK.MAX);
    else
        GLPK.set_obj_dir(lp, GLPK.MIN);
    end
end

function set_objective(lp::GLPK.Prob, objective::Number, value::Number = 1.0)
    GLPK.set_obj_coef(lp, objective, value)
end

function setup_glpk{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, S::SparseMatrixCSC, objective::String = "max")
    global lp = GLPK.Prob()
    global glpk_settings = settings_glpk()
    GLPK.set_prob_name(lp, "sample")
    GLPK.term_out(GLPK.ON)
    
    GLPK.set_obj_dir(lp, GLPK.MAX)

    if contains(lowercase(objective),"min")
        GLPK.set_obj_dir(lp, GLPK.MIN)
    end

    GLPK.add_cols(lp, length(c))
    for i = 1:length(c)
        GLPK.set_col_bnds(lp, i, GLPK.DB, lb[i],ub[i])
        if lb[i] == ub[i]
            GLPK.set_col_bnds(lp, i, GLPK.FX, lb[i],ub[i])
        end
        GLPK.set_obj_coef(lp,i,c[i])
    end

    GLPK.add_rows(lp, length(b))
    for i = 1:length(b)
        GLPK.set_row_bnds(lp,i,GLPK.FX,b[i],b[i])
    end 

    GLPK.load_matrix(lp, S)
end 

function solve(lp::GLPK.Prob)
    return GLPK.simplex(lp, settings_glpk())
end 
###################################################

############ CPLEX implementations ##############
################################################
function change_objective{T <: Number}(lp::CPLEX.Model, objective::Array{T}, value::Array = [])
    if value == []
        value = ones(Float64, length(objective))
    end

    num_colums = CPLEX.num_var(lp)

    new_objective = zeros(num_colums)

    map((x,y) -> new_objective[x] = y, objective,value)

    CPLEX.set_obj!(lp, new_objective)
end

function get_solution(lp::CPLEX.Model)
    return CPLEX.get_solution(lp)
end

function settings_cplex()
    env = CPLEX.Env()

    CPLEX.set_param!(env, "CPX_PARAM_LPMETHOD", 1) # 0=automatic
    CPLEX.set_param!(env, "CPX_PARAM_SCRIND", 1) # 1=display messages
    CPLEX.set_param!(env, "CPX_PARAM_SIMDISPLAY", 0) # 0=simplex msg off

    #
        #CPLEX.set_param!(env, CPX_PARAM_ADVIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_AGGCUTLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_AGGFILL, X)
        #CPLEX.set_param!(env, CPX_PARAM_AGGIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_ALL_MAX, X)
        #CPLEX.set_param!(env, CPX_PARAM_ALL_MIN, X)
        #CPLEX.set_param!(env, CPX_PARAM_APIENCODING, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARCOLNZ, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARCROSSALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_BAREPCOMP, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARGROWTH, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARMAXCOR, X)
        #CPLEX.set_param!(env, CPX_PARAM_BAROBJRNG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARORDER, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARQCPEPCOMP, X)
        #CPLEX.set_param!(env, CPX_PARAM_BARSTARTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_BBINTERVAL, X)
        #CPLEX.set_param!(env, CPX_PARAM_BNDSTRENIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_BQPCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_BRDIR, X)
        #CPLEX.set_param!(env, CPX_PARAM_BTTOL, X)
        #CPLEX.set_param!(env, CPX_PARAM_CALCQCPDUALS, X)
        #CPLEX.set_param!(env, CPX_PARAM_CLIQUES, X)
        #CPLEX.set_param!(env, CPX_PARAM_CLOCKTYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_CLONELOG, X)
        #CPLEX.set_param!(env, CPX_PARAM_COEREDIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_COLREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_CONFLICTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_CONFLICTDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_COVERS, X)
        #CPLEX.set_param!(env, CPX_PARAM_CPUMASK, X)
        #CPLEX.set_param!(env, CPX_PARAM_CRAIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTLO, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTPASS, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTSFACTOR, X)
        #CPLEX.set_param!(env, CPX_PARAM_CUTUP, X)
        #CPLEX.set_param!(env, CPX_PARAM_DATACHECK, X)
        #CPLEX.set_param!(env, CPX_PARAM_DEPIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_DETTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_DISJCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_DIVETYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_DPRIIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_EACHCUTLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPAGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPINT, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPLIN, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPMRK, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPOPT, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPPER, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPRELAX, X)
        #CPLEX.set_param!(env, CPX_PARAM_EPRHS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FEASOPTMODE, X)
        #CPLEX.set_param!(env, CPX_PARAM_FILEENCODING, X)
        #CPLEX.set_param!(env, CPX_PARAM_FLOWCOVERS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FLOWPATHS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FPHEUR, X)
        #CPLEX.set_param!(env, CPX_PARAM_FRACCAND, X)
        #CPLEX.set_param!(env, CPX_PARAM_FRACCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_FRACPASS, X)
        #CPLEX.set_param!(env, CPX_PARAM_GUBCOVERS, X)
        #CPLEX.set_param!(env, CPX_PARAM_HEURFREQ, X)
        #CPLEX.set_param!(env, CPX_PARAM_IMPLBD, X)
        #CPLEX.set_param!(env, CPX_PARAM_INTSOLFILEPREFIX, X)
        #CPLEX.set_param!(env, CPX_PARAM_INTSOLLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_ITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_LANDPCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_LBHEUR, X)
        #CPLEX.set_param!(env, CPX_PARAM_LOCALIMPLBD, X)
        #CPLEX.set_param!(env, CPX_PARAM_LPMETHOD, X)
        #CPLEX.set_param!(env, CPX_PARAM_MCFCUTS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MEMORYEMPHASIS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPCBREDLP, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPEMPHASIS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPINTERVAL, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPKAPPASTATS, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPORDIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPORDTYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_MIPSEARCH, X)
        #CPLEX.set_param!(env, CPX_PARAM_MPSLONGNUM, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETEPOPT, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETEPRHS, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETFIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_NETPPRIIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_NODEFILEIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_NODELIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_NODESEL, X)
        #CPLEX.set_param!(env, CPX_PARAM_NUMERICALEMPHASIS, X)
        #CPLEX.set_param!(env, CPX_PARAM_NZREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_OBJDIF, X)
        #CPLEX.set_param!(env, CPX_PARAM_OBJLLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_OBJULIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_OPTIMALITYTARGET, X)
        #CPLEX.set_param!(env, CPX_PARAM_PARALLELMODE, X)
        #CPLEX.set_param!(env, CPX_PARAM_PERIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_PERLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERDETTIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTEREPAGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTEREPGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERINTSOL, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERNODE, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHAFTERTIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_POLISHTIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_POPULATELIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_PPRIIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_PREDUAL, X)
        #CPLEX.set_param!(env, CPX_PARAM_PREIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_PREPASS, X)
        #CPLEX.set_param!(env, CPX_PARAM_PROBETIME, X)
        #CPLEX.set_param!(env, CPX_PARAM_QPMAKEPSDIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_QPMETHOD, X)
        #CPLEX.set_param!(env, CPX_PARAM_QPNZREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_QTOLININD, X)
        #CPLEX.set_param!(env, CPX_PARAM_RAMPUPDETTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_RAMPUPDURATION, X)
        #CPLEX.set_param!(env, CPX_PARAM_RAMPUPTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_REINV, X)
        #CPLEX.set_param!(env, CPX_PARAM_RELAXPREIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_RELOBJDIF, X)
        #CPLEX.set_param!(env, CPX_PARAM_REPAIRTRIES, X)
        #CPLEX.set_param!(env, CPX_PARAM_RINSHEUR, X)
        #CPLEX.set_param!(env, CPX_PARAM_ROWREADLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SCAIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_SCRIND, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIFTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIFTDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIFTITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SIMDISPLAY, X)
        #CPLEX.set_param!(env, CPX_PARAM_SINGLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SOLNPOOLAGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_SOLNPOOLGAP, X)
        #CPLEX.set_param!(env, CPX_PARAM_SOLUTIONTYPE, X)
        #CPLEX.set_param!(env, CPX_PARAM_STARTALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_STRONGITLIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_SUBALG, X)
        #CPLEX.set_param!(env, CPX_PARAM_SYMMETRY, X)
        #CPLEX.set_param!(env, CPX_PARAM_THREADS, X)
        #CPLEX.set_param!(env, CPX_PARAM_TILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_TRELIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_TUNINGREPEAT, X)
        #CPLEX.set_param!(env, CPX_PARAM_TUNINGTILIM, X)
        #CPLEX.set_param!(env, CPX_PARAM_VARSEL, X)
        #CPLEX.set_param!(env, CPX_PARAM_WORKDIR, X)
        #CPLEX.set_param!(env, CPX_PARAM_WORKMEM, X)
        #CPLEX.set_param!(env, CPX_PARAM_WRITELEVEL, X)
        #CPLEX.set_param!(env, CPX_PARAM_ZEROHALFCUTS, X)
    return env
end

function set_col_bounds(lp::CPLEX.Model, index::Number, lb::Number, ub::Number)
    new_lb = CPLEX.get_varLB(lp)
    new_lb[index] = lb
    CPLEX.set_varLB!(lp, new_lb)

    new_ub = CPLEX.get_varUB(lp)    
    new_ub[index] = ub
    CPLEX.set_varUB!(lp, new_ub)
end

function set_col_bounds(lp::CPLEX.Model, index::Number, fixed_value::Number)
    new_lb = CPLEX.get_varLB(lp)
    new_lb[index] = fixed_value
    CPLEX.set_varLB!(lp, new_lb)

    new_ub = CPLEX.get_varUB(lp)    
    new_ub[index] = fixed_value
    CPLEX.set_varUB!(lp, new_ub)
end

function set_direction(lp::CPLEX.Model, direction::AbstractString)
    if lowercase(direction) == "max"
        CPLEX.set_sense!(lp, :Max)
    else
        CPLEX.set_sense!(lp, :Min)
    end
end

function set_objective(lp::CPLEX.Model, objective::Number, value::Number = 1.0)
    new_objective = CPLEX.get_obj(lp)
    new_objective[objective] = value 

    CPLEX.set_obj!(lp, new_objective)
end

function setup_cplex{T <: Number}(c::Array{T}, lb::Array{T}, ub::Array{T}, b::Array{T}, S::SparseMatrixCSC, objective::AbstractString = "max")
  # Flux Balance Analysis usingh CPLEX.jl

    nmets = length(b)

    # Construct the CPLEX model
    env = CPLEX.Env()
    lp = CPLEX.Model(env, "CPLEX")
    CPLEX.set_param!(env, "CPX_PARAM_LPMETHOD", CPLEX.CPX_ALG_AUTOMATIC) # should?


    if objective == "max"
        set_direction(lp, "max")
    else
        set_direction(lp, "min")
    end 
    
    CPLEX.add_vars!(lp, c, lb, ub)

    CPLEX.add_constrs!(lp, S, fill('=', nmets), b)


    return lp
end

function solve(lp::CPLEX.Model)
    CPLEX.optimize!(lp)
end
###################################################

