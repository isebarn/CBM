function fba(model::Model; direction::String = "max")
    fba(model, direction)
end 

function fba(model::Model, direction::String = "max")
	lp = setup_lp(model, direction);
	solve(lp)

	return construct_lp_solution(lp)
end

# -------------------------------------------------------------------
# -------------------------------------------------------------------


function fva(model::Model; optPercentage::Number = 100, flux_matrix::Bool = false)
	fva(model, optPercentage, flux_matrix)
end

function fva(model::Model, optPercentage::Number = 100, flux_matrix::Bool = false)
	
    if length(procs()) != 1
		return parallel_fva(model, optPercentage, flux_matrix)
	end 
    

	num_rxns = model.S.n
	objective_index = find(model.c)
	progress_meter = ProgressMeter.Progress(num_rxns*2, 0.1, "\n\u1b[1Ffva:", 50)

	lp = setup_lp(model, "max")

	# fix the biomass lower bound according to "optPercentage"
	solve(lp)	
	new_lb = get_solution(lp)[objective_index]
	new_lb *= (optPercentage/100)

	map((x,y) -> set_col_bounds(lp, x, y, model.ub[x]), objective_index, new_lb)
	map(x -> set_objective(lp, x, 0.0), objective_index)

	set_direction(lp, "min")
	minFlux, minV =  variability_runner(lp, num_rxns, progress_meter, flux_matrix)

	set_direction(lp, "max")
	maxFlux, maxV =  variability_runner(lp, num_rxns, progress_meter, flux_matrix)

    if flux_matrix
       return minFlux,maxFlux, minV, maxV
    else
	   return minFlux,maxFlux
    end 
end

function variability_runner(lp, num_rxns, progress_meter, flux_matrix)
	flux = Array(Float64, num_rxns)
	v = Array(Float64, num_rxns, num_rxns)
	for i = 1:num_rxns
		
		set_objective(lp,i)
		solve(lp)

		flux[i] = get_objective_value(lp)
		
        if flux_matrix 
		  v[:,i] = get_solution(lp)
        end 

		set_objective(lp,i, 0.0)
		next!(progress_meter)
	end

	return flux, v
end


function parallel_fva(model::Model, optPercentage::Number = 100, flux_matrix::Bool = false)
    num_rxns = length(model.c)

    objective = find(model.c)
    lb = deepcopy(model.lb)
    lb[objective] = fba(model).f * optPercentage/100

    extra_cores = procs()[2:end]

    map(x -> remotecall_fetch(setup_lp, x, 
            zeros(Float64, num_rxns),
            lb,
            model.ub,
            model.b,
            model.S,
            "max"), 
            procs())

    stepsize = div(num_rxns, 10)
    checklist = map(x -> collect(x:x+stepsize-1), 1:stepsize:num_rxns)
    checklist[end] = checklist[end][checklist[end] .<= num_rxns]

    if !flux_matrix
        maxF = pmap(solve_objective_obj, checklist)
        map(core -> remotecall_fetch(set_direction, core, "min"), extra_cores)
        minF = pmap(solve_objective_obj, checklist)
        
        return vcat(minF...), vcat(maxF...)
    else 
        maxFV = pmap(solve_objective_obj_sol, checklist)
        map(core -> remotecall_fetch(set_direction, core, "min"), extra_cores)
        minFV = pmap(solve_objective_obj_sol, checklist)

        maxF = vcat(map(x -> maxFV[x][1], 1:length(maxFV))...)
        maxV = hcat(map(x -> maxFV[x][2], 1:length(maxFV))...)

        minF = vcat(map(x -> minFV[x][1], 1:length(minFV))...)
        minV = hcat(map(x -> minFV[x][2], 1:length(minFV))...)

        return minF, maxF, minV, maxV
    end 
end 



# -------------------------------------------------------------------
# -------------------------------------------------------------------



function jump_fba(model)
  # Flux balance analysis
  #    sol=fba(model) performs flux balance Analysis
  #
  # Ref.  Savinell J, Palsson . J. Theor. Biol. 1992, 154(4):421-454.

  # THINK: Add side constraints, e.g. minimize 1 and 2 norms
  # 0-norm corresponds to min_set_rxns

  nrxns, nmets, _ = statistics(model)
  m, v, constr = create_fba_model(model)
  status = solve(m)
  info=SolverInfo(cobra_lp_solver, status)
  if status == :Optimal
    return FBAsol(getobjectivevalue(m), getvalue(v), getdual(constr), getdual(v), true, info)
  end
  return FBAsol(0, [], [], [], false, info)
end

function getmaxobjective(model)
  sol = fba(model)
  if !sol.success
    error("Unable to solve to optimality")
  end
  sol.obj
end

function fba_norm2(model::Model)
  # Flux balance analysis - minimum 2-norm solution

  optobj=getmaxobjective(model)
  nrxns, nmets, _ = statistics(model)
  m = createqpmodel()
  @variable(m, model.lb[i] <= v[i=1:nrxns] <= model.ub[i])
  @objective(m, Min, sum{v[i]^2, i=1:nrxns})
  @constraint(m, model.S*v .== model.m.b)
  @constraint(m, dot(v, model.c) >= optobj)
  status = solve(m)
  info = SolverInfo(cobra_qp_solver, status)
  if status == :Optimal
    return FBAsol(optobj, getvalue(v), [], [], true, info)
  end
  return FBAsol(0, [], [], [], false, info)
end

function fba_norm1(model::Model)
  # Flux balance analysis - minimum 1-norm solution

  # THINK:A lot of this is very similar to fba_norm2
  global cobra_warning_level
  optobj=getmaxobjective(model)
  nrxns, nmets, _ = statistics(model)
  m = createlpmodel()
  @variable(m, model.r.lb[i] <= v[i=1:nrxns] <= model.r.ub[i])
  @variable(m, z[i=1:nrxns])
  @objective(m, Min, sum(z))
  @constraint(m, model.s*v .== model.m.b)
  @constraint(m, dot(v, model.r.c) >= optobj)
  @constraint(m, z .>= v)
  @constraint(m, z .>= -v)
  status = solve(m)
  info = SolverInfo(cobra_lp_solver, status)
  if status == :Optimal
    return FBAsol(optobj, getvalue(v), [], [], true, info)
  else
    if cobra_warning_level
      warn("Unable to solve to optimality, status: $status")
    end
    return FBAsol(0, [], [], [], false, info)
  end
end

#function fba_norm0(model::Model, rxns)
function fba_norm0(model::Model)
  # Flux balance analysis - minimum number of nonzero fluxes (min 0-norm)

  # THINK: TBD Implement subset of reactions

  optobj=getmaxobjective(model)
  nrxns, nmets, _ = statistics(model)
  m = createmilpmodel()
  @variable(m, v[i=1:nrxns])
  @variable(m, y[i=1:nrxns], Bin) # THINK: subset
  @objective(m, Min, sum(y))
  @constraint(m, model.S*v .== model.b)
  @constraint(m, dot(v, model.c) >= optobj)
  @constraint(m, tieub[i=1:nrxns], v[i] <= model.ub[i]*y[i])
  @constraint(m, tielb[i=1:nrxns], model.lb[i]*y[i] <= v[i])
  status = solve(m)
  info=SolverInfo(cobra_mip_solver, status)
  if status == :Optimal # THINK: Suboptimal solutions?
    return FBAsol(optobj, getvalue(v), [], [], true, info)
  else
    if cobra_warning_level
      warn("Unable to solve to optimality, status: $status")
    end
    return FBAsol(0, [], [], [], false, info)
  end
end

#function fba_norm0(model::Model)=fba_norm0(model, collect(model.r.name))
function jump_fva(model::Model, percent=100)
  # Flux variability  analysis
  #    minFlux, maxFlux = fva(model, percent) performs flux variability analysis

  # Ref. Mahadevan R, Schilling C. Metab. Eng. 2003, 5(4):264-276.

  # THINK: Side constraints?

  global cobra_warning_level, cobra_tolerance
  nrxns, nmets, _ = statistics(model)
  m, x, constr = create_fba_model(model)
  status = solve(m)
  if status != :Optimal
    if cobra_warning_level
      warn("Unable to solve to optimality, status: $status")
    end
    return [], [], false, status
  end
  if percent > 0
    # Lower bound the cellular objective
    obj_lb = floor(getobjectivevalue(m) / cobra_tolerance) * cobra_tolerance * percent / 100.0;
    @constraint(m, dot(x, model.c) >= obj_lb)
  end

  minFlux, maxFlux = zeros(Float64, nrxns), zeros(Float64, nrxns)
  # Solve all maximization problems first, then the minimization problems
  for i=1:nrxns
    c = zeros(nrxns)
    c[i] = 1
    @objective(m, Max, dot(x, c))
    status = solve(m)
    maxFlux[i] = getobjectivevalue(m)
  end
  for i=1:nrxns
    c = zeros(nrxns)
    c[i] = 1
    @objective(m, Min, dot(x, c))
    status = solve(m)
    minFlux[i] = getobjectivevalue(m)
  end
  return minFlux, maxFlux, true, status
end

function moma(model_wt::Model, model_del::Model)
  # Minimization of metabolic adjustment (MOMA)
  #    moma(model_wt, model_del)

  # Ref: Segre et al. PNAS (2002), 99(23), 15112-15117

  # This formulation is identical to the one in the COBRA 2.0 toolbox. It
  # differs from Segre et al. in such a way that the wild-type
  # and mutant flux distributions are determined simultaneously. This is done
  # to avoid dependence on solver/solver settings which arises because of the
  # existence of alternative optima in the FBA LP.

  optobj=getmaxobjective(model)

  # Solve the min norm problem
  nrxns, nmets, _ = statistics(model) # THINK: Silly function
  m = createqpmodel()
  @variable(m, model_wt.r.lb[i] <= v_wt[i=1:nrxns] <= model_wt.r.ub[i])
  @variable(m, model_del.r.lb[i] <= v_del[i=1:nrxns] <= model_del.r.ub[i])
  @objective(m, Min, sum{(v_wt[i] - v_del[i])^2, i=1:nrxns})
  @constraint(m, model_wt.s*v_wt .== model_wt.m.b)
  @constraint(m, model_del.s*v_del .== model_del.m.b)
  @constraint(m, dot(model_wt.r.c, v_wt) == optobj)
  status = solve(m)
  info=SolverInfo(cobra_qp_solver, status)
  return MOMAsol(getobjectivevalue(m), getvalue(v_wt), getvalue(v_del), (status == :Optimal), info)
end

function create_fba_model(model)
  # Create FBA linear program
  nrxns, nmets, _ = statistics(model)
  m = createlpmodel()
  @variable(m, model.lb[i] <= v[i=1:nrxns] <= model.ub[i])
  @objective(m, Max, dot(v,model.c))
  @constraint(m, constr, model.S*v .== model.b)
  solve(m)
  return m, v, constr
end

