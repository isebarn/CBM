"""
    fast_core(C, model, [eps = 1e-3])
The FASTCORE algorithm for context-specific metabolic network reconstruction
C is the core set
returns the reconstruction

### Examplex
    A = fast_core(C, model)

Return the indices of the reactions in the new model
"""
function fast_core(C::Array{Int64,1}, org_model::Model, eps::Number = 1e-3)
    model = clone(org_model)
    num_rxns = length(model.c)
    N = collect(1:num_rxns)
    I = find(model.lb .>= 0.0)

    A = Int64[]
    flipped = false
    singleton = false 

    J = sort(intersect(I,C))

    P = sort(setdiff(N, C))

    Supp = find_sparse_mode( J, P, singleton, model, eps);


    if !isempty(setdiff(J,Supp))
        println("Error: Inconsistent irreversible core reactions.\n")
        return 
    end 

    A = deepcopy(Supp)
    J = sort(setdiff(C,A))

    while !isempty(J)
        P = sort(setdiff(P,A))

        Supp = find_sparse_mode( J, P, singleton, model, eps);
        A = sort(union(A, Supp))

        if !isempty(intersect(J,A))
            J = sort(setdiff(J,A))
            flipped = false
        else 
            if singleton
                JiRev = sort(setdiff(J[1], I))
            else 
                JiRev = sort(setdiff(J, I))
            end

            if (flipped | isempty(JiRev))
                if singleton
                    return
                else 
                    flipped = false 
                    singleton = true 
                end 
            else 
                model.S[:,JiRev] = -model.S[:,JiRev]
                tmp = deepcopy(model.ub[JiRev])
                model.ub[JiRev] = -model.lb[JiRev]
                model.lb[JiRev] = -tmp
                flipped = true 
            end 
        end 
    end
    return A
end

function find_sparse_mode(J, P, singleton, model::Model, eps::Number = 1e-3)
    if isempty(J)
        return 
    end 

    if singleton
        V = LP7(J[1], model, eps)
    else 
        V = LP7(J, model, eps)
    end 
    K = sort(intersect(find(V .>= .99*eps), J))
 
    if isempty(K)
        return Int64[]
    end 

    V = LP9(K, P, model, eps)

    return find(abs(V) .>= .99 * eps)
end

function LP9(K, P, model::Model, eps::Number = 1e-3)
    scalingfactor = 1e5

    V = Int64[]

    if isempty(P) | isempty(K)
        return Int64[]
    end 

    np = length(P)    
    nk = length(K)

    m,n = size(model.S)

    f = [zeros(n); ones(np)]

    Aeq = deepcopy(model.S)
    for i in 1:np
        Aeq = add_column(Aeq)
    end 
    beq = zeros(m)

    Ip = spzeros(np,n)
    for (i,v) in enumerate(P)
        Ip[i,v] = 1
    end

    Ik = spzeros(nk,n)
    for (i,v) in enumerate(K)
        Ik[i,v] = 1
    end

    Aineq = vcat(hcat(Ip, -speye(np)), hcat(-Ip, -speye(np)), hcat(-Ik, spzeros(nk, np)))
    bineq = vcat(zeros(2*np), -ones(nk) * eps * scalingfactor)

    lb = vcat(model.lb, zeros(np)) * scalingfactor 
    ub = vcat(model.ub, max(abs(model.ub[P]), abs(model.lb)[P]))  * scalingfactor

    if settings_default_lp == "cplex"
        lp = setup_cplex(f,lb,ub,beq,Aeq,bineq,Aineq, "min")
    elseif settings_default_lp == "glpk"
        lp = setup_glpk(f,lb,ub,beq,Aeq,bineq,Aineq, "min")
    end 
    solve(lp)
    x = get_solution(lp)

    return x[1:n]
end 

function LP7(J, model::Model, eps::Number = 1e-3)
    nj = length(J)
    m,n = size(model.S)

    f = -[zeros(n); ones(nj)]
    
    Aeq = deepcopy(model.S)
    for i in 1:nj
        Aeq = add_column(Aeq)
    end 

    beq = zeros(m)

    Ij = spzeros(nj, n)
    for (i,v) in enumerate(J)
       Ij[i,v] = -1
    end


    Aineq = hcat(Ij, speye(nj))
    bineq = zeros(nj)

    lb = [model.lb; zeros(Float64,nj)]
    ub = [model.ub; ones(Float64, nj)* eps]
    if settings_default_lp == "cplex"
        lp = setup_cplex(f,lb,ub,beq,Aeq,bineq,Aineq, "min")
    elseif settings_default_lp == "glpk"
        lp = setup_glpk(f,lb,ub,beq,Aeq,bineq,Aineq, "min")
    end 
    solve(lp)
    x = get_solution(lp)
    return x[1:n]
end



"""
    fast_cc(model, [eps = 1e-3])

The FASTCC algorothm to test the consistency of a stochiometric model

### Examples
    A = fast_cc(model)

Return a boolean vector marking thos reactions having consistent flux
"""
function fast_cc(org_model::Model, eps::Number = 1e-3)
    model = deepcopy(org_model)
    num_rxns = length(model.c)

    N = collect(1:num_rxns)
    I = find(model.lb .>= 0.0)

    A = []

    J = sort(intersect(N, I))

    V = LP7(J, model, eps)

    Supp = find(abs(V) .>= 0.99*eps)

    A = deepcopy(Supp)

    incI = sort(setdiff(J,A))

    if !isempty( incI )
        println("inconsistent subset of I detected");
    end

    J = sort(setdiff( sort(setdiff(N,A)), incI))


    flipped = false
    singleton = false

    while !isempty(J)

        if singleton 
            Ji = J[1]
            V = LP3(Ji, model)
        else 
            Ji = J 
            V = LP7(Ji, model, eps)
        end 

        Supp = find( abs(V) .>= 0.99*eps)
        A = sort(union(A, Supp))

        
        if !isempty(intersect(A,J))
            J = sort(setdiff(J,A))
            flipped = false 

        else
            JiRev = sort(setdiff(Ji,I))
            if flipped | isempty(JiRev)
                flipped = false 
                if singleton 
                    J = sort(setdiff(J, Ji))
                else 
                    singleton = true 
                end 
            else 
                model.S[:,JiRev] = -model.S[:,JiRev]
                tmp = deepcopy(model.ub[JiRev])
                model.ub[JiRev] = -model.lb[JiRev]
                model.lb[JiRev] = -tmp
                flipped = true 
            end 
        end
    end 
    return A
end

function LP3(J, model::Model, eps::Number = 1e-3)
    m,n = size(model.S)

    f = zeros(n)
    f[J] = -1

    Aeq = deepcopy(model.S)
    beq = zeros(m)

    lb = deepcopy(model.lb)
    ub = deepcopy(model.ub)

    lp = setup_cplex(f,lb,ub,beq,Aeq, "min")
    solve(lp)
    x = get_solution(lp)
    return x
end 


