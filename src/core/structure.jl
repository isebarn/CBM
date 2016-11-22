type GeneDeletion
    r_g
    r_f 
    g_f 
end

function Base.show(io::IO, m::GeneDeletion)
        println("\nType: GeneDeletion")
        @printf "%30s %7s : %6s %10s\n" "" "Field" "Size" "Fieldtype"
        @printf "%30s %7s : %6s %10s\n" "Reactions => Gene-knockout" fieldnames(m)[1] length(getfield(m,1)) typeof(getfield(m,1))
        @printf "%30s %7s : %6s %10s\n" "Reactions => Flow" fieldnames(m)[2] length(getfield(m,2)) typeof(getfield(m,2))
        @printf "%30s %7s : %6s %10s\n" "Gene-knockout => Flow" fieldnames(m)[3] length(getfield(m,3)) typeof(getfield(m,3))
end
type LPSolution
	f
	x
	y
	w
	objective 
	solution 
	slack 
	rcost
	status

	function LPSolution(f,x,y,w,status)
		new(f,x,y,w,f,x,y,w,status)
	end 
end

function Base.show(io::IO, m::LPSolution)
	@printf "%1s: %3s\n"  "LPSolution" m.status

	@printf "%17s %9f\n" "objective::" m.f
	@printf "%17s %3d %2s\n" "flux::" length(m.x) "element array"
	@printf "%17s %3d %2s\n" "slack::" length(m.y) "element array"
	@printf "%17s %3d %2s\n" "rcosts::" length(m.w) "element array"
end
# THINK: Specify types wherever appropriate, remove obsolete stuff (sparsestruct?)
type Model
	rxns
	mets
	genes
	S
	lb
	ub
	c
	b
	csense
	rxn_gene_mat
	rxn_name
	rxn_rules
	rxn_subsystem
	rxn_extra
	met_formula
	met_name
	met_extra
	gene_name
	gene_extra
	description
end

function Base.show(io::IO, m::Model)
	for (i,v) in enumerate(fieldnames(m))
		@printf "%15s : %6s %10s\n" v length(getfield(m,i)) typeof(getfield(m,i))
	end
end

Base.length(m::Model) = (length(m.c), length(m.b), length(m.genes))
Base.size(m::Model) = (length(m.c), length(m.b), length(m.genes))

type Robustness
    result 
    reactions
    ranges 
end

function Base.show(io::IO, m::Robustness)
    @printf "%1s \n"  "Robustness Analysis"
    @printf "%30s %15s \n"  "result :" string(size(m.result)," Array")
    
    @printf "%30s \n"  "ranges :"
    @printf "%38s : %15s \n" "reaction" "range:"
    for (i,v) in enumerate(m.ranges)
        @printf "%38s : %15s \n"  m.reactions[i] string("(", round(v[1],5), ",", round(v[end],5), ")")
    end 
    
end

Base.range(a::Robustness, minFlux, maxFlux) = map(y -> ind2sub(a.result, y), find(x -> (x >= minFlux) & (x <= maxFlux), a.result))
Base.getindex(a::Robustness, idx...) = a.result[idx...]
Base.size(a::Robustness) = size(a.result)
