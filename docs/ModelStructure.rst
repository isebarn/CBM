Model Structure
===============

Loading a model
---------------

A model can be loaded either from a *.json* or a *.mat* file using the methods::

	load_json("file_location.json")
	load_matlab("file_location.mat")

.json
^^^^^

::

	julia> model = load_json("Models/e_coli_core.json")
	           rxns :     95 Array{String,1}
	           mets :     72 Array{String,1}
	          genes :    137 Array{String,1}
	              S :   6840 SparseMatrixCSC{Float64,Int64}
	             lb :     95 Array{Float64,1}
	             ub :     95 Array{Float64,1}
	              c :     95 Array{Float64,1}
	              b :     72 Array{Float64,1}
	         csense :      3 Dict{String,Array{Any,1}}
	   rxn_gene_mat :  13015 SparseMatrixCSC{Float64,Int64}
	       rxn_name :     95 Array{String,1}
	      rxn_rules :     95 Array{String,1}
	  rxn_subsystem :     95 Array{String,1}
	      rxn_extra :      3 Dict{String,Array{Any,1}}
	    met_formula :     72 Array{String,1}
	       met_name :     72 Array{String,1}
	      met_extra :      2 Dict{String,Array{Any,1}}
	      gene_name :    137 Array{String,1}
	     gene_extra :      1 Dict{String,Array{Any,1}}
	    description :     11     String


.mat
^^^^
::

	julia> model = load_matlab("Models/ecoli_core_model.mat")
	           rxns :     95 Array{String,1}
	           mets :     72 Array{String,1}
	          genes :    137 Array{String,1}
	              S :   6840 SparseMatrixCSC{Float64,Int64}
	             lb :     95 Array{Float64,1}
	             ub :     95 Array{Float64,1}
	              c :     95 Array{Float64,1}
	              b :     72 Array{Float64,1}
	         csense :      3 Dict{String,Array{Any,1}}
	   rxn_gene_mat :  13015 SparseMatrixCSC{Float64,Int64}
	       rxn_name :     95 Array{String,1}
	      rxn_rules :     95 Array{String,1}
	  rxn_subsystem :     95 Array{String,1}
	      rxn_extra :      3 Dict{String,Array{Any,1}}
	    met_formula :     72 Array{String,1}
	       met_name :     72 Array{String,1}
	      met_extra :      2 Dict{String,Array{Any,1}}
	      gene_name :    137 Array{String,1}
	     gene_extra :      1 Dict{String,Array{Any,1}}
	    description :     11     String

Troubleshooting
^^^^^^^^^^^^^^^

There may be some inconsistent models, so they wont be able to load properly.

.json
"""""

JSON-files are parsed as nested **dictionaries** in Julia. For a model in .json format to load properly, its structure must be consistent with the JSON-structure of models found in the BiGG model database

The following is the required structure for .json models::

	model: {} 
		reactions: []
			upper_bound        : number *
			lower_bound        : number *
			subsystem          : string
			name               : string
			gene_reaction_rule : string *
			id                 : string *
			metabolites:         {}		*
				metabolite : coefficient

		genes: []
			name : string 
			id   : string *

		metabolites: []
			formula     : string 
			compartment : string 
			name        : string
			id          : string *

Everything else goes into the fields rxn_extra, met_extra and gene_extra.

This structure must be followed for the model to load properly and the ones marked with * must be present for all functionality to be available

.mat
""""

Matlab files are saved as a dictionary of variables, different from .json in the sense that .mat files arent nested, i.e not dictionaries within dictionaries.

The following must be either present or provided in a .mat file::
	
	rxns
	mets
	genes
	S
	lb
	ub
	c
	b
	csense
	rxnGeneMat
	rxnNames
	grRules
	subSystems
	metFormulas
	metNames
	description

However, if fx in your model.mat file, "rxns" is saved as "rxn_ids" you can call::

	load_matlab(model, rxns = "rxn_ids")

To let the function know that "rxns" is called "rxn_ids" in your .mat file

There can be multiple corrections

	load_matlab(model, rxns = "rxn_ids", grRules = "rules")

csense
""""""

csense will be converted into a dictionary::

	julia> model.csense
	Dict{String,Array{Any,1}} with 3 entries:
	  "<=" => Any[]
	  "="  => Any[1,2,3 …  71,72]
	  ">=" => Any[]



csense is not strictly necessary, if csense is not present in a model, csense will be generated automatically 

Exporting a model
-----------------

A model can be exported in JSON/Matlab format::

	export_json(model, "ecoli");
	export_matlab(model, "ecoli");

This will create a .json/.mat file, "ecoli.json"
