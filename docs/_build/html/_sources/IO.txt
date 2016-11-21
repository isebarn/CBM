==
IO
==

* :ref:`loadjson`
* :ref:`loadmatlab`
* :ref:`loadtable`
* :ref:`trouble`

Loading a model
---------------

A model can be loaded from a *.json*, a *.mat* file or from a table (csv/tsv) using the methods::

	load_json(filename)
	load_matlab(filename)
	load_table(filename, [delimiter = ','])

.. _loadjson:

Loading JSON
============

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

.. _loadmatlab:

Loading Matlab
==============
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

.. _loadtable:

Load Table (csv/tsv)
====================
To load a model from a table, fx from a .csv file

::

	julia> model = load_table("./ecoli/e_coli", ',')
	           rxns :     95 Array{Any,1}
	           mets :     72 Array{Any,1}
	          genes :    137 Array{Any,1}
	              S :  13680 SparseMatrixCSC{Float64,Int64}
	             lb :     95 Array{Any,1}
	             ub :     95 Array{Any,1}
	              c :     95 Array{Any,1}
	              b :     72 Array{Any,1}
	         csense :      3 Dict{Any,Any}
	   rxn_gene_mat :      0 Array{Any,1}
	       rxn_name :     95 Array{Any,1}
	      rxn_rules :     95 Array{Any,1}
	  rxn_subsystem :     95 Array{Any,1}
	      rxn_extra :      0 Array{Any,1}
	    met_formula :     72 Array{Any,1}
	       met_name :     72 Array{Any,1}
	      met_extra :      0 Array{Any,1}
	      gene_name :    137 Array{Any,1}
	     gene_extra :      0 Array{Any,1}
	    description :      0 Array{Any,1}

**Note** that these table files need to be specifically structured:

Inside the "/ecoli/" directory there must be 4 files.

* e_coli_rxns.csv
* e_coli_mets.csv
* e_coli_genes.csv
* e_coli_constraints.csv

**Note** No additional fields are included in the model when loading from a tabular file

e_coli_rxns.csv
"""""""""""""""
With 8 columns, with the column contents as:

#. reaction id 
#. reaction lower bounds 
#. reaction upper bounds
#. objective coefficient, c, usually 0 for all, except 1 for the biomass function
#. reaction name
#. reaction gene rules
#. reaction subsystem
#. reaction formula
	**Note:** that if the ``<=``, ``=`` and ``=>`` are not wrapped with spacing, they wont be detected, so make sure they have spaces around them.

A line from ecoli_rxns.csv may look like::

	ETOHt2r,-1e3,1e3,0,ETOHt2r,,"Transport, Extracellular",etoh_e +  h_e <=>  etoh_c +  h_c



e_coli_mets.csv
"""""""""""""""

With 3 columns, with the column contents as:

1. metabolite id (as it appears in reaction formulas)
2. metabolite formula 
3. metabolite name 

A line from ecoli_mets.csv may look like::

	coa_c,C21H32N7O16P3S,Coenzyme A

e_coli_genes.csv
"""""""""""""""

With 2 columns, with the column contents as:

1. gene id (as it appears in the reaction gene rule)
2. gene name


A line from ecoli_genes.csv may look like::

	b0727,sucB


e_coli_constraints.csv
""""""""""""""""""""""

With 2 columns, with the column contents as:

1. b, usually an array of the same length as the metabolites, but filled with zeros
2. constraint sense, usually an array of the same length as the metabolites, but filled with ``=``
	
A line from ecoli_constraints.csv may look like::

	0,=

.. _trouble:

Troubleshooting
===============

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

A model can be exported in JSON/Matlab and Table format::

	export_json(model, filename);
	export_matlab(model, filename);

This will create a .json/.mat file, "ecoli.json"

To export model to a table, do::

	export_table(model, filename; [delimiter = ',', overwrite = false])

For example:: 

	export_table(model, "./ecoli/ecoli"; ',', false)

It is reccomended to store every tabular model in a private folder, as this method creates 4 files::

	./ecoli/ecoli_rxns.csv
	./ecoli/ecoli_mets.csv
	./ecoli/ecoli_genes.csv
	./ecoli/ecoli_constraints.csv

Separated by commas, ``,``. 

**Note:** By *default*, ``overwrite`` is set to ``false``, meaning that you cannot accidentally overwrite existing tabular files, the method will return an error.

However, if overwriting is desired, do::

	export_table(model, "./ecoli/ecoli"; ',', true)

