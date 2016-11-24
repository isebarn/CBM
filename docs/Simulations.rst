===========
Simulations
===========

.. highlight:: julia

CBM offers many methods inteded for model analysis.

**Note:** arguments inside brackets ``[ ]`` denote **keyword arguments**, for example 
``flux balance analysis`` is defined as ::

	fba(model; [direction = "max"])

meaning you can equivalently call::

	fba(model)
	fba(model, "max")
	fba(model, direction = "max")




* :ref:`fba`
* :ref:`fva`
* :ref:`find_blocked_reactions`
* :ref:`synthetic_lethal_genes`
* :ref:`robustness_analysis`
* :ref:`find_deadend_metabolites`
* :ref:`find_essential_genes`
* :ref:`find_exchange_reactions`
* :ref:`find_reactions_from_gene`
* :ref:`find_reactions_from_metabolite`
* :ref:`gene_deletion`
* :ref:`knockout_genes`
* :ref:`knockout_reactions`
* :ref:`print_reaction_formula`
* :ref:`reaction_info`

.. _fba:

Flux Balance Analysis
---------------------

``fba()`` is a mathematical method for simulating metabolism in genome-scale reconstructions of metabolic networks, by optimizing the network w.r.t (usually) the reactions responsible for the organisms growth::

	fba(model; [direction = "max", objective = 0])

* ``direction`` may be either ``"max"`` or ``"min"``
* ``objective`` may be either an **integer index** or **reaction name** as it appears in ``model.rxns``

**Example**

To maximize the default objective function::

	julia> fba(model)
	FBAsolution: 
	       obj::  0.873922
	         v::  95 element-array
	     slack::  72 element-array
	    rcosts::  95 element-array
	   success::  Optimal
	      info::  SolverInfo("glpk","GLP_OPT")

To maximize "ADK1"::

	julia> fba(model, objective = "ADK1")
	FBAsolution: 
	       obj:: 166.610000
	         v::  95 element-array
	     slack::  72 element-array
	    rcosts::  95 element-array
	   success::  Optimal
	      info::  SolverInfo("glpk","GLP_OPT")

To maximize the reaction number 14::

	julia> fba(model, objective = 14)
	FBAsolution: 
	       obj:: 11.104242
	         v::  95 element-array
	     slack::  72 element-array
	    rcosts::  95 element-array
	   success::  Optimal
	      info::  SolverInfo("glpk","GLP_OPT")


To minimize, use:: 

	julia> fba(model, "min")
	FBAsolution: 
	       obj::  0.000000
	         v::  95 element-array
	     slack::  72 element-array
	    rcosts::  95 element-array
	   success::  Optimal
	      info::  SolverInfo("glpk","GLP_OPT")


or ``fba(model, objective = "min")``

Where 

* ``sol.obj`` the objective value
* ``sol.v`` represents the solution vector
* ``sol.slack`` represents the slack
* ``sol.rcosts`` represents the reduced costs
* ``sol.info`` details solver information


.. _fva:

Flux Variability Analysis
-------------------------

Returns the **minimum** and **maximum** flux of every reaction in the model

**Note:** This method may run **in parallel** ::

	fva(model::Model, [optPercentage = 1, flux_matrix = false])

* ``optPercentage``: Fix the lower bound of the biomass reaction (or objective reaction) to a fraction of its maximum possible value.
* ``flux_matrix``: In addition to **minimum** and **maximum** fluxes, return the entire solution flux for every reaction.

**Example**

To calculate the minimum and maximum flux values of **every reaction** with biomass fixed
at 50% of its maximum value::

	minFlux, maxFlux = fva(model, 0.5)

To calculate the flux values of **every** reaction for every reactions minimum flux and maximum flux, call::

	minFlux, maxFlux, minFluxArray, maxFluxArray = fva(model, 1, true)


.. _find_blocked_reactions:

Find blocked reactions
----------------------

Locate every reaction that is constrianed to **zero flux** in every case::

	blocked_reactions = (model::Model, [tolerance::Number = 1e-10])

Where tolerance represents how close to absolute 0.0 the reaction flux must lie 

**Example**

For ecoli core::

	julia> find_blocked_reactions(model)
	8-element Array{Any,1}:
	 "EX_fru_e"   
	 "EX_fum_e"   
	 "EX_gln__L_e"
	 "EX_mal__L_e"
	 "FRUpts2"    
	 "FUMt2_2"    
	 "GLNabc"     
	 "MALt2_2" 

.. _synthetic_lethal_genes:

Synthetic Lethal Genes
----------------------

Check the model for essential genes, conditionally essential genes and non-essential genes

Returns a ``SLG`` type, which contains the fields ``ess``, ``cond_ess`` and ``non_ess``,::

    synthetic_lethal_genes(model, cutoff = 0.1, num_runs = 100)

where ``cutoff`` represents the minimum biomass flux as a fraction of the wild-type flux. 

``num_runs`` indicates how many times the algorithm runs, higher number gives better results, but takes longer.

Example
^^^^^^^


To find essentiality with biomass fixed at ``0.1`` in 200 runs::

	julia> slg = synthetic_lethal_genes(model, 0.1, 200)
	run 200 of 200 Number of conditionally essential genes found: 115

	Type: SLG
	     Field                  Contains    Size  
	       ess                Essential:       7  
	  cond_ess    Conditional Essential:     115  
	   non_ess            Non-Essential:      15 


``ess``, **essential genes**
""""""""""""""""""""""""""""

Those genes, that if disabled, make biomass production ratio to the wild type drop below 10%. The field slg.ess is a ``Dict()`` containing the gene indices and the resulting biomass production  after their deletion

``cond_ess``, **conditionally essential** 
"""""""""""""""""""""""""""""""""""""""""
Those genes, that if disabled along with some other genes will make the biomass production ratio to the wild type drop below 10%. The field slg.cond_ess is a ``Dict()`` containing the conditionally-essential genes and **all** the gene-combinations that were checked that reduced the biomass below 10%.
Example entry in this dictionary may be ::

	90  => Any[Any[89,74,12,79,67], Any[92,95,67]]

meaning that disabling gene 90 along with **either** genes 89,74,12,79,67 **or** 92,95,67, brought the biomass below 10%

``non_ess``, **non-essential genes** 
""""""""""""""""""""""""""""""""""""

Genes that never have the effect of making the biomass production ratio to the wild type ratio drop below 10%. This is an array containing the indices of the genes that never brought the biomass below 10%

.. _robustness_analysis:

Robustness Analysis
-------------------

Perform a robustness analysis for any number of fixed reactions, specified by indices or reaction names::

	solution = robustness_analysis(model, reactions, [objective = 0, pts = [], direction = "max"])

The method can either be called with **reaction names** or with **reaction indices.**

* ``objective`` reaction can be chosen, either as the **reaction name**, such as ``"BIOMASS_Ecoli_core_w_GAM"`` or as an index, such as ``13``. If left blank it defaults to the default objective of the model
* ``pts`` is an array used to specify how many points are tested for each reaction, if left blank, it creates 20 points for each reaction
** ``direction`` represents the optimization direction, either ``"min"`` or ``"max"`` which is the default.

Examples
^^^^^^^^

Robustness analysis for "ACONTb" at 20 points 
"""""""""""""""""""""""""""""""""""""""""""""

To see how the biomass of **e_coli_core** if the flow of "ACONTb" is fixed at 8 different points between its minimum and maximum::

julia> robust_sol = robustness_analysis(model, ["ACONTb"], "BIOMASS_Ecoli_core_w_GAM", [8])
Robustness Analysis 
                      result :      (8,) Array 
                      ranges : 
                              reaction :          range: 
                                     5 :     (-0.0,20.0) 

We see that "ACPNTb" (reaction number 5) has a minimum flux of 0.0, and maximum flux of 20.0. To view the objective flux values type ``robust_sol.result`` or robust_sol[:]::

	julia> robust_sol.result
	8-element Array{Float64,1}:
	  6.65799e-27
	  0.850845   
	  0.871775   
	  0.809412   
	  0.607059   
	  0.404706   
	  0.202353   
	 -1.52026e-16

and to plot (if ``Plots.jl`` is installed)::

	plots(robust_sol[:])


.. figure:: https://raw.githubusercontent.com/isebarn/CBM/master/docs/pics/acontb_biomass_8pts.png
	:scale: 20

.. image:: CBM/docs/pics/acontb_biomass_8pts.png
   :scale: 50

To perform a robustness analysis on reaction 13 againts reactions 5,8 and 11, where the point resolution
for reactions 5,8 and 11 is 4,2 and 10, respectively::

    robustness_analysis(model, [5,8,11], 13, [4,2,10], "max")
	Robustness Analysis 
	                      result :   (2,3,4) Array 
	                      ranges : 
	                              reaction :          range: 
	                                     5 :     (-0.0,20.0) 
	                                     8 :      (0.0,20.0) 
	                                    11 :    (8.39,175.0) 


The method returns a **Robustness** object with fields 

 * ``solution.result``

    a multidimensional matrix which contains the flux value of the objective for each point. 
 
 * ``solution.reactions``

    names of the reactions

 * ``solution.ranges``

    the flux values of the reactions at specific points ::

        solution.ranges[1]
        4-element Array{Float64,1}:
         -1.3884e-28
          6.66667   
         13.3333    
         20.0 

    So reaction 1 (reaction 5, ACONTb) is fixed at 6.6666 at point 2 

To see the value of the objective, reaction 13, when reaction 5 
is fixed at point 1, reaction 8 at point 2 and reaction 11 at 
point 10::

    solution[1,2,10]
    6.365820260756199e-16

To see every point where the objective is between 0.5 and 0.9::
    
    range(solution, 0.5, 0.9)
    4-element Array{Tuple{Int64,Int64,Int64},1}:
     (1,2,3)
     (2,2,3)
     (3,2,3)
     (4,2,3)

To get the entire vector for reaction 5 when reaction 8 is fixed at 
point 2 and reaction 11 is fixed at point 6::

    solution[:,2,6]
    4-element Array{Float64,1}:
     0.425313 
     0.314254 
     0.203194 
     0.0921351  

Return a matrix with the objective reactions optimal value at each  
value of the control reaction

**Plotting**

It is reccomended to use PyPlot to plot in Julia.

using PyPlot, to plot a 2D image of the effect of reaction 5
while reaction 8 is fixed at point 2 and reaction 11 at point 6::

    PyPlot.plot(a[:,2,6])

to plot a 3D image in PyPlot for reactions 5 and 8 while 
reaction 11 is fixed at point 6::

    PyPlot.surf(a[:,:,6])
    
.. _find_deadend_metabolites:

Find Dead End Metabolites
-------------------------

Return those metabolites that are **only produced** and **only consumed**::

	find_deadend_metabolites(model, [lbfilter = true])

By default, those metabolites with a positive lower bound are filtered out of the results

**Example** ::

	julia> idx = find_deadend_metabolites(model)
	5-element Array{Int64,1}:
	 30
	 32
	 37
	 49
	 60


.. _find_essential_genes:

Find Essential Genes
--------------------

Return a list of those genes, which would, if knocked out, bring the objective flux beneath that
specified by ``min_growth`` ::

	essential = find_essential_genes(model, min_growth)

**Example**

Find those genes that would bring biomass production beneath 0.01::

    julia> find_essential_genes(model,0.01)
    2-element Array{String,1}:
     "b2416"
     "b2415"

Note that **min_growth** represents the actual growth, not a percentage of the **wild-type** growth. For min growth of 
0.1%, do::

    julia> find_essential_genes(model,0.001 * fba(model).f)
    2-element Array{String,1}:
     "b2416"
     "b2415"


.. _find_exchange_reactions:

Find Exchange Reactions
-----------------------

Return either the indices or names of the exchange reactions in the model::

	find_exchange_reactions(model, [output = "index"])

* ``output`` can be either ``index`` or ``name``

**Example**::

	julia> find_exchange_reactions(model, "name")
	20-element Array{String,1}:
	 "EX_ac_e"    
	 "EX_acald_e" 
	 "EX_akg_e"   
	 "EX_co2_e"   
	 "EX_etoh_e"  
	 "EX_for_e"   
	 "EX_fru_e"   
	 "EX_fum_e"   
	 "EX_glc__D_e"
	 "EX_gln__L_e"
	 "EX_glu__L_e"
	 "EX_h_e"     
	 "EX_h2o_e"   
	 "EX_lac__D_e"
	 "EX_mal__L_e"
	 "EX_nh4_e"   
	 "EX_o2_e"    
	 "EX_pi_e"    
	 "EX_pyr_e"   
	 "EX_succ_e"  



.. _find_reactions_from_gene:

Find Reactions From Gene
------------------------

Return the indices of those reactions affected by a specified ``gene``, specified by name or index, and output can be either "index" or "name"::


	find_reactions_from_gene(model, gene, [output = "index"])

* ``output`` can be either ``index`` or ``name``

**Example**

To find the indices of those reactions affected by gene "s0001"::

	julia> find_reactions_from_gene(model, "s0001", "index")
	5-element Array{Int64,1}:
	  2
	 14
	 58
	 69
	 70

**Example**

Find the names of the reactions affected by gene number 3 in model.genes::

	julia> find_reactions_from_gene(model, 3, "name")
	5-element Array{String,1}:
	 "ACALDt"
	 "CO2t"  
	 "H2Ot"  
	 "NH4t"  
	 "O2t"   


.. _find_reactions_from_metabolite:

Find Reactions From Metabolite
------------------------------

Find those reactions a specified metabolite appears in, the metabolite either specified by its name
as it appears in model.mets, or by its index::


	find_reactions_from_metabolite(model, metabolite, [output = "index"])

* ``output`` can be either ``index`` or ``name`` or ``formula``

**Example** 

To find the indices of those reactions where "cit_c" appears::

	julia> find_reactions_from_metabolite(model, "cit_c")
	2-element Array{Int64,1}:
	  4
	 15


.. _gene_deletion:

Gene Deletion
-------------

Return a GeneDeletion object, which containts the result from the n-th combinatorial deletion of genes.::

	gene_deletion(model, n)

**Example**

To calculate every possible double-deletion of genes::

	julia> sol = gene_deletion(model,2)

	Type: GeneDeletion
	                                 Field :   Size  Fieldtype
	    Reactions => Gene-knockout     r_g :    684 Dict{Any,Any}
	             Reactions => Flow     r_f :    684 Dict{Any,Any}
	         Gene-knockout => Flow     g_f :   9453 Dict{Any,Any}

``sol.r_g`` is a dictionary with ``reactions`` as keys, and ``genes`` as values. An element may contain::

	Any[53,90] => Any[String["b1761","b0729"],String["b1761","b0728"]]

Which means that knocking out either ``["b1761","b0729"]`` or ``["b1761","b0728"]`` will disable reactions ``[53,90]``

``sol.r_g`` is a dictionary with ``reactions`` as keys, and the resulting biomass flow ratio to the wild-type flow after its knockout. 

An element from this dictionary may contain::

	sol.r_f[[53,90]] => 0.9533242530236902

So disabling reactions ``[53,90]`` will reduce the biomass (or whatever the objective reaction is) to 95% of the wild-type

``sol.g_f`` is a dictionary with ``gene-pairs`` as keys and the resulting biomass flow ration to the wild-tyoe flow, after the gene-pairs knockout

An element from this dictionary may contain::

	String["b3956","b3919"] => 0.670362


.. _knockout_genes:

Knockout Genes
--------------

Perform fba after knocking out a set of genes, specified either by indices or their names.

Return a GeneKnockout type with fields ``growth``, ``flux``, ``disabled`` and ``affected``

	knockout_genes(model, gene_index)
	knockout_genes(model, gene_name)

**Example**


### Example

    julia> knockout_genes(model, [model.genes[4], model.genes[16], model.genes[99]])

    Type: GeneKnockout
                            Fields 
                            growth :                      0.374230
                              flux :              Array{Float64,1}
                          disabled :                       [12,71]
                          affected ::
                                                   gene:      reactions:
                                              gene_016 :            [12]
                                              gene_004 :             [3]
                                              gene_099 :            [71]

So knocking out genes 4, 16 and 99, would **touch** reactions 12, 3 and 71, respectively, but only reactions 12 and 71 would be disabled                                              


.. _knockout_reactions:

Knockout Reactions
------------------

An algorithm that attempts to return those genes that can be knocked out to disable a specified set of reactions

The reactions can be specified as a single reaction by name or index, or as an array of multiple names or indices::

	knockout_reactions(model, reaction)
	knockout_reactions(model, reactions)


**Example**

See if reaction 5 can be disabled::

	julia> knockout_reactions(model,5)
	1-element Array{Array,1}:
	 [8; 7]

So knocking out genes 8 and 7 will disable reaction 5

.. _print_reaction_formula:

Print Reaction Formula
----------------------

Print the formula for a specified reaction specified by name or string::

	print_reaction_formul(model, reaction)

**Example**

Formula for reaction 5::

	julia> print_reaction_formula(model,5)
		1.0 acon_C_c + 1.0 h2o_c  -> 1.0 icit_c 

.. _reaction_info:

Reaction Info
-------------

Print out the information for a reaction specified by name or string::

	reaction_info(model, reaction)

**Example**

View reaction information for reaction 5 ::

	 Reaction Name:     Aconitase (half-reaction B, Isocitrate hydro-lyase)
	   Reaction ID:                                                  ACONTb
	   Lower Bound:                                                 -1000.0
	   Upper Bound:                                                  1000.0
	     Subsystem:                                       Citric Acid Cycle
	-----------------------------------------------------------------------
	     Metabolite    Coefficient
	       acon_C_c           -1.0
	          h2o_c           -1.0
	         icit_c            1.0

	----------------------------------
	b3115 || b2296 || b1849
