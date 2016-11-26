============
Modification
============

CBM offers many methods intended for model modification

* :ref:`add_metabolite`
* :ref:`add_reaction`
* :ref:`change_objective`
* :ref:`change_reaction_bounds`
* :ref:`convert_to_irreversible`
* :ref:`find_reversible`
* :ref:`convert_to_reversible`
* :ref:`reduce_model`
* :ref:`remove_gene`
* :ref:`remove_reaction`


.. _add_metabolite:

Adding a Metabolite
-------------------

Add a metabolite to the network.::

	function add_metabolite!(model, metabolite, b = 0, csense = "=")

* ``b`` stands for the **accumulation** or **release** of the metabolite, where the value 0.0 stands for conservation.
* ``csense`` relates to the chosen value of ``b``, only ``<=``, ``=`` and ``>=`` are allowed values

To add a metabolite simply call::

	add_metabolite!(model,"new_metabolite")

Note that adding a metabolite wont have any *real effect* unless the metabolite appears in a reaction.


.. _add_reaction:

Adding a Reaction
-----------------

Add a new reaction to the model::

	add_reaction!(model, reaction_name, reaction; [, lb = -1000.0, ub = 1000.0])
	new_model = add_reaction(model, reaction_name, reaction; [, lb = -1000.0, ub = 1000.0])

To add a reaction you must supply a valid **reaction** formula.

Optionally, choose

* A name for the reaction
* lower- and upper-bounds

Those metabolites that appear in the reaction formula that arent already present in the model will be added. The method will print out the names of **new** metabolites, so if the name of a metabolite gets printed out that was already present in the model, see if you wrote the correct metabolite id as it appears in ``model.mets``

Create a new reaction by typing::

	julia> new_model = add_reaction(model, "new_reaction", "met_A => met_B")
	Note: met_B is a new metabolite
	Note: met_A is a new metabolite
	           rxns :     96 Array{Any,1}
	           mets :     74 Array{Any,1}
	          genes :    137 Array{Any,1}
	              S :   7104 SparseMatrixCSC{Float64,Int64}
	             lb :     96 Array{Float64,1}
	             ub :     96 Array{Float64,1}
	              c :     96 Array{Float64,1}
	              b :     74 Array{Float64,1}
	         csense :      3 Dict{String,Array{Any,1}}
	   rxn_gene_mat :  13015 SparseVector{Float64,Int64}
	       rxn_name :     96 Array{Any,1}
	      rxn_rules :     96 Array{Any,1}
	  rxn_subsystem :     96 Array{Any,1}
	      rxn_extra :      6 Dict{Any,Any}
	    met_formula :     74 Array{Any,1}
	       met_name :     74 Array{Any,1}
	      met_extra :      5 Dict{Any,Any}
	      gene_name :      0 Array{Any,1}
	     gene_extra :      0 Dict{Any,Any}
	    description :     16     String

**Note** that the difference between ``add_reaction!()`` and ``add_reaction()`` is that the ``!`` means that the original ``model`` will be altered, but the ``add_reaction()`` method will return a copy of the original ``model`` 

.. _change_objective:

Change the Objective Reaction
-----------------------------

Change the objective of the model::

	change_objective!(model, reaction; [objective = 1.0])
	new_model = change_objective(model, reaction; [objective = 1.0])

Here, ``reaction`` and ``objective`` may be 

* The name of a reaction and an specified objective value for that reaction
* The names of reaction and specified objective values for these reactions 
* The index of a reaction and a specified objective value for that reaction
* The indices of reaction and specified objective values for these reactions 

**Note** that if the ``objective`` is left blank, if will default to ``1.0``


When performing analysis functions like ``fba()``, usually the default objective reaction will be *biomass growth*

To find a models objective reaction, call::

	objective = find(model.c)

	model.rxns[objective]

	1-element Array{AbstractString,1}:
	 "BIOMASS_Ecoli_core_w_GAM"

So calling ``fba(model,"max")`` will result in maximizing
the flow of the reaction *BIOMASS_Ecoli_core_w_GAM*::

	println(fba(model,"max").f)
	0.8739215069684305


Use either of these to change the objective by supplying the index of the reaction. You can also supply multiple indices if there are **multiple** reactions whose flow you'd like to maximize, this will maximize their **sum of flows**::

	change_objectives!(model, reaction_index; objective = 1.0)
	new_model = change_objective(model, reaction_index; objective = 1.0)

**Note**: Keep in mind that methods with names ending with ``!``, like ``change_objective!()`` **alter** the original model, but those without will return a copy of the model with an altered objective

Try calling::

	change_objective!(model,5)

	println(fba(model,"max").f)
	[20.000000000000004]

The value changes as we are now maximizing reaction 5 in ``model.rxns``, **ACONTb**

.. _change_reaction_bounds:

Change Reaction Bounds
----------------------

Alter the lower- and upper-bounds of a reaction::

	change_reaction_bounds!(model, reaction_index, value, bound = "both")
	new_model = change_reaction_bounds(model, reaction_index, value, bound = "both")

where ``bound`` should be chosen as ``l``, ``u`` or ``both`` to set either the lower bound, upper bound or both bounds to ``value``  
	
or:: 

	change_reaction_bounds(model, reaction_index, lb, ub, bound = "both")

To set both bounds

When performing analysis functions like `fba()`, the result will depend heavily on reaction **upper/lower bounds**

See as an example the result of ``Flux Balance Analysis`` focused on biomass growth (reaction 13 in ecoli_core)::

	println(fba(model,"max").f)
	0.8739215069684305

If we were to call::

	change_reaction_bounds!(model, 4, 0.0, 3.0)

The lower bound would be changed to 0.0, and the upper bound to 3.0 and performing **fba** would yield 

	println(fba(model,"max").f)
	0.8518912912634845

A slightly lower value


To change by reaction index, call either::

	change_reaction_bounds!(model, reaction_index; lb = 0.0, ub = 0.0)
	new_model = change_reaction_bounds(model, reaction_index; lb = 0.0, ub = 0.0)

To change by reaction name, call either::

	change_reaction_bounds!(model, reaction; lb = 0.0, ub = 0.0)
	new_model = change_reaction_bounds(model, reaction; lb = 0.0, ub = 0.0)

To change only either *lower* or *upper* but not both, you can call::

	change_reaction_bounds!(model, reaction_index; lb = 0.0, ub = 0.0, bound   = "both")
	new_model = change_reaction_bounds(model, reaction_index; lb = 0.0, ub = 0.0, bound   = "both")

This method has the **optional arguement** ``bound`` which may be set as **both** or as **l** for **lower** or **u** for **upper**, so::

	change_reaction_bounds!(model, 5, ub = 20.0, "u")

Will change reaction 5 upper bound to 20.0

**Note**: Keep in mind that methods with names ending with ``!``, like ``change_reaction_bounds!()`` **alter** the original model, but those without will return a copy of the model with an altered lower and upper bounds


.. _convert_to_irreversible:

Convert To Irreversible
-----------------------

Convert the model to its irreversible form, such that no reaction is **bi-directional**::

	convert_to_irreversible!(model)
	new_model = convert_to_irreversible(model)

This simple method simply finds the reactions in the model who have **negative lower bounds** and **positive upper bounds** and splits them into **two new reactions**

1. One reaction with a negative lower bound and upper bound fixed at 0.0
2. One reaction with a positive upper bound and lower bound fixed at 0.0


.. _find_reversible:

Find Reversible Reactions
-------------------------

Returns the indices of all reversible reactions::
	
	find_reversible(model)


.. _convert_to_reversible:

Convert To Reversible
---------------------

Convert the model so that duplicate reactions, but with **opposite** direction (or flow), will be joined and bounds merged

	convert_to_reversible!(model)
	new_model = convert_to_reversible(model)

This simple method finds the reactions that have only both **negative** and **positive** bounds and are the "same" reaction, and merges the lower lower-bound and the higher upper-bound


.. _reduce_model:

Reduce Model 
------------

Remove from the model those reaction that dont affect the organism's objective reaction (biomass growth), i.e those reactions that 
exhibit zero flux::

	new_model = reduce_model(model)
	reduce_model!(model)

.. _remove_gene:

Remove Gene
-----------

Remove a gene from the model::

	remove_gene!(model gene)
	new_model = remove_gene(model gene)

Where ``gene`` may be 

* the string name of the gene, such as ``"s0001"``
* its index
* an array of gene names, eg ``["s0001", "b0019"]``
* an array of indices, [1,5,6]

This method is intended to emulate the effect of knocking out genes and thus (possibly) deactivating the reactions that depend on them

Lets try running::

	println(fba(model,"max").f)
	0.8739215069684305

And then we try::

	remove_gene(model,["b0118","b1276"])
	fba(model).f
	1-element Array{Float64,1}:
	 -6.38977e-29	

So these two genes appear to be essential.

.. _remove_reaction:

Remove Reaction
---------------

Remove one or more reactions from the model::

	remove_reaction!(model, reaction)
	new_model = remove_reaction(model, reaction)

Where ``reaction`` may be:

* the name of a single reaction as they appear in model.rxns
* the index of a single reaction
* an array of reaction names
* an array of reaction indices

Lets try running::

	println(fba(model,"max").f)
	0.8739215069684305

And then we try::

	remove_reaction(model,8)
	fba(model).f
	1-element Array{Float64,1}:
	 0.858307

So removing reaction 8, **AKGDH**, makes a tiny difference in the flow of the biomass growth reaction

