===============
Getting Started
===============

The CBM package is intended for Constraint Based Analysis. 
The package includes methods to read and write models, methods for simulation such as flux balanace analysis and gene deletion and methods
to modify the model, such as changing reaction bounds, adding reactions or reducing the model

* :ref:`installation`
	* :ref:`windows`
	* :ref:`linux`

* :ref:`cobrasetup`

* :ref:`firstrun`

.. _installation:

Installation
============
Firsty, you'll need to install `Julia <http://julialang.org/downloads/platform.html>`_ (version 0.5 or above).

In addition to Julia, you will also want to install a **optimization solver**.  Free solvers are available and CBM supports both GLPK and Clp.
Supported proprietery solvers are CPLEX and Gurobi

.. _windows:

Windows
-------


.. _linux:

Linux
-----

Ubuntu 

To install Julia, run::

	sudo add-apt-repository ppa:staticfloat/juliareleases
	sudo add-apt-repository ppa:staticfloat/julia-deps
	sudo apt-get update
	sudo apt-get install julia

In addition you may need to run:: 

	sudo apt-get install hdf5-tools

.. _cobrasetup:

CBM Setup
===========

Wheter you are on Windows, Linux or MacOS, once you're in the **Julia REPL** your interface will look something similar to this::

	david@david ~ $ julia
	               _
	   _       _ _(_)_     |  A fresh approach to technical computing
	  (_)     | (_) (_)    |  Documentation: http://docs.julialang.org
	   _ _   _| |_  __ _   |  Type "?help" for help.
	  | | | | | | |/ _` |  |
	  | | |_| | | | (_| |  |  Version 0.5.0 (2016-09-19 18:14 UTC)
	 _/ |\__'_|_|_|\__'_|  |  
	|__/                   |  x86_64-linux-gnu

	julia> 

There are a few packages necessary for CBM to function properly. These can all be installed easily through the Julia REPL using the function ``Pkg.add("package name")``, fx ::

	julia> Pkg.add("JSON")


You will need at least **one** of the following solvers.

* GLPK
* CPLEX
* Gurobi
* Clp

And these packages

* JSON
* ProgressMeter
* JuMP
* MAT
* HDF5
* Combinatorics

And then do::

	julia> Pkg.clone("https://github.com/isebarn/CBM")


When CBM finishes installation, call::

	julia> using CBM
	
Now the toolbox is setup and ready for use!

.. _firstrun:

First Run
=========

You can load the example model, ``e_coli_core`` which comes with the package by calling::

	using CBM
	model = load_json(Pkg.dir() * "/CBM/Models/e_coli_core.json")

So your screen should read something like::

	               _
	   _       _ _(_)_     |  A fresh approach to technical computing
	  (_)     | (_) (_)    |  Documentation: http://docs.julialang.org
	   _ _   _| |_  __ _   |  Type "?help" for help.
	  | | | | | | |/ _` |  |
	  | | |_| | | | (_| |  |  Version 0.5.0 (2016-09-19 18:14 UTC)
	 _/ |\__'_|_|_|\__'_|  |  
	|__/                   |  x86_64-linux-gnu

	julia> using CBM
	julia> model = load_json(Pkg.dir() * "/CBM/Models/e_coli_core.json")

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

and to run flux balance analysis::

	julia> fba(model)
	LPSolution: Optimal
	      objective::  0.873922
	           flux::  95 element array
	          slack::  72 element array
	         rcosts::  95 element array



