# CBM.jl

[![Build Status](https://travis-ci.org/isebarn/CBM.svg?branch=master)](https://travis-ci.org/isebarn/CBM)

CBM.jl is a COBRA toolbox written in the [Julia](http://julialang.org/downloads/) programming language.

CBM provides a broad range of functionality found in similar COBRA toolboxes, while being easy to use and providing top of the line performance.

Install the package by calling

```julia
Pkg.clone("https://github.com/isebarn/CBM")
```

within the Julia enviroment.

Full documentation is available [here](http://cbm.readthedocs.io/en/latest/index.html)

## Installation
Firsty, version 0.5 (or above) of Julia is required

In addition to Julia, you will also want to install a optimization solver. Free solvers are available and Cobra supports both GLPK and Clp, along with the proprietery solvers CPLEX and Gurobi

## Getting Started


You can load the example model, ``e_coli_core`` which comes with the package by calling

```julia
using Cobra
model = load_model(Pkg.dir() * "/Cobra/Models/e_coli_core.json")
```

So your screen should read something like

```julia
               _
   _       _ _(_)_     |  A fresh approach to technical computing
  (_)     | (_) (_)    |  Documentation: http://docs.julialang.org
   _ _   _| |_  __ _   |  Type "?help" for help.
  | | | | | | |/ _  |  |
  | | |_| | | | (_| |  |  Version 0.5.0 (2016-09-19 18:14 UTC)
 _/ |\__ _|_|_|\__ _|  |  
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
```

and to run flux balance analysis

```julia
julia> fba(model)
FBAsolution: 
       obj::  0.873922
         v::  95 element-array
     slack::  72 element-array
    rcosts::  95 element-array
   success::  Optimal
      info::  SolverInfo("glpk","GLP_OPT")
```
## Benchmarks 

Running benchmark tests on ``flux variability``, ``single gene deletion`` and ``double gene deletion`` on the °°``iJO1366`` model demonstrates the performance capabilities of ``CBM.jl`` in comparison with ``COBRApy`` and the ``COBRA Toolbox``

|                     Simulation | CBM.jl           | COBRApy      | COBRA Toolbox      | Fast FVA      |
|-------------------------------|-----------------:|-----------------:|------------------------:|-----------------:|
| flux variablity analysis |              4.22 |   11.36          |         13.69             |     3.65           |
| single gene deletion     |              2.14 |              4.33 |                          87 | N/A                | 
 double gene deletion     |  319              |              365 |     >3600                | N/A                |



**Note**
Note: The CBM.jl package is currently under development. Please report all problems to dia2@hi.is.

## Licence

© Copyright 2016, [Davíð Isebarn](https://github.com/isebarn/) & Steinn Guðmundsson. Released under the [MIT licence](https://github.com/helpers/helper-copyright/blob/master/LICENSE).
