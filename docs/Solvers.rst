Solvers
=======

* :ref:`Supported Solvers`
* :ref:`Methods`
	* :ref:`change_objective`
	* :ref:`get_variable_bounds`
	* :ref:`get_constraint_bounds`
	* :ref:`get_constraint_matrix`
	* :ref:`get_direction`
	* :ref:`get_objective_value`
	* :ref:`get_reduced_costs`
	* :ref:`get_solution`
	* :ref:`get_solution_status`
	* :ref:`get_solution_status_code`
	* :ref:`get_slack`
	* :ref:`get_objective`
	* :ref:`set_direction`
	* :ref:`set_col_bounds`
	* :ref:`set_objective`
	* :ref:`solve`

.. _Supported Solvers:

Supported Solvers
-----------------

Currently, Cobra has support for

* GLPK
* CPLEX
* Gurobi (**Does not support parallel functions**)



All methods that have an optional arguement ``solver`` will accept GLPK, Gurobi or CPLEX, but GLPK is the default

.. _Methods:

Methods
-------

In the descriptions for all methods, LPProb stands for one of the following:

* ``GLPK.Prob``
* ``CPLEX.Model``
* ``Gurobi.Model``
* ``ClpModel``

.. _change_objective:

change objective
^^^^^^^^^^^^^^^^

Change the objective value for a single/multiple reactions.

If not specified, the value/s for the chosen reaction/s will be set to 1.0

Usage::

	change_objective(lp::LPPRob, objective::Number, value::Number = 1.0)
	change_objective(lp::LPPRob, objective::Array{Number}, value::Array=[])

**Example**
Create a lp for e_coli core, where reaction 13 is the biomass reaction. And change
the lp so that the objective reaction is reaction 11 instead of 13::

	lp = setup_lp(model)
	change_objective(lp, 11, 1.0)

**Note** that now the objective for reaction 13 is 0.0. 

To **add** an objective, see **set_objective()**


.. _get_variable_bounds:

get_variable_bounds
^^^^^^^^^^^^^^^^^^^

Return the lower and upper bounds for all reactions or a single reaction from a LPProb::

	lb, ub = get_variable_bounds(lp::LPProb)
	lb, ub = get_variable_bounds(lp::LPProb, index::Number)

**Example**
Get the bonds for reaction 13, the biomass function (in e_coli_core)::

	julia> lb, ub = get_variable_bounds(lp, 13)
	(0.0,1000.0)

.. _get_constraint_bounds:

get_constraint_bounds
^^^^^^^^^^^^^^^^^^^^^

Return the variable bounds for all reactions or a single reaction from a LPProb::

	b = get_constraint_bounds(lp::LPProb)

.. _get_constraint_matrix:

get_constraint_matrix
^^^^^^^^^^^^^^^^^^^^^

Return the coefficient matrix, S from the lp::

	S = get_constraint_matrix(lp::LPProb)

.. _get_direction:

get_direction
^^^^^^^^^^^^^

Return the optimization direction for a LPProb::

	direction = get_direction(lp::LPProb)

Returns either ``maximize`` or ``minimize``

.. _get_objective_value:

get_objective_value
^^^^^^^^^^^^^^^^^^^

Return the objective value for a **solved** LPProb::

	objective = get_objective_value(lp::LPProb)

**Example**

Lets get the biomass growth for a e_coli_core model::

	julia> lp = setup_lp(model)
			GLPK Model 
			   sense           :   maximize 
			   variables       =         95 
			   constraints     =         72 
			   coefficients    =        360 

		solve(lp);
		get_objective_value(lp)
		0.8739215069684305


.. _get_reduced_costs:

get_reduced_costs
^^^^^^^^^^^^^^^^^

Return the reduced costs for a **solved** LPProb::

	reduced_costs = get_reduced_costs(lp::LPProb)


.. _get_solution:

get_solution
^^^^^^^^^^^^

Return the entire solution vector for a **solved** LPProb::

	solution = get_solution(lp::LPPro)

.. _get_solution_status:

get_solution_status
^^^^^^^^^^^^^^^^^^^

Return the solution status. The status can be one of the following

#. Optimal
#. No feasible solution
#. Infeasible
#. Feasible
#. Unbound
#. Undefined

**Example**

Lets see the solution status code for the ecoli core model::

	julia> lp = setup_lp(model);

	solve(lp);

	get_objective_value(lp)
	0.8739215069684305

	get_solution_status(lp)
	"Optimal"


.. _get_solution_status_code:

get_solution_status_code
^^^^^^^^^^^^^^^^^^^^^^^^

Return solver-specific status code of a **solved** problem

Different solvers return different status codes. Fx for GLPK, the status code 5 means 
that the problem has an optimal solution, while for CPLEX, 1 is returned if a problem has an optimal solution::

	status_code = get_solution_status_code(lp::LPProb)

.. _get_slack:

get_slack
^^^^^^^^^

Return the slack for a **solved** LPProb::

	slack = get_slack(lp::LPProb)

.. _get_objective:

get_objective
^^^^^^^^^^^^^

Returns the objective vector from a LPProb, the same vector as is in model.c::

	objective = get_objective(lp::LPProb)


.. _set_direction:

set_direction
^^^^^^^^^^^^^

Set the optimization direction for a LPProb::

	set_direction(lp::GLPK.Prob, direction::AbstractString)

``direction`` can only be ``max`` or ``min``

.. _set_col_bounds:

set_col_bounds
^^^^^^^^^^^^^^

Set lower and upper bounds for a reaction::

	set_col_bounds(lp::GLPK.Prob, index::Number, lb::Number, ub::Number)
	set_col_bounds(lp::GLPK.Prob, index::Number, fixed_value::Number)

The second method is used to fix a reaction to a particular value, commonly used to disable a reaction by fixing in to 0.0. 

.. _set_objective:

set_objective
^^^^^^^^^^^^^

Add another objective to the LPPRob::

	set_objective(lp::GLPK.Prob, objective::Number, value::Number = 1.0)


.. _solve:

solve
^^^^^

Optimize the LPPRob::

	julia> lp = setup_lp(model)
		GLPK Model 
		   sense           :   maximize 
		   variables       =         95 
		   constraints     =         72 
		   coefficients    =        360 

	solve(lp)
	0
