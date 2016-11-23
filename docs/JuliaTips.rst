==========
Julia Tips
==========

This section is intended to familiarize the user with the tools and features of Julia



* :ref:`packages`
* :ref:`funcnames`
* :ref:`funcsignature`
* :ref:`help`
* :ref:`methods`

.. _packages:

Packages
========

Usefull commands for packages are :

* ``Pkg.add(package_name)`` to add a package from name 
* ``Pkg.clone(package_url)`` to add a package from URL, fx from Github
* ``Pkg.rm(package_name)`` to uninstall a package 
* ``Pkg.build(package_name)`` to build a package, this may sometimes be necessary 
* ``Pkg.dir()`` to get the ``path`` to all installed packages 
* ``using "package name"`` to import a package into the enviroment, fx ``Using Cobra``

For more information, see the `full documentation <http://docs.julialang.org/en/release-0.5/manual/packages/>`_ for packages


.. _funcnames:

Function names
==============

Function names ending with an exclamation point (!) modify their arguments. 

Some functions have both modifying (e.g., sort!) and non-modifying (sort) versions.

For Cobra, this means that calling::

	change_reaction_bounds!(model, 13, 0, 1000)

will change the bounds for reaction 13 in the original model, whereas calling::

	new_model = change_reaction_bounds(model, 13, 0, 1000)

returns a new copy of the model, with bounds changed for reaction 13.

.. _funcsignature:

Function argument signature
===========================

Functions in Julia may have either optional arguments **and** keyword arguments.

This makes it easier to use functions that may have multiple necessary arguments, where some of them have sensible default values.

In the documentation, **optional arguments** have a default value assigned to them, for example, the ``flux balance analysis`` function ``fba()`` has the following description ::

	fba(model; [direction = "max"])

Here, ``"max"`` is the default value given to ``direction``

**Keyword arguments** are those that appear inside brackets ``[ ]`` following a semicolon, ``;``
so ``"direction"`` is a keyword argument

meaning you can equivalently call::

	fba(model)
	fba(model, "max") 
	fba(model, direction = "max") 


.. _help:

Getting help
============

To access the help for a function, write a question mark (?) in front of the function name ::

	help?> change_reaction_bounds
	search: change_reaction_bounds find_exchange_reactions

	  change_reaction_bounds(model, reaction, lb, ub)

	  Return a duplicate of the model with the lower/upper bounds of a reaction changed

	  change_reaction_bounds(model, reaction, value[, bound])
	  
	  change_reaction_bounds(model, reaction_index, value[, bound])

	  Return a duplicate of the model after changing a bound

	     Optional arguements:
	    ––––––––––––––––––––––

	  bound:

	    •    'l' to change the lower bound
	      
	    •    'u' to change the upper bound
	      
	    •    'b' to fix lower and upper to the same value (default)


.. _methods:

Methods
=======

Methods within a module
-----------------------

To view every single method (and global variables) from a module, like ``Cobra`` or ``GLPK`` or ``JSON``, write ``JSON.`` and press the ``tab`` key twice::

	julia> JSON.
	ARRAY_BEGIN        LATIN_E             OBJECT_END          eval
	ARRAY_END          LATIN_F             PLUS_SIGN           json
	AssociativeWrapper LATIN_L             Parser              lower
	BACKSLASH          LATIN_N             RETURN              parsefile
	BACKSPACE          LATIN_R             REVERSE_ESCAPES     prefix
	DECIMAL_POINT      LATIN_S             SEPARATOR           print
	DELIMITER          LATIN_T             SOLIDUS             print_escaped
	DIGIT_NINE         LATIN_U             SPACE               printsp
	DIGIT_ZERO         LATIN_UPPER_A       STRING_DELIM        separator
	ESCAPES            LATIN_UPPER_E       State               set_state
	FORM_FEED          LATIN_UPPER_F       TAB                 start_object
	INDENT             MINUS_SIGN          _print              suffix
	JSONPrimitive      NEWLINE             _writejson
	LATIN_A            NOINDENT            end_object
	LATIN_B            OBJECT_BEGIN        escaped


Method variations
-----------------

Julia is implemented with **multiple dispatch**, which means, that methods can have the same name but take different input arguments

For example, the method ``change_reaction_bounds`` has multiple definitions, two of which are ::

	function change_reaction_bounds(model::Model, reaction::String, lb::Number, ub::Number, bound::String = "both")
	function change_reaction_bounds(model::Model, reaction::Number, lb::Number, ub::Number, bound::String = "both")

You will notice that these definitions are a bit different, namely, that in one definition, ``reaction`` must be a ``String`` and in another it should be a ``Number``. So calling::

	change_reaction_bounds!(model, "ACONTa", 0, 1000)
	change_reaction_bounds!(model, 4, 0, 1000)

will do the same thing to the model, because "ACONTa" is the 4th reaction in model.rxns

There are sometimes a few different definitions and often there are numerous definitions, to view **all** definitions, use the ``methods()`` function::

	julia> methods(change_reaction_bounds)
	# 8 methods for generic function "change_reaction_bounds":
	change_reaction_bounds(model::Cobra.Model, reaction::String, value::Number) 
	change_reaction_bounds(model::Cobra.Model, reaction::String, value::Number, bound::String) 
	change_reaction_bounds(model::Cobra.Model, reaction::String, lb::Number, ub::Number) 
	change_reaction_bounds(model::Cobra.Model, reaction::String, lb::Number, ub::Number, bound::String)
	change_reaction_bounds(model::Cobra.Model, reaction_index::Number, value::Number) 
	change_reaction_bounds(model::Cobra.Model, reaction_index::Number, value::Number, bound::String) 
	change_reaction_bounds(model::Cobra.Model, reaction_index::Number, lb::Number, ub::Number) 
	change_reaction_bounds(model::Cobra.Model, reaction_index::Number, lb::Number, ub::Number, bound::String)


