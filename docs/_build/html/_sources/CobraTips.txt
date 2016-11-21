Cobra Tips
==========

This section is intended to familiarize the user with the smaller tools and features of the cobra toolbox

Function Naming Conventions
-------------------------

Functions are declared as either::

    method(arguements)
    method!(arguements)

The "!" declares that the method **modifies** some of the arguements

**Example**::

    remove_gene!(model, "s0001")

Will remove the gene "s0001" from the model, but::

    newmodel = remove_gene(model, "s0001")

will return a new copy of the model, with the gene "s0001" removed





