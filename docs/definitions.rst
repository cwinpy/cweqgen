##################
Defining equations
##################


The equation definition information is held in the :obj:`~cweqgen.definitions.EQN_DEFINITIONS`
dictionary. This is created from a set of YAML files containing the information for each individual
equation. To define a new equation a new YAML file must be created. The file name (preceding the
YAML extension) is used to define the equation's name with which it must be referred to.

The file should contain the following information:

"description (required)":
   a string giving a brief description of the equation

"variable (required)":
   the variable that the equation is defining (this must be from the names in
   :obj:`~cweqgen.definitions.ALLOWED_VARIABLES`). This may be the same as the equation name, but it
   doesn't have to be. For example, the equation named "h0spindown" defines a "h0" variable.

"latex_string (required)":
   a LaTeX math string (excluding ``$`` symbols) with which to represent the left-hand-side of the
   equation.

"parts (required or chain)":
   pairs of values containing the all parts of the equation (constants and variables) and their
   exponents. These should both be string values. For variables the parameter names (the first
   value in the tuple) must be consistent with the names in the :obj:`~cweqgen.definitions.ALLOWED_VARIABLES`
   dictionary. Constants can be number value strings or the strings containing "G", "c" or "pi".

"chain (required or parts)":
   instead of providing the parts of the equation, it can be constructed from other equations that
   are already defined. To do this the chain of equations, rearrangments and substitutions needs to
   be set. This value should be a list with the first entry being the starting equation name.
   Subsequent entries are strings containing the words "equals", "rearrange", or "substitute"
   followed by an equation or variable name. "equals" should be followed by an equation name which
   will be set as equal to the current equation for the next "rearrangement"; "rearrange" should be
   followed by a variable name, for which the current equation should be solved for; and,
   "substitute" should be followed by the equation name into which the current equation will be
   substituted.

"default_fiducial_values (required)":
   pairs of keys and values all variable parameters in the equation. Each value gives the default
   value of that variable (the :class:`astropy.units.Unit` can be set be giving a unit string wrapped
   in the word "Unit", e.g., Unit("Hz")). The variable parameter names must be consistent with the
   names in the :obj:`~cweqgen.definitions.ALLOWED_VARIABLES` dictionary and those given in **parts**.

"alternative_variables":
   a list of variable names that can be used to derive a subset of the variables
   in the "default_fiducial_values", i.e., if the equation requires the "rotationfrequency" then
   this could contain "rotationperiod", which can instead be used to derive the rotation frequency.
   The variable parameter names must be consistent with the names in the
   :obj:`~cweqgen.definitions.ALLOWED_VARIABLES` dictionary.

"converters":
   a set of pairs of variable names and conversion functions to convert the additional values into
   the required values.

"reference":
   information about the a reference for the equation, containing:

   "short":
      a short string version of the reference

   "adsurl":
      the URL of the `NASA ADS <https://ui.adsabs.harvard.edu/>`_ entry for the paper

   "eqno":
      the equation number in the given reference

   "bibtex":
      the NASA ADS `BibTeX <http://www.bibtex.org/>`_ entry for the paper

An example YAML file for an equation is:

.. literalinclude:: eqndefinitionexample.yaml
   :language: yaml

.. automodule:: cweqgen.definitions
   :members:
