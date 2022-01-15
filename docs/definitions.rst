##################
Defining equations
##################


The equation definition information is held in the :obj:`~cweqgen.definitions.EQN_DEFINITIONS`
dictionary. If defining a new equation the information should be added into this dictionary. The
dictionary keys give the equation's "name", with which it must be referred to. The values are
dictionaries containing the following keys:

"description (required)":
   a string giving a brief description of the equation

"variable (required)":
   the variable that the equation is defining (this must be from the names in
   :obj:`~cweqgen.definitions.ALLOWED_VARIABLES`). This may be the same as the equation name, but it
   doesn't have to be. For example, the equation named "h0spindown" defines a "h0" variable.

"latex_string (required)":
   a LaTeX math string (excluding ``$`` symbols) with which to represent the left-hand-side of the
   equation.

"parts (required)":
   a list of 2-tuples containing the all parts of the equation (constants and variables) and their
   exponents. These should both be string values. For variables the parameter names (the first
   value in the tuple) must be consistent with the names in the :obj:`~cweqgen.definitions.ALLOWED_VARIABLES`
   dictionary. Constants can be number value strings or the strings containing "G", "c" or "pi".

"default_fiducial_values (required)":
   a dictionary with keys for all variable parameters in the equation. Each key
   gives contains the default value of that variable (with appropriate
   :class:`astropy.units.Unit`). The variable parameter names must be consistent with the names in
   the :obj:`~cweqgen.definitions.ALLOWED_VARIABLES` dictionary and those given in **parts**.

"alternative_variables":
   a list of variable names that can be used to derive a subset of the variables
   in the "default_fiducial_values", i.e., if the equation requires the "rotationfrequency" then
   this could contain "rotationperiod", which can instead be used to derive the rotation frequency.
   The variable parameter names must be consistent with the names in the
   :obj:`~cweqgen.definitions.ALLOWED_VARIABLES` dictionary.

"converters":
   a dictionary of conversion functions to convert the additional values into the required values.

"reference":
   a dictionary containing the following keys:

   "short":
      a short string version of the reference

   "adsurl":
      the URL of the `NASA ADS <https://ui.adsabs.harvard.edu/>`_ entry for the paper

   "eqno":
      the equation number in the given reference

   "bibtex":
      the NASA ADS `BibTeX <http://www.bibtex.org/>`_ entry for the paper

"docstring":
   a Sphinx-style docstring for the function, e.g., 

   .. literalinclude:: docstringexample.txt


.. automodule:: cweqgen.definitions
   :members:
