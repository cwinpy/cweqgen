import abc

from copy import deepcopy
from fractions import Fraction
from astropy.units.format.base import Base

import numpy as np
from numpy import pi

import astropy.units as u
from astropy.units.quantity import Quantity
from astropy.constants import G, c

from matplotlib import pyplot as plt

from sympy import (
    Add,
    Eq,
    expand_power_base,
    lambdify,
    latex,
    Mul,
    Pow,
    powsimp,
    solve,
    Symbol,
    symbols,
    sympify,
)

from .definitions import ALLOWED_VARIABLES, EQN_DEFINITIONS

#: dictionary of constants
CONSTANTS = {"G": G, "c": c, "pi": pi}


def constfunc(name):
    return CONSTANTS[name]


def equations(equation, **kwargs):
    """
    This function generates a :class:`~cweqgen.equations.EquationBase` class
    holding a requested equation. This should always be used to generate an
    equation rather than using the :class:`~cweqgen.equations.EquationBase`
    class itself.
    """

    if equation.lower() not in EQN_DEFINITIONS:
        raise KeyError(f"Equation '{equation}' is not currently defined")

    # set values for required equation
    eqinfo = EQN_DEFINITIONS[equation.lower()]
    kwargs["equation"] = equation.lower()
    kwargs["equation_variable"] = eqinfo["variable"]
    kwargs["default_fiducial_values"] = eqinfo["default_fiducial_values"]
    kwargs["alternative_variables"] = eqinfo.get("alternative_variables", [])
    kwargs["converters"] = eqinfo.get("converters", {})

    # reference information
    try:
        kwargs["reference_string"] = eqinfo["reference"].get("short", None)
        kwargs["reference_eqno"] = eqinfo["reference"].get("eqno", None)
        kwargs["reference_adsurl"] = eqinfo["reference"].get("adsurl", None)
        kwargs["reference_bibtex"] = eqinfo["reference"].get("bibtex", None)
    except KeyError:
        # no references given
        pass

    kwargs["latex_string"] = eqinfo["latex_string"]
    kwargs["description"] = eqinfo["description"]

    kwargs["rhs_latex_strings"] = {}
    for key in kwargs["default_fiducial_values"]:
        kwargs["rhs_latex_strings"][key] = (
            ALLOWED_VARIABLES[key]["latex_string"] if key in ALLOWED_VARIABLES else key
        )

    # check whether generating an equation from parts of from a chain
    parts = eqinfo.get("parts", None)

    # generate equation
    if parts is not None:
        kwargs["parts"] = parts
        
    else:
        chain = eqinfo["chain"]

        # get start equation, given by the first item in chain
        try:
            eq = equations(chain[0])
        except Exception as e:
            raise RuntimeError(f"Could not generate first equation in a chain: {e}")

        eqother = None
        for link in chain[1:]:
            # split into parts
            linkparts = link.split()

            if len(linkparts) != 2:
                raise ValueError("chain components must contain 2 values")

            if "equals" == linkparts[0].strip().lower():
                # an equality
                eqother = equations(linkparts[1].strip())
            elif "rearrange" == linkparts[0].strip().lower():
                # rearrange
                varname = linkparts[1].strip().lower()

                if eqother is not None:
                    eq = eq.rearrange(varname, equal=eqother)
                    eqother = None
                else:
                    eq = eq.rearrange(varname)
            elif "substitute" == linkparts[0].strip().lower():
                # substitute
                subeqname = linkparts[1].strip().lower()
                subeq = equations(subeqname)

                eq = subeq.substitute(eq)

        kwargs["parts"] = eq.parts

    eq = EquationBase(**kwargs)  

    # update the equation docstring
    try:
        eq.__doc__ = eqinfo["docstring"]
    except KeyError:
        pass

    # return the equation class
    return eq


class EquationBase:
    __metaclass__ = abc.ABCMeta

    def __init__(self, **kwargs):
        """
        Base class for holding equations.
        """

        self.kwargs = deepcopy(kwargs)  # store copy of initial kwargs

        self.equation_name = kwargs.pop("equation")
        self.variable = kwargs.pop("equation_variable")

        # a list of tuples containing parts of equation
        self.parts = kwargs.pop("parts")

        # dictionary to hold default fiducial values
        self.default_fiducial_values = kwargs.pop("default_fiducial_values")

        # list containing additional keyword values that can be used
        self.alternative_variables = kwargs.pop("alternative_variables", [])

        # dictionary of functions to convert from additional values into
        # required values (keyed on the required values)
        self.converters = kwargs.pop("converters", {})

        self.latex_name = kwargs.pop("latex_string")  # lhs of equation
        self.description = kwargs.pop("description")

        # LaTeX strings for RHS of equation
        self.rhs_latex_strings = kwargs.pop("rhs_latex_strings")

        self.values = self.parse_kwargs(**kwargs)

        # get simple format reference
        self.reference_string = kwargs.pop("reference_string", None)

        # equation number in reference
        self.reference_eqno = kwargs.pop("reference_eqno", None)

        # URL of reference in ADS
        self.reference_adsurl = kwargs.pop("reference_adsurl", None)

        # BibTeX for reference (from ADS)
        self.reference_bibtex = kwargs.pop("reference_bibtex", None)

        # make sure Sympy version of equation is generated
        _ = self.sympy_var
        _ = self.constant

    def parse_kwargs(self, **kwargs):
        """
        Get the required values for the specific equation.
        """

        values = {}

        for key in list(self.default_fiducial_values.keys()) + self.alternative_variables:
            value = self.check_for_alias(key, **kwargs)
            if value is not None:
                values[key] = value

        # perform conversions if required
        for key in self.converters:
            if key not in values:
                try:
                    values[key] = self.converters[key](**values)
                except ValueError:
                    pass

        return values

    @staticmethod
    def check_for_alias(key, **kwargs):
        """
        Check whether any alias of the key is given in the keyword arguments.
        """

        # check aliases
        if key in ALLOWED_VARIABLES:
            aliases = ALLOWED_VARIABLES[key]["aliases"]
        else:
            aliases = [key]

        for alias in aliases:
            if alias in kwargs:
                value = kwargs[alias]

                # check value has compatible units
                if key in ALLOWED_VARIABLES:
                    if not isinstance(value, Quantity):
                        if ALLOWED_VARIABLES[key]["units"] is not None:
                            value *= u.Unit(ALLOWED_VARIABLES[key]["units"])
                        else:
                            value = Quantity(value)
                    elif (
                        isinstance(value, Quantity)
                        and ALLOWED_VARIABLES[key]["units"] is not None
                    ):
                        try:
                            _ = value.to(ALLOWED_VARIABLES[key]["units"])
                        except (u.UnitConversionError, ValueError) as e:
                            raise IOError(
                                f"{ALLOWED_VARIABLES[key]['description']} units are not compatible:\n{e}"
                            )

                    # check value has correct sign
                    if ALLOWED_VARIABLES[key]["sign"] is not None:
                        strrep = (
                            "np.array("
                            + np.array2string(np.asarray(value.value), separator=", ")
                            + ")"
                        )

                        if not eval(
                            "np.all(" + strrep + ALLOWED_VARIABLES[key]["sign"] + ")"
                        ):
                            raise ValueError(
                                f"{ALLOWED_VARIABLES[key]['description']} does not have the correct sign"
                            )

                return value

        return None

    def calculate_fiducial(self, **kwargs):
        """
        Calculate and return the product of the fiducial components of the
        equation (excluding the constants) in SI units.

        Keyword arguments can be passed using the keys of the
        :attr:`.EquationBase.default_fiducial_values` giving the values at
        which to perform the calculation. Otherwise the default values, or
        values set at initialisation are used. These can be 1d arrays. If
        the arrays have different lengths, or the "mesh" keyword argument
        is given and is True, then a mesh grid will be created over the
        space and the values returned on that mesh. If arrays of equal
        length are given, and the "mesh" keyword is not given, then values
        will be calculated assuming for each set of equivalently positioned
        values in the array.
        """

        # fiducialunits = 1.0
        funcargs = {}
        arglens = []  # store lengths of arguments

        self._unitdic = {}

        values = self.parse_kwargs(**kwargs)

        for key in self.var_names:
            if key in values:
                # use provided value
                val = values[key]
            else:
                val = self.default_fiducial_values[key]

            if not isinstance(val, Quantity):
                val = Quantity(val)

            funcargs[key] = np.abs(val.si)

            try:
                arglens.append(len(funcargs[key]))
            except TypeError:
                arglens.append(1)

        # check whether multiple values are arrays and create a mesh if
        # necessary
        usemesh = kwargs.get("mesh", False)
        mesh = None
        if sum([length > 1 for length in arglens]) > 1:
            idx = np.argwhere(np.array(arglens) > 1).flatten()

            # create mesh is "mesh" is True or arguments have different lengths
            if usemesh or not np.all(
                (np.array([arglens[i] for i in idx]) - arglens[idx[0]]) == 0
            ):
                keys = [list(funcargs.keys())[i] for i in idx]

                mesh = np.meshgrid(*[funcargs[key] for key in keys])
                for i, key in enumerate(keys):
                    funcargs[key] = mesh[i].flatten()

        # evaluate equation (loop through each part)
        fiducial = 1.0
        for key in funcargs:
            eveq = self._sympy_var_lambda[key](**{key: funcargs[key]})
            # NOTE: can't use *= as it causes a SegFault
            fiducial = fiducial * eveq

        if mesh is not None:
            fiducial = np.reshape(fiducial, mesh[0].shape)

        if isinstance(fiducial, Quantity):
            return fiducial.si.decompose()
        else:
            return fiducial

    def equation(self, displaytype="string", **latexkwargs):
        """
        Generate the LaTeX string for the equation.

        Parameters
        ----------
        displaytype: str
            By default this will return a string containing the LaTeX equation
            text (without bounding "$" symbols). If using a Jupyter Notebook
            this will show as a formatted LaTeX equation. Alternatively, set to
            "matplotlib" to have the output returned as a Matplotlib figure
            object containing the equation.
        latexkwargs: dict
            Keyword parameters that can be passed to the
            :func:`sympy.printing.latex.latex` function. By default the
            ``fold_frac_powers`` option is set to True, the ``root_notation``
            option is set to False and the LaTeX symbol names are defined by
            the equation definition values.
        """

        # use powsimp to put values with common exponents together
        seq_const = powsimp(self.sympy_const, force=True)
        seq_var = powsimp(self.sympy_var, force=True)

        # set defaults
        symrep = {symbols(key): val for key, val in self.rhs_latex_strings.items()}
        symrep[symbols(self.variable)] = self.latex_name

        latexkwargs.setdefault("root_notation", True)
        latexkwargs.setdefault("fold_frac_powers", True)
        latexkwargs.setdefault("long_frac_ratio", 2.0)

        mode = latexkwargs.pop("mode", "plain")
        delim = "$" if displaytype == "matplotlib" or mode == "inline" else ""

        if seq_const != 1:
            latex_equation_const = latex(seq_const, **latexkwargs)
        else:
            latex_equation_const = ""

        latexkwargs["root_notation"] = False
        latexkwargs.setdefault("symbol_names", symrep)
        latex_equation_var = latex(seq_var, **latexkwargs)

        latex_equation = f"{delim}{self.latex_name} = {latex_equation_const}{latex_equation_var}{delim}"

        if displaytype.lower() == "matplotlib":
            return EquationLaTeXToImage(latex_equation)
        else:
            return EquationLaTeXString(latex_equation)

    @property
    def eqn(self):
        """
        The equation as a string.
        """

        return self.equation()

    def __str__(self):
        return str(self.equation())

    def __repr__(self):
        return str(self.equation())

    def _repr_latex_(self):
        return "$" + str(self.equation()) + "$"

    def fiducial_equation(self, dp=2, brackets="()", displaytype="string", **kwargs):
        """
        Generate the LaTeX string for the equation inserting in fiducial values.

        Parameters
        ----------
        dp: int
            The number of decimal places to use for non-integer fiducial values.
        brackets: str
            The style of brackets to use, e.g., "()", "{}" or "[]". Defaults is
            round parentheses "()".
        displaytype: string
            By default this will return a string containing the LaTeX equation
            text (without bounding "$" symbols). If using a Jupyter Notebook
            this will show as a formatted LaTeX equation. Alternatively, set to
            "matplotlib" to have the output returned as a Matplotlib figure
            object containing the equation.
        """

        latex_equation = self.latex_name + " = "

        # add in coefficient
        values = deepcopy(self.values)
        for key, val in kwargs.items():
            values[key] = val
        values = self.parse_kwargs(**values)

        coeff = self.evaluate(**values)

        if not isinstance(coeff, Quantity):
            # convert into Quantity
            coeff = Quantity(coeff)

        latex_equation += coeff.to_string(precision=(dp + 1), format="latex").replace(
            "$", ""
        )

        if brackets not in ["()", "{}", "[]", None]:
            raise ValueError(f"Bracket type {brackets} is not recognised")

        lbrace = r"\left" + brackets[0] if brackets is not None else ""
        rbrace = r"\right" + brackets[1] if brackets is not None else ""

        fiducial = ""

        for key in self.var_names:
            arg = self._sympy_var_parts[key]
            varlatex = self.rhs_latex_strings[key] if not isinstance(arg, Pow) else latex(arg.base, symbol_names={symbols(key): self.rhs_latex_strings[key]})

            # get exponent
            exp = 1 if not isinstance(arg, Pow) else Fraction(*arg.exp.as_numer_denom())

            # get value
            if key in values:
                # use provided value
                val = Quantity(values[key])
            else:
                val = Quantity(self.default_fiducial_values[key])

            if exp != 1:
                # evaluate base for cases such as args being (n + 1)^x
                val = float(arg.base.subs([(symbols(key), val.value)]).evalf()) * val.unit

            if exp < 0:
                numerator = val.to_string(precision=(dp + 1), format="latex").replace(
                    "$", ""
                )
                denominator = varlatex
            else:
                denominator = val.to_string(precision=(dp + 1), format="latex").replace(
                    "$", ""
                )
                numerator = varlatex

            if abs(exp) != 1:
                expstr = "^{" + str(abs(exp)) + "}"
            else:
                expstr = ""

            fiducial += (
                rf"{lbrace}\frac{{{numerator}}}{{{denominator}}}{rbrace}{expstr} "
            )

        latex_equation += fiducial

        if displaytype.lower() == "matplotlib":
            return EquationLaTeXToImage("$" + latex_equation + "$")
        else:
            return EquationLaTeXString(latex_equation)

    def evaluate(self, **kwargs):
        """
        Evaluate the equation using the given values.

        Keyword arguments can be passed using the keys of the
        :attr:`.EquationBase.default_fiducial_values` giving the values at
        which to perform the calculation. Otherwise the default values, or
        values set at initialisation are used. These can be 1d arrays,
        where if more that one value is an array of then a mesh grid will
        be created over the space and the values returned on that mesh.
        """

        const = self.constant
        fid = self.calculate_fiducial(**kwargs)

        value = const * fid

        if isinstance(value, Quantity):
            return value.si.decompose()
        else:
            return value

    def rearrange(self, newvar, fidval=None, equal=None):
        r"""
        Rearrange the equation so that a different variable is on the left hand side,
        i.e., solve the equation for the given variable.

        Parameters
        ----------
        newvar: str
            The variable should be rearranged to the LHS of the equation.
        fidval: float, Quantity
            The value to use for the original LHS value. If not given this
            will be based on the original fiducial values.
        equal: EquationBase
            You can pass another equation to set as equal to the current
            equation before rearranging.

        Returns
        -------
        neweq: EquationBase
            Returns a new equation. The current equation is left unchanged.

        Example
        -------
        For example rearrange the braking index equation to put the
        frequency derivative :math:`\dot{f}` on the left hand side:

        >>> eqn = equations("brakingindex")
        >>> reqn = eqn.rearrange("rotationfdot")

        which gives:

        .. math::

            \dot{f}_{\rm rot} = \frac{1}{n^{1/2}} \left(\ddot{f}_{\rm rot} f_{\rm rot}\right)^{1/2}
        """

        if newvar not in self.var_names:
            raise KeyError(f"{newvar} is not allowed")

        if fidval is None:
            # use current fiducial values to get value of parameter being
            # swapped from LHS
            curval = self.evaluate()
        else:
            curval = fidval

        # check whether equating to other equation
        eq = Eq(equal.sympy.rhs, self.sympy.rhs) if isinstance(equal, EquationBase) else self.sympy

        # rearrange using solve (use the last value in solution in case two solutions
        # from sqrt). Expand out any variables with the same powers, so that they form
        # separate parts of the new equation
        neweq = expand_power_base(solve(eq, symbols(newvar))[-1], force=True)

        newkwargs = deepcopy(self.kwargs)

        # set new equation parts
        newkwargs["parts"] = self.generate_parts(neweq)

        newfiducial = newkwargs.pop("default_fiducial_values")
        newfiducial.pop(newvar)
        newfiducial[self.equation_name] = curval

        if isinstance(equal, EquationBase):
            newfiducial.update(equal.kwargs["default_fiducial_values"])
        newkwargs["default_fiducial_values"] = newfiducial

        newkwargs["latex_string"] = ALLOWED_VARIABLES[newvar]["latex_string"]
        newkwargs["description"] = ALLOWED_VARIABLES[newvar]["description"]

        newkwargs["rhs_latex_strings"].pop(newvar)
        
        if not isinstance(equal, EquationBase):
            newkwargs["rhs_latex_strings"][newvar] = self.latex_name
        else:
            newkwargs["rhs_latex_strings"].update(equal.kwargs["rhs_latex_strings"])

        newkwargs["equation"] = newvar  # reset the equation name
        newkwargs["equation_variable"] = newvar

        # create new EquationBase with updated constants and default fiducial values
        return EquationBase(**newkwargs)

    def substitute(self, other):
        r"""
        Substitute another equation into the current equation.

        Parameters
        ----------
        other: EquationBase
            The other equation to substitute into the current equation.

        Returns
        -------
        neweq: EquationBase
            Returns a new equation. The current equation is left unchanged.

        Example
        -------
        For example put the :math:`h_0` spin-down limit in terms of the
        braking index :math:`n` and second frequency derivative
        :math:`\ddot{f}` (with the help of
        :meth:`~cweqgen.equations.EquationBase.rearrange`):

        >>> # get braking index equation
        >>> eqn = equations("brakingindex")
        >>> # rearrange to put fdot on lhs
        >>> reqn = eqn.rearrange("rotationfdot")
        >>> # get h0 spindown equation
        >>> eqnh0sd = equations("h0spindown")
        >>> # substitute in the rearranged equation
        >>> neweq = eqnh0sd.substitute(reqn)

        with gives:

        .. math::

            h_0^{\rm sd} = \frac{\sqrt{10} \sqrt{G}}{2 c^{3/2}}\frac{I_{zz}^{1/2} \ddot{f}_{\rm rot}^{1/4}}{d \left(n f_{\rm rot}\right)^{1/4}}
        """

        if not isinstance(other, EquationBase):
            raise TypeError("Other equation is not the right type")

        try:
            rhs = self.sympy.rhs
            neweq = rhs.subs([(other.sympy.lhs, other.sympy.rhs)])
        except Exception as e:
            raise RuntimeError(f"Could not perform substitution: {e}")

        newkwargs = deepcopy(self.kwargs)
        newfiducial = newkwargs.pop("default_fiducial_values")
        newfiducial.pop(other.variable)
        newfiducial.update(other.default_fiducial_values)

        newkwargs["default_fiducial_values"] = newfiducial

        newkwargs["rhs_latex_strings"].pop(other.variable)
        newkwargs["rhs_latex_strings"].update(other.rhs_latex_strings)

        newkwargs["parts"] = self.generate_parts(expand_power_base(neweq, force=True))

        return EquationBase(**newkwargs)

    def __str__(self):
        return str(self.equation())

    @property
    def sympy_const(self):
        """
        Construct and return a :class:`sympy.core.mul.Mul` containing the
        constants in the equation.
        """

        if not hasattr(self, "_sympy_const"):
            self._sympy_const = 1
            sympy_const_unit_values = {}

            for arg in self.sympy.rhs.args:
                if arg.is_constant():
                    self._sympy_const *= arg
                else:
                    if isinstance(arg, Symbol):
                        name = arg.name
                    elif isinstance(arg, Pow):
                        name = str(arg.base)
                    else:
                        name = None

                    if name in CONSTANTS:
                        sympy_const_unit_values[name] = CONSTANTS[name]
                        self._sympy_const *= arg

            # evaluate constant by creating a lamdified function
            if self._sympy_const != 1:
                constf = lambdify(
                    [symbols(name) for name in sympy_const_unit_values.keys()],
                    self._sympy_const,
                    modules=["numpy"],
                )
                constant = constf(**sympy_const_unit_values)

                if isinstance(constant, Quantity):
                    self._constant = constant.si.decompose()
                else:
                    self._constant = constant
            else:
                self._constant = self._sympy_const

        return self._sympy_const

    @property
    def constant(self):
        """
        Return the constant coefficient factor in the equation in SI units
        (if it has dimensions).
        """

        if not hasattr(self, "_constant"):
            _ = self.sympy_const

        return self._constant

    @property
    def sympy_var(self):
        """
        Construct and return a :class:`sympy.core.mul.Mul` containing the
        variables in the equation.
        """

        if not hasattr(self, "_sympy_var"):
            self._sympy_var = 1
            self._sympy_var_parts = {}
            self._sympy_var_lambda = {}

            # get variable names
            self._var_names = self._gather_var_argnames(self.sympy.rhs)

            for arg in self.sympy.rhs.args:
                if not arg.is_constant():
                    if isinstance(arg, Symbol):
                        name = arg.name
                    elif isinstance(arg, (Pow, Add)):
                        name = self._gather_var_argnames(arg)
                        if len(name) > 1:
                            raise ValueError("Currently we cannot deal with Pow or Add terms containing multiple variables")
                        elif len(name) == 1:
                            name = name[0]
                        else:
                            name = str(arg.base)
                    else:
                        raise TypeError("")

                    if name not in CONSTANTS:
                        self._sympy_var *= arg

                        # store each part
                        self._sympy_var_parts[name] = arg

                        # create functions for each part of the equation
                        self._sympy_var_lambda[name] = lambdify(
                            [symbols(name)],
                            arg,
                            modules=["numpy"],
                        )

        return self._sympy_var

    @staticmethod
    def _gather_var_argnames(args):
        """
        Find the variable names in an equation.
        """

        argnames = []
        cargs = args.args

        for arg in cargs:
            nextargs = arg.args
            if isinstance(arg, Symbol):
                if arg.name not in CONSTANTS:
                    argnames.append(arg.name)

            if len(nextargs) == 0:
                continue
            else:
                argnames.extend(EquationBase._gather_var_argnames(arg))

        return argnames

    @property
    def var_names(self):
        """
        Get variable names from Sympy representation.
        """

        if not hasattr(self, "_var_names"):
            _ = self.sympy_var

        return self._var_names

    @property
    def sympy_var_func(self):
        if not hasattr(self, "_sympy_var_func"):
            _ = self.sympy_var

        return self._sympy_var_func

    @property
    def sympy(self):
        """
        Construct and return the equation as a
        :class:`sympy.core.relational.Equality`.
        """

        if not hasattr(self, "_sympy"):
            eqstr = " * ".join([f"({c[0]})**({c[1]})" for c in self.parts])
            sympeq = sympify(eqstr)
            self._sympy = Eq(symbols(self.variable), sympeq)

        return self._sympy

    @staticmethod
    def generate_parts(eqn, pow=1):
        """
        Generate a list of "parts" from a Sympy equation.
        """

        parts = []

        for arg in eqn.args:
            if isinstance(arg, Mul):
                parts.extend(EquationBase.generate_parts(arg))
            elif isinstance(arg, Pow):
                if isinstance(arg.base, Pow):
                    exp = pow * arg.exp
                    parts.extend(EquationBase.generate_parts(arg.base, pow=exp))
                else:
                    parts.append((str(arg.base), str(arg.exp)))
            elif isinstance(arg, Symbol):
                if pow != 1 and len(eqn.args) == 2:
                    # it's got to this part from a Pow object
                    parts.append((arg.name, str(pow * eqn.args[-1])))
                    break
                else:
                    parts.append((arg.name, str(pow)))
            elif str(arg) != "1" and str(arg) != "-1":
                parts.append((str(arg), "1"))

        return parts


class EquationLaTeXString:
    def __init__(self, latexstring):
        """
        Class to hold a LaTeX equation string. It has a _repr_latex_ method to
        hook into the Jupyter notebook rich display system
        https://ipython.readthedocs.io/en/stable/config/integrating.html#rich-display
        and show the resulting string as a LaTeX equation.

        Parameters
        ----------
        latexstring: str
            The LaTeX string defining the equation.
        """

        self.text = str(latexstring)

    def __repr__(self):
        return self.text

    def __str__(self):
        return self.text

    def _repr_latex_(self):
        return "$" + self.text + "$"


class EquationLaTeXToImage:
    def __init__(self, latexstring, dpi=200):
        """
        Class to hold a LaTeX equation string and covert it to an image for
        display in a Jupyter notebook.

        Parameters
        ----------
        latexstring: str
            The LaTeX string defining the equation.
        dpi: int
            The resolution (dots per inch) of the output plot.
        """

        self.text = str(latexstring)

        self.dpi = dpi
        self.fig, self.ax = plt.subplots(1)

        try:
            t = self.ax.text(
                0.05,
                0.5,
                self.text,
                usetex=True,
            )
        except RuntimeError:
            t = self.ax.text(
                0.0,
                0.5,
                self.text,
            )
        self.ax.axis("off")

        # crop figure to tight around the text
        bb = t.get_window_extent(renderer=self.fig.canvas.get_renderer())
        transf = self.ax.transData.inverted()
        bb_datacoords = bb.transformed(transf)

        rect = bb_datacoords.get_points().flatten()
        rect[3] += 0.1 * rect[3]  # add 10% on to upper value for some reason!

        tight_params = {"rect": rect}

        self.fig.set_dpi(self.dpi)
        self.fig.set_tight_layout(tight_params)

    def savefig(self, *arg, **kwargs):
        kwargs.setdefault("dpi", self.dpi)
        return self.fig.savefig(*arg, **kwargs)
