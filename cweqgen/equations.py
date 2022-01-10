import abc

from copy import deepcopy
from fractions import Fraction

from numpy import pi

import astropy.units as u
from astropy.units.quantity import Quantity
from astropy.constants import G, c

from matplotlib import pyplot as plt

from .definitions import ALLOWED_VARIABLES, EQN_DEFINITIONS


def equations(equation, **kwargs):
    """
    This function generates a class holding a requested equation.
    """

    if equation not in EQN_DEFINITIONS:
        raise KeyError(f"Equation '{equation}' is not currently defined")

    # set values for required equation
    eqinfo = EQN_DEFINITIONS[equation]
    kwargs["equation"] = equation
    kwargs["default_fiducial_values"] = eqinfo["default_fiducial_values"]

    kwargs["constants"] = eqinfo.get("constants", [])
    kwargs["additional_values"] = eqinfo.get("additional_values", [])
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
        kwargs["rhs_latex_strings"][key] = ALLOWED_VARIABLES[key]["latex_string"] if key in ALLOWED_VARIABLES else key

    class _EquationBase:
        __metaclass__ = abc.ABCMeta

        def __init__(self, **kwargs):
            self.kwargs = deepcopy(kwargs)  # store copy of initial kwargs

            self.equation_name = kwargs.pop("equation")

            # dictionary to hold default fiducial values
            self.default_fiducial_values = kwargs.pop("default_fiducial_values")

            # list containing additional keyword values that can be used
            self.additional_values = kwargs.pop("additional_values", [])

            # dictionary of functions to convert from additional values into
            # required values (keyed on the required values)
            self.converters = kwargs.pop("converters", {})

            # a list of tuples containing equation constants
            self.constants = kwargs.pop("constants", [])

            self.latex_name = kwargs.pop("latex_string")  # lhs of equation
            self.description = kwargs.pop("description")

            # LaTeX strings for RHS of equation
            self.rhs_latex_strings = kwargs.pop("rhs_latex_strings")

            self.parse_kwargs(**kwargs)

            # get simple format reference
            self.reference_string = kwargs.pop("reference_string", None)

            # equation number in reference
            self.reference_eqno = kwargs.pop("reference_eqno", None)

            # URL of reference in ADS
            self.reference_adsurl = kwargs.pop("adsurl", None)

            # BibTeX for reference (from ADS)
            self.reference_bibtex = kwargs.pop("bibtex", None)

        def parse_kwargs(self, **kwargs):
            """
            Get the required values for the specific equation.
            """

            self.values = {}

            for key in (
                list(self.default_fiducial_values.keys()) + self.additional_values
            ):
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
                            if (
                                not isinstance(value, Quantity)
                                and ALLOWED_VARIABLES[key]["units"] is not None
                            ):
                                value *= u.Unit(ALLOWED_VARIABLES[key]["units"])
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
                                if not eval(str(value.value) + ALLOWED_VARIABLES[key]["sign"]):
                                    raise ValueError(
                                        f"{ALLOWED_VARIABLES[key]['description']} does not have the correct sign"
                                    )

                        self.values[key] = value
                        break

            # perform conversions if required
            for key in self.converters:
                if key not in self.values:
                    try:
                        self.values[key] = self.converters[key](**self.values)
                    except ValueError:
                        pass

        def calculate_constant(self):
            """
            Calculate and return the constant coefficient factor in the equation in
            SI units (if it has dimensions).
            """

            constant = 1.0

            for const in self.constants:
                constant *= eval(str(const[0])) ** Fraction(const[1])

            if isinstance(constant, Quantity):
                return constant.si.decompose()
            else:
                return constant

        def calculate_fiducial(self):
            """
            Calculate and return the product of the fiducial components of the
            equation (excluding the constants) in SI units.
            """

            fiducial = 1.0

            for key in self.default_fiducial_values:
                if key in self.values:
                    # use provided value
                    val = self.values[key]
                else:
                    val = self.default_fiducial_values[key][0]

                if isinstance(val, Quantity):
                    val = val.si

                # get exponent
                exp = Fraction(self.default_fiducial_values[key][1])

                fiducial *= abs(val) ** exp

            if isinstance(fiducial, Quantity):
                return fiducial.si.decompose()
            else:
                return fiducial

        def equation(self, displaytype="string"):
            """
            Generate the LaTeX string for the equation.

            Parameters
            ----------
            displaytype: string
                By default this will return a string containing the LaTeX equation
                text (without bounding "$" symbols). If using a Jupyter Notebook
                this will show as a formatted LaTeX equation. Alternatively, set to
                "matplotlib" to have the output returned as a Matplotlib figure
                object containing the equation.
            """

            latex_equation = self.latex_name + " = "

            constnumstr = ""  # constant numerator
            constdenstr = ""  # constant denominator

            constfractions = {}

            for const in self.constants:
                # value
                if const[0] == "pi":
                    cv = r"\pi"
                else:
                    cv = const[0]

                # exponent
                exp = abs(Fraction(const[1]))

                if str(exp) not in constfractions:
                    constfractions[str(exp)] = ["", ""]

                if Fraction(const[1]) < 0:
                    # denominator value
                    constfractions[str(exp)][1] += " " + cv
                else:
                    # numerator value
                    constfractions[str(exp)][0] += " " + cv

            # construct constant latex string
            conststr = ""
            for exp in constfractions:
                if len(constfractions[exp][0]) == 0:
                    constnumstr = "1"
                else:
                    constnumstr = constfractions[exp][0]

                constdenstr = constfractions[exp][1]

                lbrace, rbrace = (
                    ("", "") if Fraction(exp) == 1 else (r"\left(", r"\right)")
                )
                expstr = "" if exp == "1" else (f"^{{{exp}}}")

                if len(constdenstr) == 0:
                    conststr += f"{constnumstr}{expstr}"
                else:
                    conststr += rf"{lbrace}\frac{{{constnumstr}}}{{{constdenstr}}}{rbrace}{expstr}"

            latex_equation += conststr

            # use default fiducial value keys for the rest of the equation
            varnumstr = ""  # variables numerator
            vardenstr = ""  # variables denominator

            varfractions = {}

            for key in self.default_fiducial_values:
                varlatex = self.rhs_latex_strings[key]

                # get exponent
                exp = abs(Fraction(self.default_fiducial_values[key][1]))

                if str(exp) not in varfractions:
                    varfractions[str(exp)] = ["", ""]

                if Fraction(self.default_fiducial_values[key][1]) < 0:
                    # denominator value
                    varfractions[str(exp)][1] += " " + varlatex
                else:
                    # numerator value
                    varfractions[str(exp)][0] += " " + varlatex

            # construct variables latex string
            varstr = ""
            for exp in varfractions:
                if len(varfractions[exp][0]) == 0:
                    varnumstr = "1"
                else:
                    varnumstr = varfractions[exp][0]

                vardenstr = varfractions[exp][1]

                lbrace, rbrace = (
                    ("", "")
                    if Fraction(exp) == 1 or varnumstr.split() == 1
                    else (r"\left(", r"\right)")
                )
                expstr = "" if Fraction(exp) == 1 else (f"^{{{exp}}}")

                if len(vardenstr) == 0:
                    varstr += f"{lbrace}{{{varnumstr}}}{rbrace}{expstr}"
                else:
                    varstr += (
                        rf"{lbrace}\frac{{{varnumstr}}}{{{vardenstr}}}{rbrace}{expstr}"
                    )

            latex_equation += varstr

            if displaytype.lower() == "matplotlib":
                return EquationLaTeXToImage(latex_equation)
            else:
                return EquationLaTeXString(latex_equation)

        def __repr__(self):
            return str(self.equation())

        def _repr_latex_(self):
            return "$" + str(self.equation()) + "$"

        def fiducial_equation(self, dp=2, brackets="()", displaytype="string"):
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
            coeff = self.evaluate()

            if not isinstance(coeff, Quantity):
                # convert into Quantity
                coeff = Quantity(coeff)

            latex_equation += coeff.to_string(
                precision=(dp + 1), format="latex"
            ).replace("$", "")

            if brackets not in ["()", "{}", "[]", None]:
                raise ValueError(f"Bracket type {brackets} is not recognised")

            lbrace = r"\left" + brackets[0] if brackets is not None else ""
            rbrace = r"\right" + brackets[1] if brackets is not None else ""

            fiducial = ""

            for key in self.default_fiducial_values:
                varlatex = self.rhs_latex_strings[key]

                # get exponent
                # exp = self.default_fiducial_values[key][1]
                exp = Fraction(self.default_fiducial_values[key][1])

                # get value
                if key in self.values:
                    # use provided value
                    val = Quantity(self.values[key])
                else:
                    val = Quantity(self.default_fiducial_values[key][0])

                if exp < 0:
                    numerator = val.to_string(
                        precision=(dp + 1), format="latex"
                    ).replace("$", "")
                    denominator = varlatex
                else:
                    denominator = val.to_string(
                        precision=(dp + 1), format="latex"
                    ).replace("$", "")
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
                return EquationLaTeXToImage(latex_equation)
            else:
                return EquationLaTeXString(latex_equation)

        def evaluate(self):
            """
            Evaluate the equation using the given values.
            """

            const = self.calculate_constant()
            fid = self.calculate_fiducial()

            value = const * fid

            if isinstance(value, Quantity):
                return value.si.decompose()
            else:
                return value

        def rearrange(self, newval, fidval=None):
            """
            Rearrange equation so that the new values is on the left hand side.

            Parameters
            ----------
            newval: str
                The variable should be rearranged to the LHS of the equation.
            fidval: float, Quantity
                The value to use for the original LHS value. If not given this
                will be based on the original fiducial values.
            """

            if newval not in self.default_fiducial_values:
                raise KeyError(f"{newval} is not allowed")

            if fidval is None:
                # use current fiducial values to get value of parameter being
                # swapped from LHS
                curval = self.evaluate()
            else:
                curval = fidval

            exp = Fraction(self.default_fiducial_values[newval][1])

            # check whether values need inverting
            invert = -1 if exp > 0 else 1

            # check whether exponents need changing
            flipfrac = 1 / abs(exp) if abs(exp) != 1 else 1

            # set new fiducial values
            newkwargs = deepcopy(self.kwargs)

            # set new constants
            newconstants = newkwargs.pop("constants")
            for i, c in enumerate(newconstants):
                exp = Fraction(c[1]) * invert * flipfrac
                newconstants[i] = (c[0], str(exp))

            newfiducial = newkwargs.pop("default_fiducial_values")
            newfiducial.pop(newval)

            for val in newfiducial:
                exp = Fraction(newfiducial[val][1]) * invert * flipfrac
                newfiducial[val] = (newfiducial[val][0], str(exp))

            # add in current LHS value
            newfiducial[self.equation_name] = (curval, str(-1 * invert * flipfrac))

            newkwargs["default_fiducial_values"] = newfiducial
            newkwargs["constants"] = newconstants

            newkwargs["latex_string"] = ALLOWED_VARIABLES[newval]["latex_string"]
            newkwargs["description"] = ALLOWED_VARIABLES[newval]["description"]

            newkwargs["rhs_latex_strings"].pop(newval)
            newkwargs["rhs_latex_strings"][self.equation_name] = self.latex_name

            # create new _EquationBase with updated constants and default fiducial values
            return _EquationBase(**newkwargs)

        def __str__(self):
            return self.equation()

    # update the class docstring
    try:
        _EquationBase.__doc__ = eqinfo["docstring"]
    except KeyError:
        pass

    # return the equation class
    return _EquationBase(**kwargs)


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

        self.text = latexstring

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

        self.text = latexstring

        self.dpi = dpi
        self.fig, self.ax = plt.subplots(1)

        try:
            t = self.ax.text(
                0.05,
                0.5,
                "$" + self.text + "$",
                usetex=True,
            )
        except RuntimeError:
            t = self.ax.text(
                0.0,
                0.5,
                "$" + self.text + "$",
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

    def savefig(self, **kwargs):
        kwargs.setdefault("dpi", self.dpi)
        return self.fig.savefig(**kwargs)
