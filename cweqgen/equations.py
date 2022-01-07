import abc
import io

import numpy as np
from numpy import pi

import astropy.units as u
from astropy.units.quantity import Quantity
from astropy.constants import G, c

from matplotlib import pyplot as plt

from .converters import *
from .docstrings import DOCSTRINGS
from .reference import REFERENCES


ALLOWED_VALUES = {
    "h0": {
        "value": "Gravitational wave amplitude",
        "latex_string": r"h_0",
        "aliases": ["h0", "h_0", "$h_0$"],
        "units": None,
        "sign": ">= 0",
    },
    "h0spindown": {
        "value": "Gravitational wave amplitude spin-down limit",
        "latex_string": r"h_0^{\rm sd}",
        "aliases": ["h0spindown", "h0sd", "h_0sd", "h_0^{sd}"],
        "units": None,
        "sign": ">= 0",
    },
    "ellipticity": {
        "value": "Neutron star ellipticity",
        "latex_string": r"\varepsilon",
        "aliases": [
            "ellipticity",
            "ell",
            "eps",
            "epsilon",
            "$\epsilon$",
            "$\varepsilon$",
        ],
        "units": None,
        "sign": ">= 0",
    },
    "ellipticityspindown": {
        "value": "Spin-down limit for neutron star ellipticity",
        "latex_string": r"\varepsilon^{\rm sd}",
        "aliases": ["ellipticityspindown", "ellsd", "epssd", "epsilonsd"],
        "units": None,
        "sign": ">= 0",
    },
    "massquadrupole": {
        "value": "Mass quadrupole moment (l=m=2)",
        "latex_string": r"Q_{22}",
        "aliases": ["massquadrupole", "q22", "$q_{22}"],
        "units": "kg m^2",
        "sign": ">= 0",
    },
    "massquadrupolespindown": {
        "value": "Spin-down limit for the mass quadrupole moment (l=m=2)",
        "latex_string": r"Q_{22}^{\rm sd}",
        "aliases": ["massquadrupolespindown", "q22sd"],
        "units": "kg m^3",
        "sign": ">= 0",
    },
    "rotationfrequency": {
        "value": "Source rotational frequency",
        "latex_string": r"f_{\rm rot}",
        "aliases": [
            "rotationfrequency",
            "frot",
            "spinfrequency",
            "fspin",
            "f0rot",
            "f0spin",
        ],
        "units": "Hz",
        "sign": ">= 0",
    },
    "gwfrequency": {
        "value": "Gravitational-wave frequency",
        "latex_string": r"f_{\rm gw}",
        "aliases": ["gwfrequency", "fgw", "f0gw"],
        "units": "Hz",
        "sign": ">= 0",
    },
    "rotationperiod": {
        "value": "Source rotational period",
        "latex_string": r"P",
        "aliases": ["rotationperiod", "prot", "p0rot"],
        "units": "s",
        "sign": ">= 0",
    },
    "rotationfdot": {
        "value": "Source rotational frequency derivative",
        "latex_string": r"\dot{f}_{\rm rot}",
        "aliases": ["rotationfdot", "frotdot", "f1rot", "f1spin"],
        "units": "Hz / s",
        "sign": None,
    },
    "gwfdot": {
        "value": "Gravitational-wave frequency derivative",
        "latex_string": r"\dot{f}_{\rm gw}",
        "aliases": ["gwfdot", "fdotgw", "f1gw"],
        "units": "Hz / s",
        "sign": None,
    },
    "rotationpdot": {
        "value": "Source rotational period derivative",
        "latex_string": r"\dot{P}",
        "aliases": ["pdot", "p0dot", "rotationpdot"],
        "units": "s / s",
        "sign": None,
    },
    "rotationfddot": {
        "value": "Source rotational frequency second derivative",
        "latex_string": r"\ddot{f}_{\rm rot}",
        "aliases": ["rotationfddot", "frotddot", "f2rot", "f2spin"],
        "units": "Hz / s / s",
        "sign": None,
    },
    "gwfdot": {
        "value": "Gravitational-wave second frequency derivative",
        "latex_string": r"\ddot{f}_{\rm gw}",
        "aliases": ["gwfddot", "fddotgw", "f2gw"],
        "units": "Hz / s / s",
        "sign": None,
    },
    "momentofinertia": {
        "value": "Principle moment of inertia about the rotation axis",
        "latex_string": r"I_{zz}",
        "aliases": ["momentofinertia", "izz", "i38"],
        "units": "kg m^2",
        "sign": ">= 0",
    },
    "distance": {
        "value": "Distance to the source",
        "latex_string": "d",
        "aliases": ["distance", "d", "r"],
        "units": "kpc",
        "sign": ">= 0",
    },
}


def equations(equation, **kwargs):
    """
    This function holds information on the set of equations that are defined.

    The equation definition information is held in a class attribute called
    equation_info. If defining a new equation the information should be added
    into this dictionary. The dictionary keys give the equation's "name", with
    which it must be referred to. The values are dictionaries containing the
    following keys:

    "default_fiducial_values": dict - a dictionary with keys for all variable
    parameters in the equation. Each key gives a 2-tuple containing the default
    value of that variable (with appropriate :class:`astropy.unit.Unit`) and
    the variable's exponent as a string. The variable parameter names must be
    consistent with the names in the :var:`~cweqgen.equations.ALLOWED_VALUES`
    dictionary.

    "equation_constants": list - a list of 2-tuples containing the constants in
    the equation and their exponents. These should both be string values.

    "additional_values": list - a list of variable names that can be used to
    derive a subset of the variables in the "default_fiducial_values", i.e.,
    if the equation requires the "rotationfrequency" then this could contain
    "rotationperiod", which can instead be used to derive the rotation
    frequency. The variable parameter names must be consistent with the names
    in the :var:`~cweqgen.equations.ALLOWED_VALUES` dictionary.

    "converters": dict - a dictionary of coversion functions to convert the
    additional values into the required values.

    The current equations that are defined are:

    "h0": the gravitational-wave amplitude for a signal emitted from the l=m=2
    mass quadrupole mode.

    """

    equation_info = {
        "h0": {
            "default_fiducial_values": {
                "ellipticity": (1e-6, "1"),
                "momentofinertia": (1e38 * u.Unit("kg m^2"), "1"),
                "rotationfrequency": (100 * u.Hz, "2"),
                "distance": (1 * u.kpc, "-1"),
            },
            "constants": [
                ("16", "1"),
                ("pi", "2"),
                ("G", "1"),
                ("c", "-4"),
            ],
            "additional_values": ["gwfrequency", "rotationperiod"],
            "converters": {
                "rotationfrequency": convert_to_rotation_frequency,
            },
        },
        "h0spindown": {
            "default_fiducial_values": {
                "momentofinertia": (1e38 * u.Unit("kg m^2"), "1/2"),
                "rotationfrequency": (100 * u.Hz, "-1/2"),
                "rotationfdot": (-1e-11 * u.Hz / u.s, "1/2"),
                "distance": (1 * u.kpc, "-1"),
            },
            "constants": [
                ("5", "1/2"),
                ("2", "-1/2"),
                ("G", "1/2"),
                ("c", "-3/2"),
            ],
            "additional_values": [
                "gwfrequency",
                "rotationperiod",
                "gwfdot",
                "rotationpdot",
            ],
            "converters": {
                "rotationfrequency": convert_to_rotation_frequency,
                "rotationfdot": convert_to_rotation_fdot,
            },
        },
        "ellipticityspindown": {
            "default_fiducial_values": {
                "momentofinertia": (1e38 * u.Unit("kg m^2"), "1/2"),
                "rotationfrequency": (100 * u.Hz, "-1/2"),
                "rotationfdot": (-1e-11 * u.Hz / u.s, "1/2"),
            },
            "constants": [
                ("5", "1/2"),
                ("512", "-1/2"),
                ("pi", "-2"),
                ("G", "-1/2"),
                ("c", "5/2"),
            ],
            "additional_values": [
                "gwfrequency",
                "rotationperiod",
                "gwfdot",
                "rotationpdot",
            ],
            "converters": {
                "rotationfrequency": convert_to_rotation_frequency,
                "rotationfdot": convert_to_rotation_fdot,
            },
        },
    }

    if equation not in equation_info:
        raise KeyError(f"Equation '{equation}' is not currently defined")

    # set values for required equation
    eqinfo = equation_info[equation]
    kwargs["equation"] = equation
    kwargs["default_fiducial_values"] = eqinfo["default_fiducial_values"]
    kwargs["constants"] = eqinfo["constants"]
    kwargs["additional_values"] = eqinfo["additional_values"]
    kwargs["converters"] = eqinfo["converters"]

    class _EquationBase:
        __metaclass__ = abc.ABCMeta

        def __init__(self, **kwargs):
            self.equation_name = kwargs.pop("equation")
            
            # dictionary to hold default fiducial values
            self.default_fiducial_values = kwargs.pop("default_fiducial_values")

            # list containing additional keyword values that can be used
            self.additional_values = kwargs.pop("additional_values")

            # dictionary of functions to convert from additional values into
            # required values (keyed on the required values)
            self.converters = kwargs.pop("converters")

            # a list of tuples containing equation constants
            self.constants = kwargs.pop("constants")

            self.parse_kwargs(**kwargs)

            self.latex_name = ALLOWED_VALUES[self.equation_name][
                "latex_string"
            ]  # lhs of equation
            self.description = ALLOWED_VALUES[self.equation_name]["value"]

            # get simple format reference
            self.reference_string = REFERENCES[self.equation_name]["short"]

            # equation number in reference
            self.reference_eqno = REFERENCES[self.equation_name]["eqno"]

            # URL of reference in ADS
            self.reference_adsurl = REFERENCES[self.equation_name]["adsurl"]

            # BibTeX for reference (from ADS)
            self.reference_bibtex = REFERENCES[self.equation_name]["bibtex"]

        def parse_kwargs(self, **kwargs):
            """
            Get the required values for the specific equation.
            """

            self.values = {}

            for key in list(self.default_fiducial_values.keys()) + self.additional_values:
                # check aliases
                for alias in ALLOWED_VALUES[key]["aliases"]:
                    if alias in kwargs:
                        value = kwargs[alias]

                        # check value has compatible units
                        if (
                            not isinstance(value, Quantity)
                            and ALLOWED_VALUES[key]["units"] is not None
                        ):
                            value *= u.Unit(ALLOWED_VALUES[key]["units"])
                        elif (
                            isinstance(value, Quantity)
                            and ALLOWED_VALUES[key]["units"] is not None
                        ):
                            try:
                                _ = value.to(ALLOWED_VALUES[key]["units"])
                            except (u.UnitConversionError, ValueError) as e:
                                raise IOError(
                                    f"{ALLOWED_VALUES[key]['value']} units are not compatible:\n{e}"
                                )

                        # check value has correct sign
                        if ALLOWED_VALUES[key]["sign"] is not None:
                            if not eval(str(value.value) + ALLOWED_VALUES[key]["sign"]):
                                raise ValueError(
                                    f"{ALLOWED_VALUES[key]['value']} does not have the correct sign"
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
                constant *= eval(str(const[0]) + "**(" + const[1] + ")")

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
                exp = self.default_fiducial_values[key][1]

                fiducial *= abs(val) ** float(eval(exp))

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
                exp = const[1].strip("-")

                if exp not in constfractions:
                    constfractions[exp] = ["", ""]

                if eval(const[1]) < 0:
                    # denominator value
                    constfractions[exp][1] += " " + cv
                else:
                    # numerator value
                    constfractions[exp][0] += " " + cv

            # construct constant latex string
            conststr = ""
            for exp in constfractions:
                if len(constfractions[exp][0]) == 0:
                    constnumstr = "1"
                else:
                    constnumstr = constfractions[exp][0]

                constdenstr = constfractions[exp][1]

                lbrace, rbrace = ("", "") if exp == "1" else (r"\left(", r"\right)")
                expstr = "" if exp == "1" else (f"^{{{exp}}}")

                if len(constdenstr) == 0:
                    conststr += f"{constnumstr}{expstr}"
                else:
                    conststr += (
                        rf"{lbrace}\frac{{{constnumstr}}}{{{constdenstr}}}{rbrace}{expstr}"
                    )

            latex_equation += conststr

            # use default fiducial value keys for the rest of the equation
            varnumstr = ""  # variables numerator
            vardenstr = ""  # variables denominator

            varfractions = {}

            for key in self.default_fiducial_values:
                if key in ALLOWED_VALUES:
                    varlatex = ALLOWED_VALUES[key]["latex_string"]
                else:
                    varlatex = key

                # get exponent
                exp = self.default_fiducial_values[key][1].strip("-")

                if exp not in varfractions:
                    varfractions[exp] = ["", ""]

                if eval(self.default_fiducial_values[key][1]) < 0:
                    # denominator value
                    varfractions[exp][1] += " " + varlatex
                else:
                    # numerator value
                    varfractions[exp][0] += " " + varlatex

            # construct variables latex string
            varstr = ""
            for exp in varfractions:
                if len(varfractions[exp][0]) == 0:
                    varnumstr = "1"
                else:
                    varnumstr = varfractions[exp][0]

                vardenstr = varfractions[exp][1]

                lbrace, rbrace = ("", "") if exp == "1" else (r"\left(", r"\right)")
                expstr = "" if exp == "1" else (f"^{{{exp}}}")

                if len(vardenstr) == 0:
                    varstr += f"{varnumstr}{expstr}"
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

            latex_equation += coeff.to_string(precision=(dp + 1), format="latex").replace(
                "$", ""
            )

            if brackets not in ["()", "{}", "[]", None]:
                raise ValueError(f"Bracket type {brackets} is not recognised")

            lbrace = r"\left" + brackets[0] if brackets is not None else ""
            rbrace = r"\right" + brackets[1] if brackets is not None else ""

            fiducial = ""

            for key in self.default_fiducial_values:
                if key in ALLOWED_VALUES:
                    varlatex = ALLOWED_VALUES[key]["latex_string"]
                else:
                    varlatex = key

                # get exponent
                exp = self.default_fiducial_values[key][1]

                # get value
                if key in self.values:
                    # use provided value
                    val = Quantity(self.values[key])
                else:
                    val = Quantity(self.default_fiducial_values[key][0])

                if eval(exp) < 0:
                    numerator = val.to_string(precision=(dp + 1), format="latex").replace(
                        "$", ""
                    )
                    denominator = varlatex
                else:
                    denominator = val.to_string(precision=(dp + 1), format="latex").replace(
                        "$", ""
                    )
                    numerator = varlatex

                if eval(exp.strip("-")) != 1.0:
                    expstr = "^{" + exp.strip("-") + "}"
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

        def __str__(self):
            return self.equation()

    # update the class docstring
    _EquationBase.__doc__ = DOCSTRINGS[equation]

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
    def __init__(self, latexstring,  dpi=200):
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
        self.fig, ax = plt.subplots(1)
        
        try:
            t = ax.text(
                0.05,
                0.5,
                "$" + self.text + "$",
                usetex=True,
            )
        except RuntimeError:
            t = ax.text(
                0.0,
                0.5,
                "$" + self.text + "$",
            )
        ax.axis("off")

        # crop figure to tight around the text
        bb = t.get_window_extent(renderer=self.fig.canvas.get_renderer())
        transf = ax.transData.inverted()
        bb_datacoords = bb.transformed(transf)

        rect = bb_datacoords.get_points().flatten()
        rect[3] += 0.1 * rect[3]  # add 10% on to upper value for some reason!

        tight_params = {"rect": rect}

        self.fig.set_dpi(self.dpi)
        self.fig.set_tight_layout(tight_params)

    def savefig(self, **kwargs):
        kwargs.setdefault("dpi", self.dpi)
        return self.fig.savefig(**kwargs)
