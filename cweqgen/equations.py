import numpy as np
from numpy import pi

import astropy.units as u
from astropy.units.format import Latex
from astropy.units.quantity import Quantity
from astropy.constants import G, c

import matplotlib as mpl
from matplotlib import pyplot as plt


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


class EquationLaTeXString:
    def __init__(self, latexstring):
        """
        Class to hold a LaTeX equation string. It has a _repr_latex_ method to
        hook into the Jupyter notebook rich display system
        https://ipython.readthedocs.io/en/stable/config/integrating.html#rich-display
        and show the resulting string as a LaTeX equation.
        """

        self.text = latexstring

    def __repr__(self):
        return self.text

    def __str__(self):
        return self.text

    def _repr_latex_(self):
        return "$" + self.text + "$"


class EquationBase:
    def __init__(self, **kwargs):
        """
        A base class for generating equations.
        """

        # dictionary to hold default fiducial values
        self.default_fiducial_values = {}

        # list containing additional keyword values that can be used
        self.additional_values = []

        # dictionary of functions to convert from additional values into
        # required values (keyed on the required values)
        self.converters = {}

        # a list of tuples containing equation constants
        self.equation_constants = []

        self.parse_kwargs(**kwargs)

        self.latex_name = ""  # lhs of equation

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
                            qval = value.to(ALLOWED_VALUES[key]["units"])
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

        for const in self.equation_constants:
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

        for const in self.equation_constants:
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
            try:
                with mpl.rc_context({"text.usetex": True}):
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)
                    ax.text(
                        0.0,
                        0.0,
                        "$" + latex_equation + "$",
                        horizontalalignment="left",
                        verticalalignment="center",
                        transform=ax.transAxes,
                    )
                    ax.axis("off")
                    fig.tight_layout()
            except RuntimeError:
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.text(
                    0.0,
                    0.0,
                    "$" + latex_equation + "$",
                    horizontalalignment="left",
                    verticalalignment="center",
                    transform=ax.transAxes,
                )
                ax.axis("off")
                fig.tight_layout()
            return fig
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
            try:
                with mpl.rc_context({"text.usetex": True}):
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)
                    ax.text(
                        0.0,
                        0.0,
                        "$" + latex_equation + "$",
                        horizontalalignment="left",
                        verticalalignment="center",
                        transform=ax.transAxes,
                    )
                    ax.axis("off")
                    fig.tight_layout()
            except RuntimeError:
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.text(
                    0.0,
                    0.0,
                    "$" + latex_equation + "$",
                    horizontalalignment="left",
                    verticalalignment="center",
                    transform=ax.transAxes,
                )
                ax.axis("off")
                fig.tight_layout()
            return fig
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


class Equationh0(EquationBase):
    def __init__(self, **kwargs):
        """
        Generate the equation for the gravitational wave amplitude for a signal
        emitted from the l=m=2 mass quadrupole mode.

        For the optional input keyword parameters below a range of aliases, as
        given in ALLOWED_ALIASES, can be used instead.

        Parameters
        ----------
        ellipticity: float
            The ellipticity of the source with which the calculate the GW
            amplitude. If given as a float units of kg m^2 are assumed. The
            default value used if none is given is 1e-6.
        momentofinertia: float, Quantity
            The principle moment of inertia with which the calculate the GW
            amplitude. If given as a float units of kg m^2 are assumed. The
            default value used if none is given is 1e38 kg m^2.
        rotationfrequency: float, Quantity
            The rotation frequency of the source. If given as a float units of
            Hz are assumed. The default value if none is given is 100 Hz. If the
            rotational period or gravitational wave frequency are given instead
            then they will be converted into rotational frequency (for GW
            frequency it is assumed that this is twice the rotational
            frequency).
        distance: float, Quantity
            The distance to the source. If given as a float units of kpc are
            assumed. The default value used if none is given is 1 kpc.
        """

        # some default fiducial values (tuples of value and exponenent)
        self.default_fiducial_values = {
            "ellipticity": (1e-6, "1"),
            "momentofinertia": (1e38 * u.Unit("kg m^2"), "1"),
            "rotationfrequency": (100 * u.Hz, "2"),
            "distance": (1 * u.kpc, "-1"),
        }

        # dictionary containing additional values that can be used
        self.additional_values = ["gwfrequency", "rotationperiod"]

        # dictionary of functions to convert from additional values into require
        # values
        self.converters = {
            "rotationfrequency": self._convert_to_rotation_frequency,
        }

        self.parse_kwargs(**kwargs)

        # constants in the equation (tuples with the value, its exponent)
        self.equation_constants = [
            ("16", "1"),
            ("pi", "2"),
            ("G", "1"),
            ("c", "-4"),
        ]

        self.latex_name = ALLOWED_VALUES["h0"]["latex_string"]

        # get simple format reference
        self.reference_string = (
            "Jaranowski, P., Krolak, A., & Schutz, B. F. 1998, PhRvD, 58, 063001",
        )
        self.reference_eqno = "23"  # equation number in reference

        # URL of reference in ADS
        self.reference_ads = (
            "https://ui.adsabs.harvard.edu/abs/1998PhRvD..58f3001J/abstract"
        )

        # BibTeX for reference (from ADS)
        self.reference_bibtex = r"""\
@ARTICLE{1998PhRvD..58f3001J,
       author = {{Jaranowski}, Piotr and {Kr{\'o}lak}, Andrzej and {Schutz}, Bernard F.},
        title = "{Data analysis of gravitational-wave signals from spinning neutron stars: The signal and its detection}",
      journal = {\prd},
     keywords = {95.55.Ym, 04.80.Nn, 95.75.Pq, 97.60.Gb, Gravitational radiation detectors, mass spectrometers, and other instrumentation and techniques, Gravitational wave detectors and experiments, Mathematical procedures and computer techniques, Pulsars, General Relativity and Quantum Cosmology},
         year = 1998,
        month = sep,
       volume = {58},
       number = {6},
          eid = {063001},
        pages = {063001},
          doi = {10.1103/PhysRevD.58.063001},
archivePrefix = {arXiv},
       eprint = {gr-qc/9804014},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1998PhRvD..58f3001J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}"""

    @staticmethod
    def _convert_to_rotation_frequency(gwfrequency=None, rotationperiod=None, **kwargs):
        """
        Convert the GW frequency (assumed to be twice the rotation frequency) or
        the rotation period into the rotation frequency.
        """

        if gwfrequency is not None:
            return gwfrequency / 2.0
        elif rotationperiod is not None:
            return 1.0 / rotationperiod
        else:
            raise ValueError("Required conversion parameters are not present")


class Equationh0spindown(EquationBase):
    def __init__(self, **kwargs):
        """
        Generate the equation for the spin-down limit on the gravitational wave
        amplitude for a signal emitted from the l=m=2 mass quadrupole mode.

        For the optional input keyword parameters below a range of aliases, as
        given in ALLOWED_ALIASES, can be used instead.

        Parameters
        ----------
        momentofinertia: float, Quantity
            The principle moment of inertia with which the calculate the
            spin-down limit. If given as a float units of kg m^2 are assumed.
            The default value used if none is given is 1e38 kg m^2.
        rotationfrequency: float, Quantity
            The rotation frequency of the source. If given as a float units of
            Hz are assumed. The default value if none is given is 100 Hz. If the
            rotational period or gravitational wave frequency are given instead
            then they will be converted into rotational frequency (for GW
            frequency it is assumed that this is twice the rotational
            frequency).
        rotationfdot: float, Quantity
            The first rotational frequency derivative (i.e. the spin-down). If
            given as a float units of Hz/s are assumed. The default value used
            if none is given is -1e-11 Hz/s. If the rotational period derivative
            or gravitational wave frequency derivative is given instead then
            they will be converted into rotational frequency derivative.
        distance: float, Quantity
            The distance to the source. If given as a float units of kpc are
            assumed. The default value used if none is given is 1 kpc.
        """

        # some default fiducial values (tuples of value and exponenent)
        self.default_fiducial_values = {
            "momentofinertia": (1e38 * u.Unit("kg m^2"), "1/2"),
            "rotationfrequency": (100 * u.Hz, "-1/2"),
            "rotationfdot": (-1e-11 * u.Hz / u.s, "1/2"),
            "distance": (1 * u.kpc, "-1"),
        }

        # dictionary containing additional values that can be used
        self.additional_values = [
            "gwfrequency",
            "rotationperiod",
            "gwfdot",
            "rotationpdot",
        ]

        # dictionary of functions to convert from additional values into require
        # values
        self.converters = {
            "rotationfrequency": self._convert_to_rotation_frequency,
            "rotationfdot": self._convert_to_rotation_fdot,
        }

        self.parse_kwargs(**kwargs)

        # constants in the equation (tuples with the value, its exponent)
        self.equation_constants = [
            ("5", "1/2"),
            ("2", "-1/2"),
            ("G", "1/2"),
            ("c", "-3/2"),
        ]

        self.latex_name = ALLOWED_VALUES["h0spindown"]["latex_string"]

        # get simple format reference
        self.reference_string = ("Aasi, A., et al. 2014, ApJ, 785, 119",)
        self.reference_eqno = "5"  # equation number in reference

        # URL of reference in ADS
        self.reference_ads = (
            "https://ui.adsabs.harvard.edu/abs/2014ApJ...785..119A/abstract"
        )

        # BibTeX for reference (from ADS)
        self.reference_bibtex = r"""\
        @ARTICLE{2014ApJ...785..119A,
       author = {{Aasi}, J. and others},
       title = "{Gravitational Waves from Known Pulsars: Results from the Initial Detector Era}",
      journal = {\apj},
     keywords = {gravitational waves, pulsars: general, Astrophysics - High Energy Astrophysical Phenomena, General Relativity and Quantum Cosmology},
         year = 2014,
        month = apr,
       volume = {785},
       number = {2},
          eid = {119},
        pages = {119},
          doi = {10.1088/0004-637X/785/2/119},
archivePrefix = {arXiv},
       eprint = {1309.4027},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2014ApJ...785..119A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}"""

    @staticmethod
    def _convert_to_rotation_frequency(gwfrequency=None, rotationperiod=None, **kwargs):
        """
        Convert the GW frequency (assumed to be twice the rotation frequency) or
        the rotation period into the rotation frequency.
        """

        if gwfrequency is not None:
            return gwfrequency / 2.0
        elif rotationperiod is not None:
            return 1.0 / rotationperiod
        else:
            raise ValueError("Required conversion parameters are not present")

    @staticmethod
    def _convert_to_rotation_fdot(
        gwfrequency=None,
        rotationfrequency=None,
        rotationperiod=None,
        gwfdot=None,
        rotationpdot=None,
        **kwargs,
    ):
        """
        Convert the GW frequency (assumed to be twice the rotation frequency) or
        the rotation period and GW rotation frequency deritaive or rotation
        period derivative into rotation frequency derivative.
        """

        freq = (
            gwfrequency / 2.0
            if gwfrequency is not None
            else (
                (1.0 / rotationperiod)
                if rotationperiod is not None
                else rotationfrequency
            )
        )

        if freq is not None:
            if gwfdot is not None:
                fdot = gwfdot / 2.0
            elif rotationpdot is not None:
                fdot = -rotationpdot * freq ** 2
            else:
                fdot = None

        if freq is None or fdot is None:
            raise ValueError("Required conversion parameters are not present")

        return fdot


class Equationellipticityspindown(Equationh0spindown):
    def __init__(self, **kwargs):
        """
        Generate the equation for the spin-down limit on the source ellipticity
        for a signal emitted from the l=m=2 mass quadrupole mode.

        For the optional input keyword parameters below a range of aliases, as
        given in ALLOWED_ALIASES, can be used instead.

        Parameters
        ----------
        momentofinertia: float, Quantity
            The principle moment of inertia with which the calculate the
            spin-down limit. If given as a float units of kg m^2 are assumed.
            The default value used if none is given is 1e38 kg m^2.
        rotationfrequency: float, Quantity
            The rotation frequency of the source. If given as a float units of
            Hz are assumed. The default value if none is given is 100 Hz. If the
            rotational period or gravitational wave frequency are given instead
            then they will be converted into rotational frequency (for GW
            frequency it is assumed that this is twice the rotational
            frequency).
        rotationfdot: float, Quantity
            The first rotational frequency derivative (i.e. the spin-down). If
            given as a float units of Hz/s are assumed. The default value used
            if none is given is -1e-11 Hz/s. If the rotational period derivative
            or gravitational wave frequency derivative is given instead then
            they will be converted into rotational frequency derivative.
        """

        # some default fiducial values (tuples of value and exponenent)
        self.default_fiducial_values = {
            "momentofinertia": (1e38 * u.Unit("kg m^2"), "-1/2"),
            "rotationfrequency": (100 * u.Hz, "-5/2"),
            "rotationfdot": (-1e-11 * u.Hz / u.s, "1/2"),
        }

        # dictionary containing additional values that can be used
        self.additional_values = [
            "gwfrequency",
            "rotationperiod",
            "gwfdot",
            "rotationpdot",
        ]

        # dictionary of functions to convert from additional values into require
        # values
        self.converters = {
            "rotationfrequency": self._convert_to_rotation_frequency,
            "rotationfdot": self._convert_to_rotation_fdot,
        }

        self.parse_kwargs(**kwargs)

        # constants in the equation (tuples with the value, its exponent)
        self.equation_constants = [
            ("5", "1/2"),
            ("512", "-1/2"),
            ("pi", "-2"),
            ("G", "-1/2"),
            ("c", "5/2"),
        ]

        self.latex_name = ALLOWED_VALUES["ellipticityspindown"]["latex_string"]

        # get simple format reference
        self.reference_string = ("Abbott, B. P., et al. 2019, ApJ, 879, 10",)
        self.reference_eqno = "A9"  # equation number in reference

        # URL of reference in ADS
        self.reference_ads = (
            "https://ui.adsabs.harvard.edu/abs/2019ApJ...879...10A/abstract"
        )

        # BibTeX for reference (from ADS)
        self.reference_bibtex = r"""\
        @ARTICLE{2019ApJ...879...10A,
       author = {{Aasi}, J. and others},
        title = "{Searches for Gravitational Waves from Known Pulsars at Two Harmonics in 2015-2017 LIGO Data}",
      journal = {\apj},
     keywords = {gravitational waves, pulsars: general, stars: neutron, Astrophysics - High Energy Astrophysical Phenomena, General Relativity and Quantum Cosmology},
         year = 2019,
        month = jul,
       volume = {879},
       number = {1},
          eid = {10},
        pages = {10},
          doi = {10.3847/1538-4357/ab20cb},
archivePrefix = {arXiv},
       eprint = {1902.08507},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019ApJ...879...10A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}"""
