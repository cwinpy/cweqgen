import astropy.units as u
from astropy.units.quantity import Quantity

from .converters import *


#: allowed variables for equations
ALLOWED_VARIABLES = {
    "h0": {
        "description": "Gravitational wave amplitude",
        "latex_string": r"h_0",
        "aliases": ["h0", "h_0"],
        "units": None,
        "sign": ">= 0",
    },
    "ellipticity": {
        "description": "Neutron star ellipticity",
        "latex_string": r"\varepsilon",
        "aliases": [
            "ellipticity",
            "ell",
            "eps",
            "epsilon",
            "𝜀",
        ],
        "units": None,
        "sign": ">= 0",
    },
    "massquadrupole": {
        "description": "Mass quadrupole moment (l=m=2)",
        "latex_string": r"Q_{22}",
        "aliases": ["massquadrupole", "q22", "q_22"],
        "units": "kg m^2",
        "sign": ">= 0",
    },
    "rotationfrequency": {
        "description": "Source rotational frequency",
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
        "description": "Gravitational-wave frequency",
        "latex_string": r"f_{\rm gw}",
        "aliases": ["gwfrequency", "fgw", "f0gw"],
        "units": "Hz",
        "sign": ">= 0",
    },
    "rotationperiod": {
        "description": "Source rotational period",
        "latex_string": r"P",
        "aliases": ["rotationperiod", "prot", "p0rot"],
        "units": "s",
        "sign": ">= 0",
    },
    "rotationfdot": {
        "description": "Source rotational frequency derivative",
        "latex_string": r"\dot{f}_{\rm rot}",
        "aliases": ["rotationfdot", "frotdot", "f1rot", "f1spin"],
        "units": "Hz / s",
        "sign": None,
    },
    "gwfdot": {
        "description": "Gravitational-wave frequency derivative",
        "latex_string": r"\dot{f}_{\rm gw}",
        "aliases": ["gwfdot", "fdotgw", "f1gw"],
        "units": "Hz / s",
        "sign": None,
    },
    "rotationpdot": {
        "description": "Source rotational period derivative",
        "latex_string": r"\dot{P}",
        "aliases": ["pdot", "p0dot", "rotationpdot"],
        "units": "s / s",
        "sign": None,
    },
    "rotationfddot": {
        "description": "Source rotational frequency second derivative",
        "latex_string": r"\ddot{f}_{\rm rot}",
        "aliases": ["rotationfddot", "frotddot", "f2rot", "f2spin"],
        "units": "Hz / s / s",
        "sign": None,
    },
    "gwfdot": {
        "description": "Gravitational-wave second frequency derivative",
        "latex_string": r"\ddot{f}_{\rm gw}",
        "aliases": ["gwfddot", "fddotgw", "f2gw"],
        "units": "Hz / s / s",
        "sign": None,
    },
    "momentofinertia": {
        "description": "Principal moment of inertia about the rotation axis",
        "latex_string": r"I_{zz}",
        "aliases": ["momentofinertia", "izz", "i38"],
        "units": "kg m^2",
        "sign": ">= 0",
    },
    "distance": {
        "description": "Distance to the source",
        "latex_string": "d",
        "aliases": ["distance", "d", "r"],
        "units": "kpc",
        "sign": ">= 0",
    },
    "brakingindex": {
        "description": "The braking index of a pulsar",
        "latex_string": "n",
        "aliases": ["brakingindex", "n"],
        "units": None,
        "sign": None,
    },
    "characteristicage": {
        "description": "The characteristic age of a pulsar",
        "latex_string": r"\tau",
        "aliases": ["characteristicage", "tau", "𝜏"],
        "units": "yr",
        "sign": ">= 0",
    },
    "luminosity": {
        "description": "The luminosity of a source",
        "latex_string": "L",
        "aliases": ["luminosity", "l", "spindownluminosity", "gwluminosity", "lgw", "lsd"],
        "units": "W",
        "sign": ">= 0",
    }
}


class EqDict(dict):
    """
    Dictionary class to hold equation definitions.
    """

    def __setitem__(self, key, subdict):
        super(EqDict, self).__setitem__(key, {})

        try:
            self[key]["description"] = subdict["description"]
        except KeyError:
            raise KeyError("Equation dictionary must contain a 'description'")

        try:
            self[key]["latex_string"] = subdict["latex_string"]
        except KeyError:
            raise KeyError("Equation dictionary must contain a 'latex_string'")

        try:
            self[key]["default_fiducial_values"] = subdict["default_fiducial_values"]
        except KeyError:
            raise KeyError("Equation dictionary must contain 'default_fiducial_values'")

        try:
            self[key]["parts"] = subdict["parts"]
        except KeyError:
            raise KeyError("Equation dictionary must contain 'parts'")

        self[key]["alternative_variables"] = subdict.get("alternative_variables", [])
        self[key]["converters"] = subdict.get("converters", {})

        if "reference" in subdict:
            self[key]["reference"] = subdict.get("reference")

        if "docstring" in subdict:
            self[key]["docstring"] = subdict.get("docstring").format(
                **self.get_defaults(key)
            )

    def get_defaults(self, key):
        """
        Get the default values for a given equation.

        Parameters
        ----------
        key: str
            The equation name.

        Returns
        -------
        defaults: dict
            A dictionary of default values.
        """

        defaults = {
            key: Quantity(item)._repr_latex_().replace("$", "")
            for key, item in self[key]["default_fiducial_values"].items()
        }

        defaults["name"] = key

        return defaults


#: equation definitions
EQN_DEFINITIONS = EqDict()

EQN_DEFINITIONS["h0"] = {
    "description": "Gravitational wave amplitude",
    "variable": "h0",
    "latex_string": "h_0",
    "default_fiducial_values": {
        "ellipticity": 1e-6,
        "momentofinertia": 1e38 * u.Unit("kg m^2"),
        "rotationfrequency": 100 * u.Hz,
        "distance": 1 * u.kpc,
    },
    "parts": [
        ("16", "1"),
        ("pi", "2"),
        ("G", "1"),
        ("c", "-4"),
        ("ellipticity", "1"),
        ("momentofinertia", "1"),
        ("rotationfrequency", "2"),
        ("distance", "-1"),
    ],
    "alternative_variables": ["gwfrequency", "rotationperiod"],
    "converters": {
        "rotationfrequency": convert_to_rotation_frequency,
    },
    "reference": {
        "short": "Jaranowski, P., Krolak, A., & Schutz, B. F. 1998, PhRvD, 58, 063001",
        "adsurl": "https://ui.adsabs.harvard.edu/abs/1998PhRvD..58f3001J/abstract",
        "eqno": "23",
        "bibtex": r"""\
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
}""",
    },
    "docstring": """
Generate the gravitational-wave amplitude for a signal emitted from the l=m=2
mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as given in
:obj:`~cweqgen.definitions.ALLOWED_VARIABLES`, can be used instead.

:param str equation: "{name}"
:keyword float or ~astropy.units.quantity.Quantity ellipticity: The ellipticity of the source with which the calculate the GW amplitude. The default value is :math:`{ellipticity}`.
:keyword float or ~astropy.units.quantity.Quantity momentofinertia: The principal moment of inertia with which the calculate the GW amplitude. If given as a float units of kg m^2 are assumed. The default value is :math:`{momentofinertia}`.
:keyword float or ~astropy.units.quantity.Quantity rotationfrequency: The rotation frequency of the source. If given as a float units of Hz are assumed. The default value is :math:`{rotationfrequency}`. If the rotational period or gravitational wave frequency are given instead then they will be converted into rotational frequency (for GW frequency it is assumed that this is twice the rotational frequency).
:keyword float or ~astropy.units.quantity.Quantity distance: The distance to the source. If given as a float units of kpc are assumed. The default value is :math:`{distance}`.
""",
}

EQN_DEFINITIONS["spindownluminosity"] = {
    "description": "The spin-down luminosity of a pulsar",
    "variable": "luminosity",
    "latex_string": r"L_{\rm sd}",
    "default_fiducial_values": {
        "momentofinertia": 1e38 * u.Unit("kg m^2"),
        "rotationfrequency": 100 * u.Hz,
        "rotationfdot": -1e-11 * u.Hz / u.s,
    },
    "parts": [
        ("4", "1"),
        ("pi", "2"),
        ("momentofinertia", "1"),
        ("rotationfrequency", "1"),
        ("rotationfdot", "1"),
    ],
    "alternative_variables": [
        "gwfrequency",
        "rotationperiod",
        "gwfdot",
        "rotationpdot",
    ],
    "converters": {
        "rotationfrequency": convert_to_rotation_frequency,
        "rotationfdot": convert_to_rotation_fdot,
    },
    "reference": {  # this is just an example reference for the braking index (there will be earlier references!)
        "short": "Condon, J. J. and Ransom, S. M., 2016, Essential Radio Astronomy",
        "adsurl": "https://ui.adsabs.harvard.edu/abs/2016era..book.....C/abstract",
        "eqno": "6.35",
        "bibtex": r"""\
@BOOK{2016era..book.....C,
       author = {{Condon}, James J. and {Ransom}, Scott M.},
        title = "{Essential Radio Astronomy}",
         year = 2016,
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016era..book.....C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}""",
    },
    "docstring": """
Generate the spin-down luminosity of a pulsar.

For the optional input keyword parameters below a range of aliases, as given in
:obj:`~cweqgen.definitions.ALLOWED_VARIABLES`, can be used instead.

:param str equation: "{name}"
:keyword float or ~astropy.units.quantity.Quantity momentofinertia: The principal moment of inertia with which the calculate the GW amplitude. If given as a float units of kg m^2 are assumed. The default value is :math:`{momentofinertia}`.
:keyword float or ~astropy.units.quantity.Quantity rotationfrequency: The rotation frequency of the source. If given as a float units of Hz are assumed. The default value is :math:`{rotationfrequency}`. If the rotational period or gravitational wave frequency are given instead then they will be converted into rotational frequency (for GW frequency it is assumed that this is twice the rotational frequency).
:keyword float or ~astropy.units.quantity.Quantity rotationfdot: The first rotational frequency derivative (i.e. the spin-down). If given as a float units of Hz/s are assumed. The default value is :math:`{rotationfdot}`. If the rotational period derivative or gravitational wave frequency derivative is given instead then they will be converted into rotational frequency derivative.
""",
}

EQN_DEFINITIONS["gwluminosity"] = {
    "description": "The gravitational-wave luminosity of a pulsar",
    "variable": "luminosity",
    "latex_string": r"L_{\rm gw}",
    "default_fiducial_values": {
        "momentofinertia": 1e38 * u.Unit("kg m^2"),
        "rotationfrequency": 100 * u.Hz,
        "ellipticity": 1e-6,
    },
    "parts": [
        ("2048/5", "1"),
        ("pi", "6"),
        ("G", "1"),
        ("c", "-5"),
        ("momentofinertia", "2"),
        ("ellipticity", "2"),
        ("rotationfrequency", "6"),
    ],
    "alternative_variables": [
        "gwfrequency",
        "rotationperiod",
    ],
    "converters": {
        "rotationfrequency": convert_to_rotation_frequency,
    },
    "reference": {
        "short": "Aasi, A., et al. 2014, ApJ, 785, 119",
        "adsurl": "https://ui.adsabs.harvard.edu/abs/2014ApJ...785..119A/abstract",
        "eqno": "4",
        "bibtex": r"""\
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
}""",
    },
    "docstring": """
Generate the gravitational-wave luminosity of a pulsar.

For the optional input keyword parameters below a range of aliases, as given in
:obj:`~cweqgen.definitions.ALLOWED_VARIABLES`, can be used instead.

:param str equation: "{name}"
:keyword float or ~astropy.units.quantity.Quantity ellipticity: The ellipticity of the source with which the calculate the GW amplitude. The default value is :math:`{ellipticity}`.
:keyword float or ~astropy.units.quantity.Quantity momentofinertia: The principal moment of inertia with which the calculate the GW amplitude. If given as a float units of kg m^2 are assumed. The default value is :math:`{momentofinertia}`.
:keyword float or ~astropy.units.quantity.Quantity rotationfrequency: The rotation frequency of the source. If given as a float units of Hz are assumed. The default value is :math:`{rotationfrequency}`. If the rotational period or gravitational wave frequency are given instead then they will be converted into rotational frequency (for GW frequency it is assumed that this is twice the rotational frequency).
""",
}

EQN_DEFINITIONS["h0spindown"] = {
    "description": "Gravitational wave amplitude spin-down limit",
    "variable": "h0",
    "latex_string": r"h_0^{\rm sd}",
    "default_fiducial_values": {
        "momentofinertia": 1e38 * u.Unit("kg m^2"),
        "rotationfrequency": 100 * u.Hz,
        "rotationfdot": -1e-11 * u.Hz / u.s,
        "distance": 1 * u.kpc,
    },
    "parts": [
        ("5/2", "1/2"),
        ("G", "1/2"),
        ("c", "-3/2"),
        ("momentofinertia", "1/2"),
        ("rotationfrequency", "-1/2"),
        ("rotationfdot", "1/2"),
        ("distance", "-1"),
    ],
    "alternative_variables": [
        "gwfrequency",
        "rotationperiod",
        "gwfdot",
        "rotationpdot",
    ],
    "converters": {
        "rotationfrequency": convert_to_rotation_frequency,
        "rotationfdot": convert_to_rotation_fdot,
    },
    "reference": {
        "short": "Aasi, A., et al. 2014, ApJ, 785, 119",
        "adsurl": "https://ui.adsabs.harvard.edu/abs/2014ApJ...785..119A/abstract",
        "eqno": "5",
        "bibtex": r"""\
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
}""",
    },
    "docstring": """
Generate the equation for the spin-down limit on the gravitational wave
amplitude for a signal emitted from the l=m=2 mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as given in
:obj:`~cweqgen.definitions.ALLOWED_VARIABLES`, can be used instead.

:param str equation: "{name}"
:keyword float or ~astropy.units.quantity.Quantity momentofinertia: The principal moment of inertia with which the calculate the GW amplitude. If given as a float units of kg m^2 are assumed. The default value is :math:`{momentofinertia}`.
:keyword float or ~astropy.units.quantity.Quantity rotationfrequency: The rotation frequency of the source. If given as a float units of Hz are assumed. The default value is :math:`{rotationfrequency}`. If the rotational period or gravitational wave frequency are given instead then they will be converted into rotational frequency (for GW frequency it is assumed that this is twice the rotational frequency).
:keyword float or ~astropy.units.quantity.Quantity distance: The distance to the source. If given as a float units of kpc are assumed. The default value is :math:`{distance}`.
:keyword float or ~astropy.units.quantity.Quantity rotationfdot: The first rotational frequency derivative (i.e. the spin-down). If given as a float units of Hz/s are assumed. The default value is :math:`{rotationfdot}`. If the rotational period derivative or gravitational wave frequency derivative is given instead then they will be converted into rotational frequency derivative.
""",
}

EQN_DEFINITIONS["ellipticityspindown"] = {
    "description": "Spin-down limit for neutron star ellipticity",
    "variable": "ellipticity",
    "latex_string": r"\varepsilon^{\rm sd}",
    "default_fiducial_values": {
        "momentofinertia": 1e38 * u.Unit("kg m^2"),
        "rotationfrequency": 100 * u.Hz,
        "rotationfdot": -1e-11 * u.Hz / u.s,
    },
    "parts": [
        ("5/512", "1/2"),
        ("pi", "-2"),
        ("G", "-1/2"),
        ("c", "5/2"),
        ("momentofinertia", "-1/2"),
        ("rotationfrequency", "-5/2"),
        ("rotationfdot", "1/2"),
    ],
    "alternative_variables": [
        "gwfrequency",
        "rotationperiod",
        "gwfdot",
        "rotationpdot",
    ],
    "converters": {
        "rotationfrequency": convert_to_rotation_frequency,
        "rotationfdot": convert_to_rotation_fdot,
    },
    "reference": {
        "short": "Abbott, B. P., et al. 2019, ApJ, 879, 10",
        "adsurl": "https://ui.adsabs.harvard.edu/abs/2019ApJ...879...10A/abstract",
        "eqno": "A9",
        "bibtex": r"""\
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
}""",
    },
    "docstring": """
Generate the equation for the spin-down limit on the source ellipticity
for a signal emitted from the l=m=2 mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as
given in :obj:`cweqgen.definitions.ALLOWED_VARIABLES`, can be used instead.

:param str equation: "{name}"
:keyword float or ~astropy.units.quantity.Quantity momentofinertia: The principle moment of inertia with which the calculate the spin-down limit. If given as a float units of kg m^2 are assumed. The default value is :math:`{momentofinertia}`.
:keyword float or ~astropy.units.quantity.Quantity rotationfrequency: The rotation frequency of the source. If given as a float units of Hz are assumed. The default value is :math:`{rotationfrequency}`. If the rotational period or gravitational wave frequency are given instead then they will be converted into rotational frequency (for GW frequency it is assumed that this is twice the rotational frequency).
:keyword float or ~astropy.units.quantity.Quantity rotationfdot: The first rotational frequency derivative (i.e. the spin-down). If given as a float units of Hz/s are assumed. The default value is :math:`{rotationfdot}`. If the rotational period derivative or gravitational wave frequency derivative is given instead then they will be converted into rotational frequency derivative.
""",
}

EQN_DEFINITIONS["brakingindex"] = {
    "description": "The braking index of a pulsar",
    "variable": "brakingindex",
    "latex_string": "n",
    "default_fiducial_values": {
        "rotationfrequency": 50 * u.Hz,
        "rotationfddot": 1e-23 * u.Hz / (u.s ** 2),
        "rotationfdot": -1e-11 * u.Hz / u.s,
    },
    "parts": [
        ("rotationfrequency", "1"),
        ("rotationfddot", "1"),
        ("rotationfdot", "-2"),
    ],
    "alternative_variables": [
        "gwfrequency",
        "rotationperiod",
        "gwfdot",
        "rotationpdot",
    ],
    "converters": {
        "rotationfrequency": convert_to_rotation_frequency,
        "rotationfdot": convert_to_rotation_fdot,
    },
    "reference": {  # this is just an example reference for the braking index (there will be earlier references!)
        "short": "Condon, J. J. and Ransom, S. M., 2016, Essential Radio Astronomy",
        "adsurl": "https://ui.adsabs.harvard.edu/abs/2016era..book.....C/abstract",
        "eqno": "6.35",
        "bibtex": r"""\
@BOOK{2016era..book.....C,
       author = {{Condon}, James J. and {Ransom}, Scott M.},
        title = "{Essential Radio Astronomy}",
         year = 2016,
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016era..book.....C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}""",
    },
    "docstring": """
Generate the equation for the braking index of a pulsar.

For the optional input keyword parameters below a range of aliases, as
given in :obj:`cweqgen.definitions.ALLOWED_VARIABLES`, can be used instead.

:param str equation: "{name}"
:keyword float or ~astropy.units.quantity.Quantity rotationfrequency: The rotation frequency of the source. If given as a float units of Hz are assumed. The default value is :math:`{rotationfrequency}`. If the rotational period or gravitational wave frequency are given instead then they will be converted into rotational frequency (for GW frequency it is assumed that this is twice the rotational frequency).
:keyword float or ~astropy.units.quantity.Quantity rotationfdot: The first rotational frequency derivative (i.e. the spin-down). If given as a float units of Hz/s are assumed. The default value is :math:`{rotationfdot}`. If the rotational period derivative or gravitational wave frequency derivative is given instead then they will be converted into rotational frequency derivative.
:keyword float or ~astropy.units.quantity.Quantity rotationfddot: The second rotational frequency derivative. If given as a float units of Hz/s^2 are assumed. The default value is :math:`{rotationfddot}`.
""",
}

EQN_DEFINITIONS["characteristicage"] = {
    "description": "The characteristic age of a pulsar",
    "variable": "characteristicage",
    "latex_string": r"\tau",
    "default_fiducial_values": {
        "brakingindex": 3,
        "rotationperiod": 0.01 * u.s,
        "rotationpdot": 1e-15 * u.s / u.s,
    },
    "parts": [
        ("rotationperiod", "1"),
        ("rotationpdot", "-1"),
        ("brakingindex - 1", "-1"),
    ],
    "alternative_variables": ["gwfrequency", "rotationfrequency", "rotationfdot"],
    "converters": {
        "rotationperiod": convert_to_rotation_period,
        "rotationpdot": convert_to_rotation_pdot,
    },
    "reference": {
        "short": "Condon, J. J. and Ransom, S. M., 2016, Essential Radio Astronomy",
        "adsurl": "https://ui.adsabs.harvard.edu/abs/2016era..book.....C/abstract",
        "eqno": "6.31",
        "bibtex": r"""\
@BOOK{2016era..book.....C,
       author = {{Condon}, James J. and {Ransom}, Scott M.},
        title = "{Essential Radio Astronomy}",
         year = 2016,
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016era..book.....C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}""",
    },
    "docstring": """
Generate the pulsar's characteristic age.

For the optional input keyword parameters below a range of aliases, as given in
:obj:`~cweqgen.definitions.ALLOWED_VARIABLES`, can be used instead.

:param str equation: "{name}"
:keyword float or ~astropy.units.quantity.Quantity rotationperiod: The rotation period of a pulsar. If given as a float units of s are assumed. The default value is :math:`{rotationperiod}`.
:keyword float or ~astropy.units.quantity.Quantity rotationpdot: The first derivative of the rotation period of a pulsar. If given as a float units of s/s are assumed. The default value is :math:`{rotationpdot}`.
:keyword float or ~astropy.units.quantity.Quantity brakingindex: The braking index of the pulsar. The default value is :math:`{brakingindex}`.
""",
}