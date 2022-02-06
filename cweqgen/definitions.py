from pathlib import Path
import pkg_resources
import yaml

import astropy.units as u
from astropy.units import Unit
from astropy.units.quantity import Quantity
from numpy import pi


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
    "angularrotationfrequency": {
        "description": "Source angular rotational frequency",
        "latex_string": r"\Omega_{\rm rot}",
        "aliases": [
            "angularrotationfrequency",
            "omegarot",
            "omrot",
            "omega0rot",
            "om0rot",
            "Ω",
            "Ωrot",
            "Ω0rot",
        ],
        "units": "rad / s",
        "sign": ">= 0",
    },
    "gwfrequency": {
        "description": "Gravitational-wave frequency",
        "latex_string": r"f_{\rm gw}",
        "aliases": ["gwfrequency", "fgw", "f0gw"],
        "units": "Hz",
        "sign": ">= 0",
    },
    "angulargwfrequency": {
        "description": "Angular gravitational-wave frequency",
        "latex_string": r"\Omega_{\rm gw}",
        "aliases": [
            "angulargwfrequency",
            "omegagw",
            "omgw",
            "omega0gw",
            "om0gw",
            "Ωgw",
            "Ω0gw",
        ],
        "units": "rad / s",
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
    "angularrotationfdot": {
        "description": "Source angular rotational frequency derivative",
        "latex_string": r"\dot{\Omega}_{\rm rot}",
        "aliases": ["angularrotationfdot", "omrotdot", "om1rot", "om1spin"],
        "units": "rad / s^2",
        "sign": None,
    },
    "gwfdot": {
        "description": "Gravitational-wave frequency derivative",
        "latex_string": r"\dot{f}_{\rm gw}",
        "aliases": ["gwfdot", "fdotgw", "f1gw"],
        "units": "Hz / s",
        "sign": None,
    },
    "angulargwfdot": {
        "description": "Gravitational-wave angular frequency derivative",
        "latex_string": r"\dot{\Omega}_{\rm gw}",
        "aliases": ["angulargwfdot", "omgwdot", "om1gw"],
        "units": "rad / s^2",
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
        "latex_string": r"\dot{f}_{\rm gw}",
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
        "aliases": [
            "luminosity",
            "l",
            "spindownluminosity",
            "gwluminosity",
            "lgw",
            "lsd",
        ],
        "units": "W",
        "sign": ">= 0",
    },
}


class EqDict(dict):
    """
    Dictionary class to hold equation definitions.
    """

    def __init__(self):
        """
        Full dictionary with equations from eqnfiles directory
        """

        # get all YAML equation files
        for eqfile in Path(pkg_resources.resource_filename("cweqgen", "eqnfiles")).glob(
            "*.yaml"
        ):
            # read in information
            with open(eqfile, "r") as fp:
                eqdata = yaml.safe_load(fp.read())

            # use file name as the equation key
            key = eqfile.name.split(".")[0]

            self[key] = eqdata

    def __setitem__(self, key, subdict):
        super(EqDict, self).__setitem__(key, {})

        try:
            self[key]["description"] = subdict["description"]
        except KeyError:
            raise KeyError("Equation dictionary must contain a 'description'")

        try:
            self[key]["variable"] = subdict["variable"]
        except KeyError:
            raise KeyError("Equation dictionary must contain a 'variable'")

        try:
            self[key]["latex_string"] = subdict["latex_string"]
        except KeyError:
            raise KeyError("Equation dictionary must contain a 'latex_string'")

        try:
            self[key]["default_fiducial_values"] = {}
            for k, v in subdict["default_fiducial_values"].items():
                # use eval so that astropy units are evaluated
                try:
                    self[key]["default_fiducial_values"][k] = eval(v)
                except (SyntaxError, TypeError, ValueError):
                    self[key]["default_fiducial_values"][k] = v
        except KeyError:
            raise KeyError("Equation dictionary must contain 'default_fiducial_values'")

        try:
            self[key]["parts"] = subdict["parts"]
        except KeyError:
            try:
                self[key]["chain"] = subdict["chain"]
            except KeyError:
                raise KeyError("Equation dictionary must contain 'parts' or 'chain'")

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

SUPPLEMENTAL_EQUATIONS = EqDict()

SUPPLEMENTAL_EQUATIONS["rotationfrequency_to_period"] = {
    "description": "The rotation frequency of the pulsar in terms of the rotation period",
    "variable": "rotationfrequency",
    "latex_string": r"f_{\rm rot}",
    "default_fiducial_values": {
        "rotationperiod": 0.01 * u.s,
    },
    "parts": [
        ("rotationperiod", "-1"),
    ],
}

SUPPLEMENTAL_EQUATIONS["gwfrequency_to_rotationfrequency"] = {
    "description": "The GW frequency of the pulsar in terms of the rotation frequency",
    "variable": "gwfrequency",
    "latex_string": r"f_{\rm gw}",
    "default_fiducial_values": {
        "rotationfrequency": 100 * u.Hz,
    },
    "parts": [
        ("2", "1"),
        ("rotationfrequency", "1"),
    ],
}

SUPPLEMENTAL_EQUATIONS["angulargwfrequency_to_gwfrequency"] = {
    "description": "The angular GW frequency of the pulsar in terms of the GW frequency",
    "variable": "angulargwfrequency",
    "latex_string": r"\Omega_{\rm gw}",
    "default_fiducial_values": {
        "gwfrequency": 200 * u.Hz,
    },
    "parts": [
        ("2", "1"),
        ("pi", "1"),
        ("gwfrequency", "1"),
    ],
}

SUPPLEMENTAL_EQUATIONS["angularrotationfrequency_to_angulargwfrequency"] = {
    "description": "The angular rotation frequency of the pulsar in terms of the angular GW frequency",
    "variable": "angularrotationfrequency",
    "latex_string": r"\Omega_{\rm rot}",
    "default_fiducial_values": {
        "angulargwfrequency": 2 * pi * 200 * u.rad / u.s,
    },
    "parts": [
        ("1/2", "1"),
        ("angulargwfrequency", "1"),
    ],
}

SUPPLEMENTAL_EQUATIONS["rotationperiod_to_angularrotationfrequency"] = {
    "description": "The rotation period of the pulsar in terms of the angular rotation frequency",
    "variable": "rotationperiod",
    "latex_string": "P",
    "default_fiducial_values": {
        "angularrotationfrequency": 2 * pi * 100 * u.rad / u.s,
    },
    "parts": [
        ("2", "1"),
        ("pi", "1"),
        ("angularrotationfrequency", "-1"),
    ],
}

SUPPLEMENTAL_EQUATIONS["rotationfdot_to_period"] = {
    "description": "The rotation frequency derivative in terms of period and period derivative",
    "variable": "rotationfdot",
    "latex_string": r"\dot{f}_{\rm rot}",
    "default_fiducial_values": {
        "rotationperiod": 0.01 * u.s,
        "rotationpdot": 1e-15 * u.s / u.s,
    },
    "parts": [
        ("-1", "1"),
        ("rotationpdot", "1"),
        ("rotationperiod", "-2"),
    ],
}

SUPPLEMENTAL_EQUATIONS["rotationpdot_to_angularrotationfdot"] = {
    "description": "The rotation period derivative in terms of angular rotation frequency and its derivative",
    "variable": "rotationpdot",
    "latex_string": r"\dot{P}",
    "default_fiducial_values": {
        "angularrotationfrequency": 2 * pi * 100 * u.rad / u.s,
        "angularrotationfdot": -2 * pi * 1e-11 * u.rad / (u.s ** 2),
    },
    "parts": [
        ("-2", "1"),
        ("pi", "1"),
        ("angularrotationfdot", "1"),
        ("angularrotationfrequency", "-2"),
    ],
}

SUPPLEMENTAL_EQUATIONS["angularrotationfdot_to_angulargwfdot"] = {
    "description": "The angular rotation frequency derivative in terms of angular gravitational-wave frequency derivative",
    "variable": "angularrotationfdot",
    "latex_string": r"\dot{\Omega}_{\rm rot}",
    "default_fiducial_values": {
        "angulargwfdot": 4 * pi * 100 * u.rad / (u.s ** 2),
    },
    "parts": [
        ("1/2", "1"),
        ("angulargwfdot", "1"),
    ],
}

SUPPLEMENTAL_EQUATIONS["angulargwfdot_to_gwfdot"] = {
    "description": "The angular gravitational-wave frequency derivative in terms of gravitational-wave frequency derivative",
    "variable": "angulargwfdot",
    "latex_string": r"\dot{\Omega}_{\rm gw}",
    "default_fiducial_values": {
        "gwfdot": 2 * 100 * u.Hz / u.s,
    },
    "parts": [
        ("2", "1"),
        ("pi", "1"),
        ("gwfdot", "1"),
    ],
}

SUPPLEMENTAL_EQUATIONS["gwfdot_to_rotationfdot"] = {
    "description": "The angular rotation frequency derivative in terms of angular gravitational-wave frequency derivative",
    "variable": "gwfdot",
    "latex_string": r"\dot{f}_{\rm gw}",
    "default_fiducial_values": {
        "rotationfdot": 100 * u.Hz / u.s,
    },
    "parts": [
        ("2", "1"),
        ("rotationfdot", "1"),
    ],
}
