"""
The equation definition information is held in the :obj:`~cweqgen.definitions.EQN_DEFINITIONS`
dictionary. If defining a new equation the information should be added into this dictionary. The
dictionary keys give the equation's "name", with which it must be referred to. The values are
dictionaries containing the following keys:

"default_fiducial_values":
   a dictionary with keys for all variable parameters in the equation. Each key gives a
   2-tuple containing the default value of that variable (with appropriate
   :class:`astropy.units.Unit`) and the variable's exponent as a string. The variable parameter names
   must be consistent with the names in the :obj:`~cweqgen.equations.ALLOWED_VALUES` dictionary.

"equation_constants":
   a list of 2-tuples containing the constants in the equation and their
   exponents. These should both be string values.

"additional_values":
   a list of variable names that can be used to derive a subset of the
   variables in the "default_fiducial_values", i.e., if the equation requires the "rotationfrequency"
   then this could contain "rotationperiod", which can instead be used to derive the rotation
   frequency. The variable parameter names must be consistent with the names in the
   :obj:`~cweqgen.equations.ALLOWED_VALUES` dictionary.

"converters":
   a dictionary of coversion functions to convert the additional values into the
   required values.
"""

import astropy.units as u

from .converters import *

#: equation definitions
EQN_DEFINITIONS = {
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
    "brakingindex": {
        "default_fiducial_values": {
            "rotationfrequency": (50 * u.Hz, "1"),
            "rotationfddot": (1e-23 * u.Hz / (u.s ** 2), "1"),
            "rotationfdot": (-1e-11 * u.Hz / u.s, "-2"),
        },
        "constants": [],
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
    }
}