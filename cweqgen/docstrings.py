"""
Doc strings for different equation classes.
"""

from .definitions import EQN_DEFINITIONS


def get_defaults(eqn):
    """
    Get the default values for a given equation.

    Parameters
    ----------
    eqn: str
        The equation name.

    Returns
    -------
    defaults: dict
        A dictionary of default values.
    """

    defaults = {
        key: (item[0].value if hasattr(item[0], "value") else item[0])
        for key, item in EQN_DEFINITIONS[eqn]["default_fiducial_values"].items()
    }

    return defaults


# dictionart of docstrings for equation functions
DOCSTRINGS = {}

# add doctring for h0
DOCSTRINGS[
    "h0"
] = """
Generate the gravitational-wave amplitude for a signal emitted from the l=m=2
mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as given in
:obj:`~cweqgen.equations.ALLOWED_VALUES`, can be used instead.

Parameters
----------
ellipticity: float
    The ellipticity of the source with which the calculate the GW amplitude. If
    given as a float units of kg m^2 are assumed. The default value used if
    none is given is {ellipticity}.
momentofinertia: float, Quantity
    The principle moment of inertia with which the calculate the GW amplitude.
    If given as a float units of kg m^2 are assumed. The default value used if
    none is given is {momentofinertia} kg m^2.
rotationfrequency: float, Quantity
    The rotation frequency of the source. If given as a float units of Hz are
    assumed. The default value if none is given is {rotationfrequency} Hz. If
    the rotational period or gravitational wave frequency are given instead
    then they will be converted into rotational frequency (for GW frequency
    it is assumed that this is twice the rotational frequency).
distance: float, Quantity
    The distance to the source. If given as a float units of kpc are assumed.
    The default value used if none is given is {distance} kpc.
""".format(
    **get_defaults("h0")
)

# add doctring for h0 spin-down limit
DOCSTRINGS[
    "h0spindown"
] = """
Generate the equation for the spin-down limit on the gravitational wave
amplitude for a signal emitted from the l=m=2 mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as given in
:obj:`~cweqgen.equations.ALLOWED_VALUES`, can be used instead.

Parameters
----------
momentofinertia: float, Quantity
    The principle moment of inertia with which the calculate the spin-down
    limit. If given as a float units of kg m^2 are assumed. The default value
    used if none is given is {momentofinertia} kg m^2.
rotationfrequency: float, Quantity
    The rotation frequency of the source. If given as a float units of Hz are
    assumed. The default value if none is given is {rotationfrequency} Hz. If
    the rotational period or gravitational wave frequency are given instead
    then they will be converted into rotational frequency (for GW frequency it
    is assumed that this is twice the rotational frequency).
rotationfdot: float, Quantity
    The first rotational frequency derivative (i.e. the spin-down). If given as
    a float units of Hz/s are assumed. The default value used if none is given
    is {rotationfdot} Hz/s. If the rotational period derivative or
    gravitational wave frequency derivative is given instead then they will be
    converted into rotational frequency derivative.
distance: float, Quantity
    The distance to the source. If given as a float units of kpc are assumed.
    The default value used if none is given is {distance} kpc.
""".format(
    **get_defaults("h0spindown")
)

# add doctring for ellipticity spin-down limit
DOCSTRINGS[
    "ellipticityspindown"
] = """
Generate the equation for the spin-down limit on the source ellipticity
for a signal emitted from the l=m=2 mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as
given in :obj:`cweqgen.equations.ALLOWED_VALUES`, can be used instead.

Parameters
----------
momentofinertia: float, Quantity
    The principle moment of inertia with which the calculate the spin-down
    limit. If given as a float units of kg m^2 are assumed. The default value
    used if none is given is {momentofinertia} kg m^2.
rotationfrequency: float, Quantity
    The rotation frequency of the source. If given as a float units of Hz are
    assumed. The default value if none is given is {rotationfrequency} Hz. If
    the rotational period or gravitational wave frequency are given instead
    then they will be converted into rotational frequency (for GW frequency it
    is assumed that this is twice the rotational frequency).
rotationfdot: float, Quantity
    The first rotational frequency derivative (i.e. the spin-down). If given as
    a float units of Hz/s are assumed. The default value used if none is given
    is {rotationfdot} Hz/s. If the rotational period derivative or
    gravitational wave frequency derivative is given instead then they will be
    converted into rotational frequency derivative.
""".format(
    **get_defaults("ellipticityspindown")
)


# add doctring for braking index
DOCSTRINGS[
    "brakingindex"
] = """
Generate the equation for the braking index of a pulsar.

For the optional input keyword parameters below a range of aliases, as
given in :obj:`cweqgen.equations.ALLOWED_VALUES`, can be used instead.

Parameters
----------
rotationfrequency: float, Quantity
    The rotation frequency of the source. If given as a float units of Hz are
    assumed. The default value if none is given is {rotationfrequency} Hz. If
    the rotational period or gravitational wave frequency are given instead
    then they will be converted into rotational frequency (for GW frequency it
    is assumed that this is twice the rotational frequency).
rotationfdot: float, Quantity
    The first rotational frequency derivative (i.e. the spin-down). If given as
    a float units of Hz/s are assumed. The default value used if none is given
    is {rotationfdot} Hz/s. If the rotational period derivative or
    gravitational wave frequency derivative is given instead then they will be
    converted into rotational frequency derivative.
rotationfddot: float, Quantity
    The second rotational frequency derivative. If given as a float units of
    Hz/s^2 are assumed. The default value used if none is given is
    {rotationfddot} Hz/s^2.
""".format(
    **get_defaults("brakingindex")
)
