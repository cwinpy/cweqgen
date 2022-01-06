"""
Doc strings for different equation classes.
"""


DOCSTRINGS = {
    "h0": """
Generate the gravitational-wave amplitude for a signal emitted from the l=m=2
mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as
given in :var:`~cweqgen.equations.ALLOWED_VALUES`, can be used instead.

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
""",
    "h0spindown": """
Generate the equation for the spin-down limit on the gravitational wave
amplitude for a signal emitted from the l=m=2 mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as
given in :var:`~cweqgen.equations.ALLOWED_VALUES`, can be used instead.

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
""",
    "ellipticityspindown": """
Generate the equation for the spin-down limit on the source ellipticity
for a signal emitted from the l=m=2 mass quadrupole mode.

For the optional input keyword parameters below a range of aliases, as
given in :var:`cweqgen.equations.ALLOWED_VALUES`, can be used instead.

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
}
