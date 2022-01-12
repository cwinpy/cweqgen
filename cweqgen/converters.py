"""
Functions to convert between equivalent parameters.
"""


def convert_to_rotation_frequency(gwfrequency=None, rotationperiod=None, **kwargs):
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


def convert_to_rotation_fdot(
    gwfrequency=None,
    rotationfrequency=None,
    rotationperiod=None,
    gwfdot=None,
    rotationpdot=None,
    **kwargs,
):
    """
    Convert the GW frequency (assumed to be twice the rotation frequency) or
    the rotation period and GW rotation frequency derivative or rotation
    period derivative into rotation frequency derivative.
    """

    freq = (
        gwfrequency / 2.0
        if gwfrequency is not None
        else (
            (1.0 / rotationperiod) if rotationperiod is not None else rotationfrequency
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
