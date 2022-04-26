import os
import pkg_resources
import shutil

from ..equations import equations
from ..definitions import EQN_DEFINITIONS


def generate_equations_doc(docfile):
    """
    Helper function to automatically generate a documentation page containing
    all the available equations within cweqgen.

    Parameters
    ----------
    docfile: str:
        The output file for the documentation.
    """

    doccontents = """
#########
Equations
#########

The currently implemented equations are:

"""

    references = """\
References
----------
"""

    usedreferences = []
    usedrefurls = []
    refcount = 1

    for eqn in EQN_DEFINITIONS:
        eqstr = ""

        # create equation
        eq = equations(eqn)

        if eq.reference_string in usedreferences:
            refnum = usedreferences.index(eq.reference_string) + 1
        else:
            usedreferences.append(eq.reference_string)
            usedrefurls.append(eq.reference_adsurl)
            refnum = refcount
            refcount += 1

        eqstr += """\
{0}
{1}
""".format(
            eq.description, "-" * len(eq.description)
        )

        eqstr += """
This equation can be accessed from the :func:`~cweqgen.equations.equations` function using
the name ``{}``.
""".format(
            eqn
        )

        eqno = "" if eq.reference_eqno is None else f"Eqn. {eq.reference_eqno} in "
        eqstr += """
The generated equation ({0} [{1}]_) is:

.. math::

    {2}
""".format(
            eqno, refnum, eq.equation(nocomment=True)
        )

        eqstr += """
The fiducial values defined for this equation are:

.. math::

    {}
""".format(
            eq.fiducial_equation(nocomment=True)
        )

        eqstr += """
.. note::

    These fiducial values are just those defined within this package and may not be representative
    of fiducial values used elsewhere in the literature.
"""

        eqstr += """
To generate the equation as calculated at particular values, the
:func:`~cweqgen.equations.equations` can be used as

.. py:function:: equations("{0}", {1})
    :noindex:
""".format(
            eq.equation_name,
            ", ".join(
                [
                    "{}={}".format(fid, str(val))
                    for fid, val in eq.default_fiducial_values.items()
                ]
            ),
        )

        # add doc string lines
        for line in eq.__doc__.split("\n"):
            eqstr += f"    {line}\n"

        doccontents += eqstr + "\n"

    # add in list of references
    for i in range(len(usedreferences)):
        references += """
.. [{0}] {1} [`ADS URL <{2}>`_]
""".format(
            (i + 1), usedreferences[i], usedrefurls[i]
        )

    with open(docfile, "w") as fp:
        fp.write(doccontents + references)


def generate_yamlexample_doc(docfile, eqn="h0"):
    """
    Output an example YAML file.

    Parameters
    ----------
    docfile: str
        The output file for the documentation.
    eqn: str
        The name of the equation for which the docstring is required.
    """

    src = os.path.join(pkg_resources.resource_filename("cweqgen", "eqnfiles"), f"{eqn}.yaml")
    shutil.copyfile(src, docfile)

