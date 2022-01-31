The CW Equation Generator
=========================

This package, "cweqgen" (pronouced like "*Queck-Jen*"), generates a range of common equations used
by `continuous gravitational-wave <https://www.ligo.org/science/GW-Continuous.php>`_ community. As
well as providing the LaTeX string for the equation, it can evaluate the equations based on
specified values, and create a LaTeX string based on supplied fiducial values.

The current set of equations that the package contains are given in :ref:`Equations`. A selection
of examples of using the package are shown in :ref:`Examples`.

Installation
------------

Installation from PyPI
^^^^^^^^^^^^^^^^^^^^^^

cweqgen is available on `PyPI <https://pypi.org/project/cweqgen/>`_ and can be installed using:

.. code-block:: console

   $ pip install cweqgen

Installation from conda
^^^^^^^^^^^^^^^^^^^^^^^

cweqgen can be install in a `conda <https://docs.conda.io/en/latest/>`_ environment using:

.. code-block:: console

   $ conda install -c conda-forge cweqgen

Installation from source
^^^^^^^^^^^^^^^^^^^^^^^^

cweqgen can be installed from its source `git <https://git-scm.com/>`_ `repository
<https://github.com/cwinpy/cweqgen>`_ using the supplied Python setup script.

First, clone the repository

.. tabbed:: HTTPS

   .. code-block:: console

      $ git clone https://github.com/cwinpy/cweqgen.git

.. tabbed:: ssh

   .. code-block:: console

      $ git clone git@github.com:cwinpy/cweqgen.git

.. tabbed:: GitHub CLI

   Using the GitHub `CLI <https://cli.github.com/>`_

   .. code-block:: console

      $ gh repo clone cwinpy/cweqgen

then install the requirements and the software using:

.. code-block:: console

   $ cd cwinpy/
   $ pip install .


.. toctree::
   :maxdepth: 2
   :hidden:

   Equations <equations>
   examples
   Defining equations <definitions>
   API <api>
