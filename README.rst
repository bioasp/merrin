Merrin: MEtabolic Regulation Rule INference from time series data
=================================================================

``merrin`` is a *Python3* tool to compute metabolic regulatory rules
from time series observations. This implementation rely on
```merrinasp`` <https://github.com/kthuillier/merrinasp>`__, extension
of the *Answer Set Programming* (ASP) solver
```clingo`` <https://arxiv.org/abs/1705.09811>`__ with quantified linear
constraints.

Quick install
-------------

To install the ``merrin`` package from the GitHub repository, run the
pip command:

.. code:: sh

   python3.X -m pip install git+https://github.com/bioasp/merrin

Usage
-----

``merrin`` can be used in the terminal as follows:

.. code:: sh

   merrin  [-h] -sbml SBML -pkn PKN -obj OBJ -obs OBS [-out OUTPUT] [--lpsolver {glpk,gurobi}] [--timelimit TIMELIMIT]
           [--optimization {all,subsetmin}] [--projection {network,node}]

**Mandatory arguments:**

.. code:: sh

   -sbml SBML, --SBML SBML
                           Metabolic network in SBML file format.
     -pkn PKN, --PKN PKN   Prior Knowledge Network.
     -obj OBJ, --objective-reaction OBJ
                           Objective reaction.
     -obs OBS, --observations OBS
                           JSON file describing the input timeseries.

**Optional arguments:**

.. code:: sh

     -out OUTPUT, --output-file OUTPUT
                           Output CSV file (default: ./merrin-<optimization>-<projection>-<timestamp>.csv)
     --lpsolver {glpk,gurobi}
                           Linear solver to use (default: glpk)
     --timelimit TIMELIMIT
                           Timelimit for each resolution, -1 if none (default: -1)
     --optimization {all,subsetmin}
                           Select optimization mode: all networks or subset minimal ones only (default: subsetmin)
     --projection {network,node}
                           Select project mode (default: network):
                           - node: enumerate the candidate rules for each node;
                           - network: enumerate all the rules of each network

Input files
-----------

Metabolic network
~~~~~~~~~~~~~~~~~

Metabolic network should be in ``SBML`` (*Systems Biology Markup
Language*) version 3 format.

Prior Knowledge Network (PKN)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prior Knowledge Network (PKN) is a text file where each line is such
that:

.. code:: text

   node_1  sign    node_2

with: \* ``node_1`` and ``node_2`` are two components of the regulatory
or metabolic systems. \* ``sign`` in (``0``, ``-1``, ``1``) such that: -
``-1`` is an **inhibition effect** of ``node_1`` on ``node_2``; - ``1``
is an **activation effect** of ``node_1`` on ``node_2``; - ``0`` is an
**unknown effect** (either activation or inhibition effect) of
``node_1`` on ``node_2``;

**Examples**

.. code:: text

   Carbon1 0   RPcl
   RPcl    1   Tc2
   Tc2 -1  RPcl

In this example, ``RPcl`` regulatory rule depends on an unknown
interaction with ``Carbon1`` and an inhibition effect of ``Tc2``.

Observations
~~~~~~~~~~~~

``merrin`` is compatible with any combination of the following
datatypes: *kinetics*, *fluxomics* and *transcriptomics*.

The observations can be noisy. Note that it is preferable not to enter
observations that are not certain.

Observations are described in a ``json`` file. Each time series
observation is defined as follows:

.. code:: json

   {
       "file": "path/to/the/csv/file",
       "type": ["Kinetics","Fluxomics","Transcriptomics"], <- any non-empty subset
       "constraints": {
           "mutations": {
               "node_u": true, <- forced activation
               "node_v": false, <- forced inhibition
           },
           "bounds": {
               "reaction": [lower_bound, upper_bound]
           }
       }
   }

The ``csv`` file describing the observation needs to have a ``Time``
column with an integer timestamp for each observed time step.

**For kinetics and fluxomics data types:** - *Metabolites*: real-values,
modeling the metabolite concentration in the substrate. - Need to
contain a ``biomass`` column with the measured value of the biomass.

**For fluxomics data types:** - *Reaction*: real-values, modeling the
reaction activity rates in the metabolic network.

**For transcriptomics data types:** - All values are binary (``0`` or
``1``), modeling the activity (``1``) or inactivity (``0``) of a
component (metabolite, reaction, regulatory nodes).

Output format
-------------

``merrin`` generates a ``CSV`` file describing the inferred regulatory
networks. A rule set to ``1`` represents a constant value (*i.e.* always
activated) for which no regulatory rules are necessary to explain the
component dynamics.

**Remarks 1:** If *no regulatory networks are returned*, then the
instance is *unsatisfiable*. Try to change the ``max_gap`` and
``max_error`` variables before launching ``merrin`` again.

**Remarks 2:** For unsatisfiable instances with *kinetics* and/or
*fluxomics* data, launching ``merrin`` with the observation declared as
*transcriptomics* data only can sometimes allow inferring some
regulatory networks.

Rule syntax and semantics
-------------------------

| Regulatory rules are returned in **disjunctive normal form** (DNF)
  with the following syntax: > R := 1 \|\| C \|\| (C_1 \| … \| C_n)
| > C := L \|\| (L_1 & … & L_m)
| > L := N \|\| !N
| > N := regulatory component name

with ``!`` denoting the negation, ``&`` the logical and, and ``|`` the
logical or.

Examples
--------

An example is provided in ``./examples``. The instance
``./examples/ecoli-small`` has been generated from the regulatory
metabolic network and the experiments described in `(Covert et al.,
2001) <https://www.sciencedirect.com/science/article/pii/S0022519301924051?via%3Dihub>`__.

To solve the instance using the console command, see the bash file:
``./examples/run-merrin.sh``. It can be executed with:

.. code:: bash

   sh ./examples/run-merrin.sh

To solve the instance using a *Python* script using ``merrin``, check
the *jupyter* notebook: ``./examples/notebook-merrin.ipynb``.

Inferred rules on the example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Network projection:** Infer regulatory networks. Each row of the
output ``CSV`` is a regulatory network and each column is the rules for
a given regulatory component.

**Example 1:** Network projection + All optimization

.. code:: csv

   R2a,R2b,R5a,R5b,R7,R8a,RPO2,RPb,RPcl,RPh,Rres,Tc2
   !RPb,1,1,!RPO2,1,!RPh,!Oxygen,R2b,Carbon1,Hext,1,!RPcl
   !RPb,1,1,!RPO2,!RPb,!RPh,!Oxygen,R2b,Carbon1,Hext,1,!RPcl
   !RPb,1,!RPO2,!RPO2,!RPb,!RPh,!Oxygen,R2b,Carbon1,Hext,1,!RPcl
   !RPb,1,!RPO2,!RPO2,!RPb,!RPh,!Oxygen,R2b,Carbon1,Hext,!RPO2,!RPcl
   ...

Only the first 4 inferred regulatory networks are shown. The node
``R2b`` is always set to ``1``, it does not have any regulatory rules,
and so, is always activated.

**Example 2:** Network projection + Subset minimal optimization

.. code:: csv

   R2a,R2b,R5a,R5b,R7,R8a,RPO2,RPb,RPcl,RPh,Rres,Tc2
   !RPb,1,1,1,1,!RPh,!Oxygen,R2b,Carbon1,Hext,1,!RPcl

**Node projection:** Infer possible regulatory rules for each regulatory
component. Output file will only contain 1 row. Each cell contains a set
of compatible regulatory rules separated by ‘;’.

**Example 3:** Node projection + All optimization

.. code:: csv

   R2a,R2b,R5a,R5b,R7,R8a,RPO2,RPb,RPcl,RPh,Rres,Tc2
   !RPb,1,!RPO2;1,!RPO2;1;RPO2,!RPb;1,!RPh,!Oxygen,R2b,Carbon1,Hext,!RPO2;1,!RPcl

The node ``R5a`` has 2 possible regulatory rules: ``!RPO2`` or ``1``
(unregulated).

**Example 4:** Node projection + Subset minimal optimization

.. code:: csv

   R2a,R2b,R5a,R5b,R7,R8a,RPO2,RPb,RPcl,RPh,Rres,Tc2
   !RPb,1,1,1,1,!RPh,!Oxygen,R2b,Carbon1,Hext,1,!RPcl

References
----------

To cite this tool: > Kerian Thuillier, Caroline Baroukh, Alexander
Bockmayr, Ludovic Cottret, Loïc Paulevé, Anne Siegel, MERRIN: MEtabolic
regulation rule INference from time series data, Bioinformatics, Volume
38, Issue Supplement_2, September 2022, Pages ii127–ii133,
https://doi.org/10.1093/bioinformatics/btac479
[`pdf <https://hal.science/hal-03207589v3/document>`__]
