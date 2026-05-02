.. BFEE documentation master file.

Binding Free Energy Estimator 3
================================

:Release: |version|
:Date: |today|

Binding Free Energy Estimator (BFEE) is a Python-based software package that
automates absolute binding free-energy calculations through alchemical and
geometrical routes using molecular dynamics simulations.

Theoretical Background
----------------------

The degrees of freedom of the protein-ligand or host-guest system are described
by a series of geometric variables, or collective variables, as first described
by the `Karplus group`_. In BFEE, generalized best-fit-rotation-based geometric
variables are used, making the method in principle applicable to any
protein-ligand complex. See `this paper`_ for an introduction to these
variables.

In the `geometric route`_, the degrees of freedom are investigated one by one
through one-dimensional free-energy calculations. BFEE supports WTM-eABF and
plain eABF for these calculations. The `alchemical route`_ is a variant of the
double `decoupling method (DDM)`_. It uses a thermodynamic cycle in which the
ligand and geometric restraints are decoupled independently to improve
simulation convergence.

Recent BFEE3 releases also include high-efficiency methods such as Lucid DDM
(LDDM) and WTM-lambdaABF, plus a streamlined geometrical route for
protein-protein binding free-energy calculations.

.. External Links

.. _`Karplus group`: https://pubs.acs.org/doi/abs/10.1021/jp0217839
.. _`this paper`: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791
.. _`geometric route`: https://pubs.acs.org/doi/10.1021/ct3008099
.. _`alchemical route`: https://pubs.acs.org/doi/10.1021/ct3008099
.. _`decoupling method (DDM)`: https://www.sciencedirect.com/science/article/pii/S0006349597787563

Features
--------

- Generate input files for absolute binding free-energy calculations.
- Support protein-protein and protein-ligand complexes.
- Support NAMD for alchemical and geometrical routes, and GROMACS for the
  geometrical route.
- Support PSF/PDB/PRM, PRM7/RST7, and TOP/PDB input systems.
- Support DDM, LDDM, WTM-eABF, WTM-lambdaABF, and GaWTM-eABF workflows.
- Provide quick setup actions for common calculation types.
- Provide an OpenAI-compatible AI assistant for setup guidance.
- Perform post-treatment and quick plotting inside the GUI.

Requirements
------------

- Python 3.6+
- PySide6, appdirs, MDAnalysis, matplotlib, numpy<2.3, scipy, parmed, requests
- NAMD 3.0 or later
- GROMACS 2024 or later

Since both NAMD and GROMACS include Colvars in their recent binaries, each BFEE3
release corresponds to a specific tested NAMD/GROMACS version. Use the matching
or a later MD engine version for free-energy calculations.

Installation
------------

We suggest installing BFEE through conda in a new environment.

.. code-block:: bash

    conda create --name bfee   (optional)
    conda activate bfee        (optional)
    conda install -c conda-forge BFEE2

Please force ``numpy<2.3`` if your environment pulls a newer NumPy release,
because recent NumPy versions can break ParmEd compatibility.

Usage
-----

Run ``BFEE2Gui.py`` in a terminal or PowerShell. On Microsoft Windows, you may
need to use the absolute path.

Citations
---------

- BFEE main protocol: `Fu et al. Nat. Protoc. 2022, 17, 1114-1141`_
- BFEE2: `Fu et al. J. Chem. Inf. Model. 2021, 61, 2116-2123`_
- Protein-protein BFEE: `Fu et al. J. Chem. Inf. Model. 2023, 63, 2512-2519`_
- LDDM route: `Bian et al. Nat. Comput. Sci. 2025, 5, 621-626`_
- WTM-lambdaABF: `Zhou et al. Acc. Chem. Res. 2026, 59, 90-102`_
- Alchemical and geometrical routes: `Gumbart et al. J. Chem. Theory Comput. 2013, 9, 794-802`_
- WTM-eABF: `Fu et al. Acc. Chem. Res. 2019, 52, 3254-3264`_ and `Fu et al. J. Phys. Chem. Lett. 2018, 9, 4738-4745`_
- Collective variables: `Fu et al. J. Chem. Theory Comput. 2017, 13, 5173-5178`_

.. External links

.. _`Fu et al. Nat. Protoc. 2022, 17, 1114-1141`: https://www.nature.com/articles/s41596-021-00676-1
.. _`Fu et al. J. Chem. Inf. Model. 2021, 61, 2116-2123`: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.1c00269
.. _`Fu et al. J. Chem. Inf. Model. 2023, 63, 2512-2519`: https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00487
.. _`Bian et al. Nat. Comput. Sci. 2025, 5, 621-626`: https://www.nature.com/articles/s43588-025-00821-w
.. _`Zhou et al. Acc. Chem. Res. 2026, 59, 90-102`: https://pubs.acs.org/doi/10.1021/acs.accounts.5c00666
.. _`Gumbart et al. J. Chem. Theory Comput. 2013, 9, 794-802`: https://pubs.acs.org/doi/abs/10.1021/ct3008099
.. _`Fu et al. Acc. Chem. Res. 2019, 52, 3254-3264`: https://pubs.acs.org/doi/abs/10.1021/acs.accounts.9b00473
.. _`Fu et al. J. Phys. Chem. Lett. 2018, 9, 4738-4745`: https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01994
.. _`Fu et al. J. Chem. Theory Comput. 2017, 13, 5173-5178`: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791

Contact Us
----------

This software is distributed under the `GPLv3`_ license. For more information
about BFEE, contact Haohao Fu (fhh2626@nankai.edu.cn) and Haochuan Chen
(yjcoshc@mail.nankai.edu.cn).

.. _`GPLv3`: https://www.gnu.org/licenses/gpl-3.0.en.html

.. toctree::
   :maxdepth: 4
   :caption: Contents:
   :hidden:

   python_api

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
