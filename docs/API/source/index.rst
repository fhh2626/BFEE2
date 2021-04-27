.. BFEE2 documentation master file, created by
   sphinx-quickstart on Mon Apr 26 23:53:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Binding Free Energy Estimator 2
=================================

Binding free energy estimator (BFEE) is a python-based software that automates absolute binding free energy calculations through either the alchemical or geometric route by molecular dynamics simulations.

Theoretical backgrounds
-----------------------

The degrees of freedom of the protein-ligand (or host-guest) system are described by a series of geometric variables (or collective variables), as firstly described by the `Karplus group`_. In BFEE, a generalized, best-fit-rotation-based geometric variables are used, making it in principle available to any protein-ligand complex. See `this paper`_ for introduction of these variables.

In the `geometric route`_, the degrees of freedom is investigated one by one, through one-dimensional free-energy calculations. In BFEE, `WTM-eABF`_ is used, while other importance-sampling algorithms such as `plain eABF`_ are also acceptable. The `alchemical route`_ is a variants of the double `decoupling method (DDM)`_. A thermodynamic cycle, in which the ligand and the geometric restraints are decoupled independently to guarantee the convergence of the simulations.
`这里`_\ 是标准结合自由能计算方法的中文介绍.

.. External Links

.. _`这里`: http://sioc-journal.cn/Jwk_hxxb/CN/10.6023/A20100489
.. _`Karplus group`: https://pubs.acs.org/doi/abs/10.1021/jp0217839
.. _`this paper`: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791
.. _`geometric route`: https://pubs.acs.org/doi/10.1021/ct3008099
.. _`WTM-eABF`: https://pubs.acs.org/doi/abs/10.1021/acs.accounts.9b00473
.. _`plain eABF`: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00447
.. _`alchemical route`: https://pubs.acs.org/doi/10.1021/ct3008099
.. _`decoupling method (DDM)`: https://www.sciencedirect.com/science/article/pii/S0006349597787563

.. Hide the contents from the front page because they are already in
.. the side bar in the Alabaster sphinx style; requires Alabaster
.. config sidebar_includehidden=True (default)

Features
--------

- Generates all the input files for absolute binding free energy calculations;
- Perform post-treatment automatedly;
- Support NAMD (alchemical and geometric route) and Gromacs (geometric route) as molecular dynamics engines;
- Support many file formats for the input complex structure (PSF/PDB/PRM, PRM7/RST7, TOP/PDB);
- ...

Requirements
------------

- Python 3.6+, PySide 2, numpy, matplotlib and MDAnalysis.
- NAMD 3.0 or later / Colvars patched Gromacs.

Installation
------------

We suggest to install BFEE2 through conda. It will be safe if conda is install in a new environment.

.. code-block:: bash

    conda create --name bfee   (optional)
    conda activate bfee        (optional)
    conda install -c conda-forge BFEE2

Usage
-----

Simply run BFEE2Gui.py in terminal or PowerShell. One may need to use the absolute path on MS Windows.


.. Contents
.. ========

Citations
---------

- BFEE: `Fu et al. J. Chem. Inf. Model. 2018, 58, 556–560`_
- Alchemical and geometric routes: `Gumbart et al. J. Chem. Theory Comput. 2013, 9, 794–802`_
- WTM-eABF: `Fu et al. Acc. Chem. Res. 2019, 52, 3254–3264`_ and `Fu et al. J. Phys. Chem. Lett. 2018, 9, 4738–4745`_
- Collective variables: `Fu et al. J. Chem. Theory Comput. 2017, 13, 5173–5178`_

.. External links

.. _`Fu et al. J. Chem. Inf. Model. 2018, 58, 556–560`: https://pubs.acs.org/doi/10.1021/acs.jcim.7b00695
.. _`Gumbart et al. J. Chem. Theory Comput. 2013, 9, 794–802`: https://pubs.acs.org/doi/abs/10.1021/ct3008099
.. _`Fu et al. Acc. Chem. Res. 2019, 52, 3254–3264`: https://pubs.acs.org/doi/abs/10.1021/acs.accounts.9b00473
.. _`Fu et al. J. Phys. Chem. Lett. 2018, 9, 4738–4745`: https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01994
.. _`Fu et al. J. Chem. Theory Comput. 2017, 13, 5173–5178`: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791

Contact us
----------

Technique issues: Haohao Fu (fhh2626@mail.nankai.edu.cn) and Haochuan Chen (yjcoshc@mail.nankai.edu.cn)

This software is under the `GPLv3`_ license. For more information about the copyright of BFEE, contact the corresponding authors of the aforementioned papers (wscai@nankai.edu.cn, Christophe.Chipot@univ-lorraine.fr).

.. _`GPLv3`: https://www.gnu.org/licenses/gpl-3.0.en.html

.. toctree::
   :maxdepth: 4
   :caption: Contents:
   :hidden:

   python_api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
