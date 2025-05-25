# Binding Free Energy Estimator 3
[![DOI](https://zenodo.org/badge/322234705.svg)](https://zenodo.org/badge/latestdoi/322234705)
[![Downloads](https://static.pepy.tech/badge/bfee2)](https://pepy.tech/project/bfee2)

**Binding free energy estimator 3 (BFEE3) has come! There are three major upgrades in BFEE3: (1) LDDM, a high-throughput alchemical route for absolute binding free-energy calculations; (2) streamlined geometrical route for protein-protien binding free-energy calculations and (3) quick settings for common calculations.**

BFEE is a python-based software that automates absolute binding free energy calculations through either the alchemical or geometric route by molecular dynamics simulations.<br>

## Theoretical backgrounds
The degrees of freedom of the protein-ligand (or host-guest) system are described by a series of geometric variables (or collective variables), as firstly described by the [Karplus group](https://pubs.acs.org/doi/abs/10.1021/jp0217839). In BFEE, a generalized, best-fit-rotation-based geometric variables are used, making it in principle available to any protein-ligand complex. See [this paper](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791) for introduction of these variables.<br>

In the [geometric route](https://pubs.acs.org/doi/10.1021/ct3008099), the degrees of freedom is investigated one by one, through one-dimensional free-energy calculations. In BFEE, [WTM-eABF](https://pubs.acs.org/doi/abs/10.1021/acs.accounts.9b00473) is used, while other importance-sampling algorithms such as [plain eABF](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00447) are also acceptable.
The [alchemical route](https://pubs.acs.org/doi/10.1021/ct3008099) is a variants of the [double decoupling method (DDM)](https://www.sciencedirect.com/science/article/pii/S0006349597787563). A thermodynamic cycle, in which the ligand and the geometric restraints are decoupled independently to guarantee the convergence of the simulations.<br>
[这里](http://sioc-journal.cn/Jwk_hxxb/CN/10.6023/A20100489)是标准结合自由能计算方法的中文介绍.<br>

## Features
Generates all the input files for absolute binding free energy calculations;<br>
Support Protein-protein and protein-ligand complexes;<br>
Perform post-treatment automatedly;<br>
Support NAMD (alchemical and geometric route) and Gromacs (geometric route) as molecular dynamics engines;<br>
Support many file formats for the input complex structure (PSF/PDB/PRM, PRM7/RST7, TOP/PDB);<br>
Support both rigid (exclude RMSD CV) and flexible (include RMSD CV) ligands and protein-protein complexes (streamlined geometrical route);<br>
...<br>

## Requirements
Python 3.6+, PySide 2, numpy, scipy, matplotlib, parmed and MDAnalysis.<br>
[NAMD 3.0 or later](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD) / [GROMACS 2024 or later](https://manual.gromacs.org/).<br>
**Note: Since both NAMD and Gromacs have incorporated Colvars in their latest binaries. Each release of BFEE3 corresponds to a version of NAMD/Gromacs. Please always use
the corresponding or later versions of the MD engines for free-energy calculations!**

## Installation
We suggest to install BFEE through conda. It will be safe if conda is install in a new environment<br>
```
conda create --name bfee   (optional)
conda activate bfee        (optional)
conda install -c conda-forge BFEE2
```

## Usage
Simply run BFEE2Gui.py in terminal or PowerShell. One may need to use the absolute path on MS Windows.<br>
A step-by-step tutorial is provided [here](https://www.nature.com/articles/s41596-021-00676-1).<br>
A tutorial about new streamlined geometrical route is provided in the SI of [this paper](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00487).<br>

## Test files
One can download the Supplement Data [here](https://www.nature.com/articles/s41596-021-00676-1#Sec47) for testing BFEE3.

## Citations
When possible, please consider mentioning [Fu et al. Nat. Protoc. 2022, doi:10.1038/s41596-021-00676-1](https://www.nature.com/articles/s41596-021-00676-1#citeas) when BFEE is used in your project.


Additional references:<br>
BFEE2: [Fu et al. J. Chem. Inf. Model. 2021, 61, 2116–2123](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.1c00269)<br>
BFEE2 for protein-protein binding free energy calculations [Fu et al. J. Chem. Inf. Model. 2023, 63, 2512–2519](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00487)<br>
Alchemical and geometric routes [Gumbart et al. J. Chem. Theory Comput. 2013, 9, 794–802](https://pubs.acs.org/doi/abs/10.1021/ct3008099)<br>
WTM-eABF: [Fu et al. Acc. Chem. Res. 2019, 52, 3254–3264](https://pubs.acs.org/doi/abs/10.1021/acs.accounts.9b00473) and [Fu et al. J. Phys. Chem. Lett. 2018, 9, 4738–4745](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01994)<br>
Collective variables: [Fu et al. J. Chem. Theory Comput. 2017, 13, 5173–5178](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791)<br>
Colvars module: [Fiorin et al. Mol. Phys. 2013 111, 3345-3362](https://www.tandfonline.com/doi/full/10.1080/00268976.2013.813594)<br>
"Mother" of all restraint-based binding free-energy calculations: [Hermans et al. Isr. J. Chem. 1986, 27, 225–227](https://onlinelibrary.wiley.com/doi/abs/10.1002/ijch.198600032)<br>

## Contact us
This software is under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) license. For more information about BFEE, contact Haohao Fu (fhh2626@nankai.edu.cn) and Haochuan Chen (yjcoshc@mail.nankai.edu.cn).
