# Binding Free Energy Estimator 3
[![DOI](https://zenodo.org/badge/322234705.svg)](https://zenodo.org/badge/latestdoi/322234705)
[![Downloads](https://static.pepy.tech/badge/bfee2)](https://pepy.tech/project/bfee2)

**Binding Free Energy Estimator 3 (BFEE3) is here! BFEE3 includes many upgrades: <br>
(1) LDDM, a high-throughput alchemical route for absolute binding free-energy calculations; <br>
(2) WTM-λABF, an efficient algorithm for alchemical transformations; <br>
(3) a streamlined geometrical route for protein-protein binding free-energy calculations;<br>
(4) quick setup options for common calculations;<br>
(5) an AI assistant for generating input files and answering questions.**

BFEE is a Python-based software package that automates absolute binding free-energy calculations through either the alchemical or geometric route using molecular dynamics simulations.<br>

## Theoretical background
The degrees of freedom of the protein-ligand or host-guest system are described by a series of geometric variables, or collective variables, as first described by the [Karplus group](https://pubs.acs.org/doi/abs/10.1021/jp0217839). In BFEE, generalized geometric variables based on best-fit rotation are used, making the method, in principle, applicable to any protein-ligand complex. See [this paper](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791) for an introduction to these variables.<br>

In the [geometric route](https://pubs.acs.org/doi/10.1021/ct3008099), the degrees of freedom are investigated one by one through one-dimensional free-energy calculations. In BFEE, [WTM-eABF](https://pubs.acs.org/doi/abs/10.1021/acs.accounts.9b00473) is used, while other importance-sampling algorithms, such as [plain eABF](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00447), are also supported.

The [alchemical route](https://pubs.acs.org/doi/10.1021/ct3008099) is a variant of the [double decoupling method (DDM)](https://www.sciencedirect.com/science/article/pii/S0006349597787563). It uses a thermodynamic cycle in which the ligand and the geometric restraints are decoupled independently to ensure convergence of the simulations.

More recently, high-efficiency methods such as [Lucid DDM](https://www.nature.com/articles/s43588-025-00821-w) and [WTM-λABF](https://pubs.acs.org/doi/10.1021/acs.accounts.5c00666) have also been implemented in BFEE3. These methods are generally recommended as the first choice for absolute binding free-energy calculations.
<br>
[这里](http://sioc-journal.cn/Jwk_hxxb/CN/10.6023/A20100489)是标准结合自由能计算方法的中文介绍。<br>

## Features
Generates all input files for absolute binding free-energy calculations;<br>
Supports protein-protein and protein-ligand complexes;<br>
Performs post-processing automatically;<br>
Supports NAMD for both the alchemical and geometric routes, and GROMACS for the geometric route, as molecular dynamics engines;<br>
Supports many file formats for the input complex structure, including PSF/PDB/PRM, PRM7/RST7, and TOP/PDB;<br>
Supports both rigid ligands and protein-protein complexes, which exclude the RMSD CV, and flexible ligands and protein-protein complexes, which include the RMSD CV, through the streamlined geometrical route;<br>
...<br>

## Requirements
Python 3.6+, PySide 2, numpy, scipy, matplotlib, parmed, and MDAnalysis.<br>
[NAMD 3.0 or later](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD) / [GROMACS 2024 or later](https://manual.gromacs.org/).<br>
**Note: Since both NAMD and GROMACS have incorporated Colvars into their latest binaries, each release of BFEE3 corresponds to a specific version of NAMD/GROMACS. Please always use the corresponding or later versions of the MD engines for free-energy calculations!**

## Installation
We suggest installing BFEE through conda. It is safer to install it in a new conda environment.<br>
```
conda create --name bfee   (optional)
conda activate bfee        (optional)
conda install -c conda-forge BFEE2
```
**IMPORTANT: Please force numpy<2.3 when installing BFEE, because a recent update of numpy has broken ParmEd (https://github.com/ParmEd/ParmEd/issues/1406).**

## Usage
Simply run BFEE2Gui.py in a terminal or PowerShell. On Microsoft Windows, you may need to use the absolute path.<br>
A step-by-step tutorial is provided [here](https://www.nature.com/articles/s41596-021-00676-1).<br>
A tutorial on the new streamlined geometrical route is provided in the Supporting Information of [this paper](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00487).<br>

## Test files
You can download the supplementary data [here](https://www.nature.com/articles/s41596-021-00676-1#Sec47) to test BFEE3.

## Citations
When possible, please consider citing [Fu et al. Nat. Protoc. 2022, doi:10.1038/s41596-021-00676-1](https://www.nature.com/articles/s41596-021-00676-1#citeas) when BFEE is used in your project.

Additional references:<br>
BFEE2: [Fu et al. J. Chem. Inf. Model. 2021, 61, 2116–2123](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.1c00269)<br>
BFEE2 for protein-protein binding free-energy calculations: [Fu et al. J. Chem. Inf. Model. 2023, 63, 2512–2519](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00487)<br>
Alchemical and geometric routes: [Gumbart et al. J. Chem. Theory Comput. 2013, 9, 794–802](https://pubs.acs.org/doi/abs/10.1021/ct3008099)<br>
Lucid DDM method: [Bian et al. Nat. Comput. Sci. 2025, 5, 621–626](https://www.nature.com/articles/s43588-025-00821-w)<br>
WTM-λABF: [Zhou et al. Acc. Chem. Res. 2026, 59, 90–102](https://pubs.acs.org/doi/10.1021/acs.accounts.5c00666)<br>
WTM-eABF: [Fu et al. Acc. Chem. Res. 2019, 52, 3254–3264](https://pubs.acs.org/doi/abs/10.1021/acs.accounts.9b00473) and [Fu et al. J. Phys. Chem. Lett. 2018, 9, 4738–4745](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01994)<br>
Collective variables: [Fu et al. J. Chem. Theory Comput. 2017, 13, 5173–5178](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00791)<br>
Colvars module: [Fiorin et al. Mol. Phys. 2013, 111, 3345–3362](https://www.tandfonline.com/doi/full/10.1080/00268976.2013.813594) and [Fiorin et al. J. Phys. Chem. B 2024, 128, 11108–11123](https://pubs.acs.org/doi/10.1021/acs.jpcb.4c05604)<br>
The "mother" of all restraint-based binding free-energy calculations: [Hermans et al. Isr. J. Chem. 1986, 27, 225–227](https://onlinelibrary.wiley.com/doi/abs/10.1002/ijch.198600032)<br>

## Contact us
This software is distributed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) license. For more information about BFEE, please contact Haohao Fu (fhh2626@nankai.edu.cn) and Haochuan Chen (yjcoshc@mail.nankai.edu.cn).
