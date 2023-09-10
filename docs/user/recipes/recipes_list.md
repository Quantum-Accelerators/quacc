# List of Recipes

## DFTB+

!!! Note

    [DFTB+](https://dftbplus.org/) is especially useful for periodic GFN-xTB calculations and the DFTB+ method based on Slater-Koster parameters.

<center>

| Name         | Decorator       | Documentation                          |
| ------------ | --------------- | -------------------------------------- |
| DFTB+ Static | `#!Python @job` | [quacc.recipes.dftb.core.static_job][] |
| DFTB+ Relax  | `#!Python @job` | [quacc.recipes.dftb.core.relax_job][]  |

</center>

## EMT

!!! Note

    [Effective medium theory (EMT)](https://doi.org/10.1016/0039-6028(96)00816-3) is a semi-empirical method for modeling solids that is predominantly used for prototyping workflows. Because it is solely for demonstration purposes, it only supports the following metals: Al, Ni, Cu, Pd, Ag, Pt, and Au.

<center>

| Name                | Decorator        | Documentation                                      |
| ------------------- | ---------------- | -------------------------------------------------- |
| EMT Static          | `#!Python @job`  | [quacc.recipes.emt.core.static_job][]              |
| EMT Relax           | `#!Python @job`  | [quacc.recipes.emt.core.relax_job][]               |
| EMT Bulk to Defects | `#!Python @flow` | [quacc.recipes.emt.defects.bulk_to_defects_flow][] |
| EMT Bulk to Slabs   | `#!Python @flow` | [quacc.recipes.emt.slabs.bulk_to_slabs_flow][]     |

</center>

## Gaussian

!!! Note

    [Gaussian](https://gaussian.com/) is an extremely popular molecular DFT code that is quite robust and easy to use.

<center>

| Name            | Decorator       | Documentation                              |
| --------------- | --------------- | ------------------------------------------ |
| Gaussian Static | `#!Python @job` | [quacc.recipes.gaussian.core.static_job][] |
| Gaussian Relax  | `#!Python @job` | [quacc.recipes.gaussian.core.relax_job][]  |

</center>

## GULP

!!! Note

    [GULP](https://gulp.curtin.edu.au/) is especially useful for periodic GFN-FF calculations and force field methods. GULP can be downloaded and installed [here](https://gulp.curtin.edu.au/download.html).

<center>

| Name        | Decorator       | Documentation                          |
| ----------- | --------------- | -------------------------------------- |
| GULP Static | `#!Python @job` | [quacc.recipes.gulp.core.static_job][] |
| GULP Relax  | `#!Python @job` | [quacc.recipes.gulp.core.relax_job][]  |

</center>

## Lennard-Jones Potential

!!! Note

    [Lennard Jones (LJ)](https://en.wikipedia.org/wiki/Lennard-Jones_potential) is an empirical potential that is predominantly used for prototyping workflows for molecules.

<center>

| Name         | Decorator       | Documentation                        |
| ------------ | --------------- | ------------------------------------ |
| LJ Static    | `#!Python @job` | [quacc.recipes.lj.core.static_job][] |
| LJ Relax     | `#!Python @job` | [quacc.recipes.lj.core.relax_job][]  |
| LJ Frequency | `#!Python @job` | [quacc.recipes.lj.core.freq_job][]   |

</center>

## NewtonNet

!!! Note

    [NewtonNet](https://github.com/ericyuan00000/NewtonNet) is a message passing network for deep learning of interatomic potentials and forces, as described [here](https://pubs.rsc.org/en/content/articlehtml/2022/dd/d2dd00008c).

<center>

| Name                | Decorator       | Documentation                                |
| ------------------- | --------------- | -------------------------------------------- |
| NewtonNet Static    | `#!Python @job` | [quacc.recipes.newtonnet.core.static_job][]  |
| NewtonNet Relax     | `#!Python @job` | [quacc.recipes.newtonnet.core.relax_job][]   |
| NewtonNet Frequency | `#!Python @job` | [quacc.recipes.newtonnet.core.freq_job][]    |
| NewtonNet Frequency | `#!Python @job` | [quacc.recipes.newtonnet.core.freq_job][]    |
| NewtonNet TS        | `#!Python @job` | [quacc.recipes.newtonnet.ts.ts_job][]        |
| NewtonNet IRC       | `#!Python @job` | [quacc.recipes.newtonnet.ts.irc_job][]       |
| NewtonNet Quasi IRC | `#!Python @job` | [quacc.recipes.newtonnet.ts.quasi_irc_job][] |

</center>

## ORCA

!!! Note

    [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal) is a free code that is especially useful for molecular DFT calculations with recently developed methods. ORCA can be downloaded and installed [here](https://orcaforum.kofo.mpg.de/app.php/dlext/).

<center>

| Name        | Decorator       | Documentation                          |
| ----------- | --------------- | -------------------------------------- |
| ORCA Static | `#!Python @job` | [quacc.recipes.orca.core.static_job][] |
| ORCA Relax  | `#!Python @job` | [quacc.recipes.orca.core.relax_job][]  |

</center>

## Psi4

!!! Note

    [Psi4](https://github.com/psi4/psi4) is an open-source quantum chemistry electronic structure package.

<center>

| Name        | Decorator       | Documentation                          |
| ----------- | --------------- | -------------------------------------- |
| Psi4 Static | `#!Python @job` | [quacc.recipes.psi4.core.static_job][] |

</center>

## Q-Chem

!!! Note

    [Q-Chem](https://www.q-chem.com/) is a powerful, general-purpose molecular DFT code with a variety of features.

<center>

| Name             | Decorator       | Documentation                            |
| ---------------- | --------------- | ---------------------------------------- |
| Q-Chem Static    | `#!Python @job` | [quacc.recipes.qchem.core.static_job][]  |
| Q-Chem Relax     | `#!Python @job` | [quacc.recipes.qchem.core.relax_job][]   |
| Q-Chem TS        | `#!Python @job` | [quacc.recipes.qchem.ts.ts_job][]        |
| Q-Chem IRC       | `#!Python @job` | [quacc.recipes.qchem.ts.irc_job][]       |
| Q-Chem Quasi IRC | `#!Python @job` | [quacc.recipes.qchem.ts.quasi_irc_job][] |

</center>

## TBLite

!!! Note

    [tblite](https://github.com/tblite/tblite) is a code that interfaces with the xtb package for running GFN-xTB calculations.

<center>

| Name             | Decorator       | Documentation                            |
| ---------------- | --------------- | ---------------------------------------- |
| TBLite Static    | `#!Python @job` | [quacc.recipes.tblite.core.static_job][] |
| TBLite Relax     | `#!Python @job` | [quacc.recipes.tblite.core.relax_job]    |
| TBLite Frequency | `#!Python @job` | [quacc.recipes.tblite.core.freq_job][]   |

</center>

## VASP

!!! Note

    [VASP](https://www.vasp.at/) is a very widely used code for plane-wave, periodic DFT calculations. Quacc has built-in support for automatically fixing failed VASP jobs via [Custodian](https://github.com/materialsproject/custodian).

<center>

| Name                    | Decorator        | Documentation                                   |
| ----------------------- | ---------------- | ----------------------------------------------- |
| VASP Static             | `#!Python @job`  | [quacc.recipes.vasp.core.static_job][]          |
| VASP Relax              | `#!Python @job`  | [quacc.recipes.vasp.core.relax_job][]           |
| VASP Double Relax       | `#!Python @job`  | [quacc.recipes.vasp.core.double_relax_job][]    |
| VASP Slab Static        | `#!Python @job`  | [quacc.recipes.vasp.slabs.slab_static_job][]    |
| VASP Slab Relax         | `#!Python @job`  | [quacc.recipes.vasp.slabs.slab_relax_job][]     |
| VASP Bulk to Slabs      | `#!Python @flow` | [quacc.recipes.vasp.slabs.bulk_to_slabs_flow][] |
| VASP Slab to Adsorbates | `#!Python @flow` | [quacc.recipes.vasp.slabs.slab_to_ads_flow][]   |
| VASP MP-Prerelax        | `#!Python @job`  | [quacc.recipes.vasp.mp.mp_relax_job][]          |
| VASP MP-Relax           | `#!Python @job`  | [quacc.recipes.vasp.mp.mp_relax_job][]          |
| VASP MP Workflow        | `#!Python @flow` | [quacc.recipes.vasp.mp.mp_relax_flow][]         |
| VASP QMOF Relax         | `#!Python @job`  | [quacc.recipes.vasp.qmof.qmof_relax_job][]      |

</center>
