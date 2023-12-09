from __future__ import annotations

from io import StringIO
from pathlib import Path
import numpy as np

from ase.io.espresso import construct_namelist, read_fortran_namelist

from quacc.calculators.espresso.io import (
    write_espresso_io,
    write_espresso_ph,
    read_espresso_ph,
)
from quacc.calculators.espresso.keys import ALL_KEYS


def test_write_espresso_io():
    input_data = {
        "degauss": 2,
        "Emin": True,
        "Emax": None,
        "nbnd": 10,
        "ecutwfc": 20,
        "irmin(2, 6)": 0.1,
        "outdir(5)": "/path/to/outdir",
    }

    input_data = construct_namelist(input_data, ALL_KEYS["projwfc"])

    string_io = StringIO()

    write_espresso_io(string_io, input_data=input_data)

    expected = (
        "&PROJWFC\n"
        "   outdir(5)        = '/path/to/outdir'\n"
        "   degauss          = 2\n"
        "   emin             = .true.\n"
        "   irmin(2, 6)      = 0.1\n"
        "/\n"
        "EOF"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


def test_write_espresso_ph_single():
    input_data = {
        "amass(1)": 1.0,
        "amass(2)": 2.0,
        "prefix": "prefix",
        "outdir": "/path/to/outdir",
        "eth_rps": 0.1,
    }

    qpts = (0.5, -0.1, 1 / 3)

    input_data = construct_namelist(input_data, ALL_KEYS["ph"])

    string_io = StringIO()

    write_espresso_ph(string_io, input_data=input_data, qpts=qpts)

    expected = (
        "&INPUTPH\n"
        "   amass(1)         = 1.0\n"
        "   amass(2)         = 2.0\n"
        "   outdir           = '/path/to/outdir'\n"
        "   prefix           = 'prefix'\n"
        "   eth_rps          = 0.1\n"
        "/\n"
        "0.50000000 -0.10000000 0.33333333\n"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


def test_write_espresso_ph_list():
    input_data = {
        "amass(1)": 1.0,
        "amass(2)": 2.0,
        "prefix": "prefix",
        "outdir": "/path/to/outdir",
        "eth_rps": 0.1,
        "qplot": True,
        "ldisp": True,
    }

    qpts = [(0.5, -0.1, 1 / 3, 2), (0.1, 0.2, 0.3, 10), (0.2, 0.3, 0.4, 1)]

    input_data = construct_namelist(input_data, ALL_KEYS["ph"])

    string_io = StringIO()

    write_espresso_ph(string_io, input_data=input_data, qpts=qpts)

    expected = (
        "&INPUTPH\n"
        "   amass(1)         = 1.0\n"
        "   amass(2)         = 2.0\n"
        "   outdir           = '/path/to/outdir'\n"
        "   prefix           = 'prefix'\n"
        "   eth_rps          = 0.1\n"
        "   qplot            = .true.\n"
        "   ldisp            = .true.\n"
        "/\n"
        "3\n"
        "0.50000000 -0.10000000 0.33333333 2\n"
        "0.10000000 0.20000000 0.30000000 10\n"
        "0.20000000 0.30000000 0.40000000 1\n"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


def test_write_espresso_ph_nat_todo():
    input_data = {
        "amass(1)": 1.0,
        "amass(2)": 2.0,
        "prefix": "prefix",
        "outdir": "/path/to/outdir",
        "eth_rps": 0.1,
        "qplot": True,
        "nat_todo": True,
        "ldisp": True,
    }

    qpts = [(0.5, -0.1, 1 / 3, 1), (0.1, 0.2, 0.3, -1), (0.2, 0.3, 0.4, 4)]

    input_data = construct_namelist(input_data, ALL_KEYS["ph"])

    string_io = StringIO()

    write_espresso_ph(string_io, input_data=input_data, qpts=qpts, nat_todo=[1, 2, 3])

    expected = (
        "&INPUTPH\n"
        "   amass(1)         = 1.0\n"
        "   amass(2)         = 2.0\n"
        "   outdir           = '/path/to/outdir'\n"
        "   prefix           = 'prefix'\n"
        "   eth_rps          = 0.1\n"
        "   qplot            = .true.\n"
        "   ldisp            = .true.\n"
        "   nat_todo         = .true.\n"
        "/\n"
        "3\n"
        "0.50000000 -0.10000000 0.33333333 1\n"
        "0.10000000 0.20000000 0.30000000 -1\n"
        "0.20000000 0.30000000 0.40000000 4\n"
        "1 2 3\n"
    )

    string_io.seek(0)

    recycled_input_data = read_fortran_namelist(string_io)[0]

    assert recycled_input_data == input_data
    assert string_io.getvalue() == expected


test_1 = """
     Program PHONON v.6.0 (svn rev. 13188M) starts on  7Dec2016 at 10:43: 7 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     R & G space division:  proc/nbgrp/npool/nimage =       4

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/aluminum.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     10                  216      216      45
     Max          31      31     11                  218      218      46
     Sum         121     121     43                  869      869     181



     Dynamical matrices for ( 4, 4, 4)  uniform grid of q-points
     (   8q-points):
       N         xq(1)         xq(2)         xq(3) 
       1   0.000000000   0.000000000   0.000000000
       2  -0.250000000   0.250000000  -0.250000000
       3   0.500000000  -0.500000000   0.500000000
       4   0.000000000   0.500000000   0.000000000
       5   0.750000000  -0.250000000   0.750000000
       6   0.500000000   0.000000000   0.500000000
       7   0.000000000  -1.000000000   0.000000000
       8  -0.500000000  -1.000000000   0.000000000

     Calculation of q =    0.0000000   0.0000000   0.0000000

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   0.0000000 )

     49 Sym.Ops. (with q -> -q+G )


     G cutoff =   85.4897  (    217 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=    29  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, O_h (m-3m)  point group:


     Atomic displacements:
     There are   1 irreducible representations

     Representation     1      3 modes -T_1u G_15  G_4-  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     0.17s CPU         0.19s WALL



     Representation #  1 modes #   1  2  3

     Self-consistent Calculation

     Pert. #  1: Fermi energy shift (Ry) =     5.5145E-25    -2.5077E-37
     Pert. #  2: Fermi energy shift (Ry) =    -2.3437E-24     3.6048E-37
     Pert. #  3: Fermi energy shift (Ry) =    -1.3097E-24     3.1347E-38

      iter #   1 total cpu time :     0.3 secs   av.it.:   3.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.257E-08

     Pert. #  1: Fermi energy shift (Ry) =    -3.3087E-24     1.3469E-39
     Pert. #  2: Fermi energy shift (Ry) =    -2.7573E-25     6.7346E-40
     Pert. #  3: Fermi energy shift (Ry) =     3.2398E-24     6.1224E-40

      iter #   2 total cpu time :     0.4 secs   av.it.:   5.5
      thresh= 1.121E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.314E-09

     Pert. #  1: Fermi energy shift (Ry) =    -1.3786E-24    -1.6224E-39
     Pert. #  2: Fermi energy shift (Ry) =     6.8932E-25     1.1020E-39
     Pert. #  3: Fermi energy shift (Ry) =     2.0680E-25    -4.2857E-40

      iter #   3 total cpu time :     0.5 secs   av.it.:   5.3
      thresh= 3.625E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.570E-13

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       0.173268 [THz] =       5.779601 [cm-1]
     freq (    2) =       0.173268 [THz] =       5.779601 [cm-1]
     freq (    3) =       0.173268 [THz] =       5.779601 [cm-1]
 **************************************************************************

     Mode symmetry, O_h (m-3m)  point group:

     freq (  1 -  3) =          5.8  [cm-1]   --> T_1u G_15  G_4- I  
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0000   gamma=    0.00 GHz
     lambda(    2)=  0.0000   gamma=    0.00 GHz
     lambda(    3)=  0.0000   gamma=    0.00 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0000   gamma=    0.00 GHz
     lambda(    2)=  0.0000   gamma=    0.00 GHz
     lambda(    3)=  0.0000   gamma=    0.00 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.0000   gamma=    0.00 GHz
     lambda(    2)=  0.0000   gamma=    0.00 GHz
     lambda(    3)=  0.0000   gamma=    0.00 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.0000   gamma=    0.02 GHz
     lambda(    2)=  0.0000   gamma=    0.03 GHz
     lambda(    3)=  0.0000   gamma=    0.03 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.0000   gamma=    0.08 GHz
     lambda(    2)=  0.0000   gamma=    0.09 GHz
     lambda(    3)=  0.0000   gamma=    0.09 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.0000   gamma=    0.16 GHz
     lambda(    2)=  0.0000   gamma=    0.18 GHz
     lambda(    3)=  0.0000   gamma=    0.18 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.0000   gamma=    0.25 GHz
     lambda(    2)=  0.0000   gamma=    0.27 GHz
     lambda(    3)=  0.0000   gamma=    0.27 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.0000   gamma=    0.35 GHz
     lambda(    2)=  0.0000   gamma=    0.38 GHz
     lambda(    3)=  0.0000   gamma=    0.38 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.0000   gamma=    0.48 GHz
     lambda(    2)=  0.0000   gamma=    0.50 GHz
     lambda(    3)=  0.0000   gamma=    0.50 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.0000   gamma=    0.61 GHz
     lambda(    2)=  0.0000   gamma=    0.63 GHz
     lambda(    3)=  0.0000   gamma=    0.64 GHz


     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

     Calculation of q =   -0.2500000   0.2500000  -0.2500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     13                  216      216      64
     Max          31      31     14                  218      218      65
     Sum         121     121     55                  869      869     259


     Title: 
     Electron-phonon coefficients for Al                                        


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   240  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.12Mb

     Estimated total allocated dynamical RAM >       0.47Mb

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_2/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.7

     total cpu time spent up to now is        1.2 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.2500000   0.2500000  -0.2500000 )

      6 Sym.Ops. (no q -> -q+G )


     G cutoff =   85.4897  (    218 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   240  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      2 modes -E    L_3  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     3.67s CPU         3.94s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     4.0 secs   av.it.:   4.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.094E-02

      iter #   2 total cpu time :     4.1 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  9.107E-01

      iter #   3 total cpu time :     4.2 secs   av.it.:   4.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.162E-07

      iter #   4 total cpu time :     4.3 secs   av.it.:   5.2
      thresh= 7.185E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.353E-09

      iter #   5 total cpu time :     4.4 secs   av.it.:   5.4
      thresh= 4.851E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.600E-10

      iter #   6 total cpu time :     4.5 secs   av.it.:   5.2
      thresh= 1.265E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.187E-11

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   2  3

     Self-consistent Calculation

      iter #   1 total cpu time :     4.7 secs   av.it.:   3.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.275E-08

      iter #   2 total cpu time :     4.9 secs   av.it.:   6.0
      thresh= 1.810E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.070E-09

      iter #   3 total cpu time :     5.1 secs   av.it.:   5.7
      thresh= 5.541E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.011E-11

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    8
     List of q in the star:
          1  -0.250000000   0.250000000  -0.250000000
          2   0.250000000  -0.250000000  -0.250000000
          3   0.250000000  -0.250000000   0.250000000
          4   0.250000000   0.250000000   0.250000000
          5  -0.250000000  -0.250000000  -0.250000000
          6  -0.250000000  -0.250000000   0.250000000
          7  -0.250000000   0.250000000   0.250000000
          8   0.250000000   0.250000000  -0.250000000

     Diagonalizing the dynamical matrix

     q = (   -0.250000000   0.250000000  -0.250000000 ) 

 **************************************************************************
     freq (    1) =       3.512771 [THz] =     117.173427 [cm-1]
     freq (    2) =       3.512771 [THz] =     117.173427 [cm-1]
     freq (    3) =       6.338040 [THz] =     211.414258 [cm-1]
 **************************************************************************

     Mode symmetry, C_3v (3m)   point group:

     freq (  1 -  2) =        117.2  [cm-1]   --> E    L_3           
     freq (  3 -  3) =        211.4  [cm-1]   --> A_1  L_1           
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0023   gamma=    0.04 GHz
     lambda(    2)=  0.0023   gamma=    0.04 GHz
     lambda(    3)=  0.0285   gamma=    1.47 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0204   gamma=    0.45 GHz
     lambda(    2)=  0.0207   gamma=    0.46 GHz
     lambda(    3)=  0.2321   gamma=   16.75 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.0250   gamma=    0.63 GHz
     lambda(    2)=  0.0251   gamma=    0.63 GHz
     lambda(    3)=  0.2280   gamma=   18.57 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.0283   gamma=    0.75 GHz
     lambda(    2)=  0.0282   gamma=    0.75 GHz
     lambda(    3)=  0.2027   gamma=   17.50 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.0323   gamma=    0.89 GHz
     lambda(    2)=  0.0322   gamma=    0.88 GHz
     lambda(    3)=  0.1880   gamma=   16.81 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.0366   gamma=    1.03 GHz
     lambda(    2)=  0.0365   gamma=    1.03 GHz
     lambda(    3)=  0.1841   gamma=   16.92 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.0408   gamma=    1.18 GHz
     lambda(    2)=  0.0408   gamma=    1.18 GHz
     lambda(    3)=  0.1873   gamma=   17.64 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.0448   gamma=    1.33 GHz
     lambda(    2)=  0.0449   gamma=    1.33 GHz
     lambda(    3)=  0.1946   gamma=   18.72 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.0485   gamma=    1.46 GHz
     lambda(    2)=  0.0485   gamma=    1.46 GHz
     lambda(    3)=  0.2039   gamma=   19.97 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.0517   gamma=    1.58 GHz
     lambda(    2)=  0.0516   gamma=    1.57 GHz
     lambda(    3)=  0.2137   gamma=   21.23 GHz


     Number of q in the star =    8
     List of q in the star:
          1  -0.250000000   0.250000000  -0.250000000
          2   0.250000000  -0.250000000  -0.250000000
          3   0.250000000  -0.250000000   0.250000000
          4   0.250000000   0.250000000   0.250000000
          5  -0.250000000  -0.250000000  -0.250000000
          6  -0.250000000  -0.250000000   0.250000000
          7  -0.250000000   0.250000000   0.250000000
          8   0.250000000   0.250000000  -0.250000000

     Calculation of q =    0.5000000  -0.5000000   0.5000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     15                  216      216      82
     Max          31      31     16                  218      218      83
     Sum         121     121     61                  869      869     331


     Title: 
     Electron-phonon coefficients for Al                                        


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   130  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.12Mb

     Estimated total allocated dynamical RAM >       0.47Mb

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_3/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.8

     total cpu time spent up to now is        2.3 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.5000000  -0.5000000   0.5000000 )

     13 Sym.Ops. (with q -> -q+G )


     G cutoff =   85.4897  (    218 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   130  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, D_3d (-3m)  point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_2u L_2'  To be done

     Representation     2      2 modes -E_u  L_3'  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     7.08s CPU         7.74s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     7.8 secs   av.it.:   4.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.394E-04

      iter #   2 total cpu time :     7.9 secs   av.it.:   5.5
      thresh= 1.547E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.813E-04

      iter #   3 total cpu time :     7.9 secs   av.it.:   5.0
      thresh= 1.677E-03 alpha_mix =  0.700 |ddv_scf|^2 =  6.318E-09

      iter #   4 total cpu time :     8.0 secs   av.it.:   5.5
      thresh= 7.949E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.940E-10

      iter #   5 total cpu time :     8.0 secs   av.it.:   5.1
      thresh= 1.715E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.672E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   2  3

     Self-consistent Calculation

      iter #   1 total cpu time :     8.2 secs   av.it.:   3.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.601E-08

      iter #   2 total cpu time :     8.3 secs   av.it.:   5.9
      thresh= 1.898E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.165E-09

      iter #   3 total cpu time :     8.4 secs   av.it.:   5.5
      thresh= 5.626E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.781E-11

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    4
     List of q in the star:
          1   0.500000000  -0.500000000   0.500000000
          2   0.500000000   0.500000000   0.500000000
          3  -0.500000000   0.500000000   0.500000000
          4   0.500000000   0.500000000  -0.500000000

     Diagonalizing the dynamical matrix

     q = (    0.500000000  -0.500000000   0.500000000 ) 

 **************************************************************************
     freq (    1) =       4.438882 [THz] =     148.065163 [cm-1]
     freq (    2) =       4.438882 [THz] =     148.065163 [cm-1]
     freq (    3) =       9.422553 [THz] =     314.302524 [cm-1]
 **************************************************************************

     Mode symmetry, D_3d (-3m)  point group:

     freq (  1 -  2) =        148.1  [cm-1]   --> E_u  L_3'          
     freq (  3 -  3) =        314.3  [cm-1]   --> A_2u L_2'          
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0000   gamma=    0.00 GHz
     lambda(    2)=  0.0000   gamma=    0.00 GHz
     lambda(    3)=  0.0000   gamma=    0.00 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0290   gamma=    1.03 GHz
     lambda(    2)=  0.0262   gamma=    0.93 GHz
     lambda(    3)=  0.0411   gamma=    6.56 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.0668   gamma=    2.67 GHz
     lambda(    2)=  0.0610   gamma=    2.44 GHz
     lambda(    3)=  0.1030   gamma=   18.53 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.0779   gamma=    3.30 GHz
     lambda(    2)=  0.0722   gamma=    3.06 GHz
     lambda(    3)=  0.1296   gamma=   24.73 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.0826   gamma=    3.62 GHz
     lambda(    2)=  0.0782   gamma=    3.43 GHz
     lambda(    3)=  0.1457   gamma=   28.78 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.0853   gamma=    3.84 GHz
     lambda(    2)=  0.0824   gamma=    3.71 GHz
     lambda(    3)=  0.1554   gamma=   31.56 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.0869   gamma=    4.01 GHz
     lambda(    2)=  0.0853   gamma=    3.94 GHz
     lambda(    3)=  0.1609   gamma=   33.49 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.0881   gamma=    4.16 GHz
     lambda(    2)=  0.0876   gamma=    4.13 GHz
     lambda(    3)=  0.1647   gamma=   35.02 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.0891   gamma=    4.28 GHz
     lambda(    2)=  0.0894   gamma=    4.29 GHz
     lambda(    3)=  0.1676   gamma=   36.28 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.0897   gamma=    4.37 GHz
     lambda(    2)=  0.0906   gamma=    4.42 GHz
     lambda(    3)=  0.1699   gamma=   37.30 GHz


     Number of q in the star =    4
     List of q in the star:
          1   0.500000000  -0.500000000   0.500000000
          2   0.500000000   0.500000000   0.500000000
          3  -0.500000000   0.500000000   0.500000000
          4   0.500000000   0.500000000  -0.500000000

     Calculation of q =    0.0000000   0.5000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     13                  216      216      64
     Max          31      31     14                  218      218      65
     Sum         121     121     55                  869      869     259


     Title: 
     Electron-phonon coefficients for Al                                        


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   200  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.12Mb

     Estimated total allocated dynamical RAM >       0.47Mb

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_4/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.4

     total cpu time spent up to now is        3.5 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.5000000   0.0000000 )

      8 Sym.Ops. (no q -> -q+G )


     G cutoff =   85.4897  (    218 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   200  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_4v (4mm)  point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_1  G_1 D_1  To be done

     Representation     2      2 modes -E    G_5 D_5  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :    10.34s CPU        11.36s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :    11.4 secs   av.it.:   3.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.373E-03

      iter #   2 total cpu time :    11.5 secs   av.it.:   4.5
      thresh= 9.151E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.998E-01

      iter #   3 total cpu time :    11.6 secs   av.it.:   4.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.925E-08

      iter #   4 total cpu time :    11.7 secs   av.it.:   5.5
      thresh= 2.434E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.217E-09

      iter #   5 total cpu time :    11.7 secs   av.it.:   5.0
      thresh= 4.709E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.721E-10

      iter #   6 total cpu time :    11.8 secs   av.it.:   4.3
      thresh= 1.312E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.106E-12

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   2  3

     Self-consistent Calculation

      iter #   1 total cpu time :    12.0 secs   av.it.:   3.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.929E-08

      iter #   2 total cpu time :    12.2 secs   av.it.:   6.1
      thresh= 2.988E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.127E-09

      iter #   3 total cpu time :    12.4 secs   av.it.:   5.6
      thresh= 5.592E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.752E-10

      iter #   4 total cpu time :    12.5 secs   av.it.:   5.4
      thresh= 1.324E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.767E-14

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    6
     List of q in the star:
          1   0.000000000   0.500000000   0.000000000
          2   0.000000000  -0.500000000   0.000000000
          3   0.500000000   0.000000000   0.000000000
          4   0.000000000   0.000000000   0.500000000
          5   0.000000000   0.000000000  -0.500000000
          6  -0.500000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.500000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       4.200435 [THz] =     140.111422 [cm-1]
     freq (    2) =       4.200435 [THz] =     140.111422 [cm-1]
     freq (    3) =       6.478556 [THz] =     216.101363 [cm-1]
 **************************************************************************

     Mode symmetry, C_4v (4mm)  point group:

     freq (  1 -  2) =        140.1  [cm-1]   --> E    G_5 D_5       
     freq (  3 -  3) =        216.1  [cm-1]   --> A_1  G_1 D_1       
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0004   gamma=    0.01 GHz
     lambda(    2)=  0.0004   gamma=    0.01 GHz
     lambda(    3)=  0.0021   gamma=    0.11 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0593   gamma=    1.88 GHz
     lambda(    2)=  0.0593   gamma=    1.88 GHz
     lambda(    3)=  0.0605   gamma=    4.56 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.1028   gamma=    3.68 GHz
     lambda(    2)=  0.1028   gamma=    3.68 GHz
     lambda(    3)=  0.0888   gamma=    7.56 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.1112   gamma=    4.21 GHz
     lambda(    2)=  0.1112   gamma=    4.21 GHz
     lambda(    3)=  0.1107   gamma=    9.98 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.1150   gamma=    4.51 GHz
     lambda(    2)=  0.1150   gamma=    4.51 GHz
     lambda(    3)=  0.1419   gamma=   13.25 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.1210   gamma=    4.88 GHz
     lambda(    2)=  0.1210   gamma=    4.88 GHz
     lambda(    3)=  0.1719   gamma=   16.50 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.1287   gamma=    5.32 GHz
     lambda(    2)=  0.1287   gamma=    5.32 GHz
     lambda(    3)=  0.1953   gamma=   19.22 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.1366   gamma=    5.77 GHz
     lambda(    2)=  0.1366   gamma=    5.77 GHz
     lambda(    3)=  0.2129   gamma=   21.40 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.1439   gamma=    6.19 GHz
     lambda(    2)=  0.1439   gamma=    6.19 GHz
     lambda(    3)=  0.2263   gamma=   23.15 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.1500   gamma=    6.54 GHz
     lambda(    2)=  0.1500   gamma=    6.54 GHz
     lambda(    3)=  0.2365   gamma=   24.55 GHz


     Number of q in the star =    6
     List of q in the star:
          1   0.000000000   0.500000000   0.000000000
          2   0.000000000  -0.500000000   0.000000000
          3   0.500000000   0.000000000   0.000000000
          4   0.000000000   0.000000000   0.500000000
          5   0.000000000   0.000000000  -0.500000000
          6  -0.500000000   0.000000000   0.000000000

     Calculation of q =    0.7500000  -0.2500000   0.7500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     15                  216      216      84
     Max          31      31     16                  218      218      87
     Sum         121     121     61                  869      869     339


     Title: 
     Electron-phonon coefficients for Al                                        


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   576  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.12Mb

     Estimated total allocated dynamical RAM >       0.47Mb

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_5/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.5

     total cpu time spent up to now is        6.7 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.7500000  -0.2500000   0.7500000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =   85.4897  (    218 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   576  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A''  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :    16.16s CPU        17.80s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :    18.0 secs   av.it.:   4.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.561E-04

      iter #   2 total cpu time :    18.3 secs   av.it.:   5.4
      thresh= 1.250E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.316E-04

      iter #   3 total cpu time :    18.5 secs   av.it.:   4.7
      thresh= 1.522E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.117E-07

      iter #   4 total cpu time :    18.7 secs   av.it.:   5.7
      thresh= 3.343E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.666E-09

      iter #   5 total cpu time :    18.9 secs   av.it.:   5.6
      thresh= 5.163E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.501E-10

      iter #   6 total cpu time :    19.1 secs   av.it.:   5.6
      thresh= 1.225E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.883E-12

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :    19.6 secs   av.it.:   4.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.262E-05

      iter #   2 total cpu time :    19.8 secs   av.it.:   5.7
      thresh= 5.711E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.558E-05

      iter #   3 total cpu time :    20.0 secs   av.it.:   5.0
      thresh= 5.965E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.597E-07

      iter #   4 total cpu time :    20.3 secs   av.it.:   5.4
      thresh= 6.780E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.987E-09

      iter #   5 total cpu time :    20.5 secs   av.it.:   5.7
      thresh= 4.457E-06 alpha_mix =  0.700 |ddv_scf|^2 =  8.539E-11

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :    21.0 secs   av.it.:   3.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.876E-07

      iter #   2 total cpu time :    21.2 secs   av.it.:   5.4
      thresh= 6.983E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.225E-08

      iter #   3 total cpu time :    21.4 secs   av.it.:   4.8
      thresh= 1.107E-05 alpha_mix =  0.700 |ddv_scf|^2 =  9.339E-10

      iter #   4 total cpu time :    21.6 secs   av.it.:   4.9
      thresh= 3.056E-06 alpha_mix =  0.700 |ddv_scf|^2 =  8.374E-14

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =   24
     List of q in the star:
          1   0.750000000  -0.250000000   0.750000000
          2   0.750000000  -0.250000000  -0.750000000
          3  -0.750000000  -0.250000000  -0.750000000
          4  -0.750000000  -0.250000000   0.750000000
          5  -0.750000000   0.250000000  -0.750000000
          6  -0.250000000   0.750000000  -0.750000000
          7  -0.750000000   0.750000000  -0.250000000
          8   0.750000000   0.250000000   0.750000000
          9  -0.750000000   0.250000000   0.750000000
         10   0.750000000   0.250000000  -0.750000000
         11  -0.750000000   0.750000000   0.250000000
         12  -0.250000000   0.750000000   0.750000000
         13   0.250000000   0.750000000  -0.750000000
         14  -0.250000000  -0.750000000  -0.750000000
         15   0.750000000   0.750000000  -0.250000000
         16   0.750000000  -0.750000000   0.250000000
         17  -0.750000000  -0.750000000  -0.250000000
         18   0.250000000  -0.750000000   0.750000000
         19  -0.750000000  -0.750000000   0.250000000
         20   0.250000000   0.750000000   0.750000000
         21  -0.250000000  -0.750000000   0.750000000
         22   0.750000000   0.750000000   0.250000000
         23   0.250000000  -0.750000000  -0.750000000
         24   0.750000000  -0.750000000  -0.250000000

     Diagonalizing the dynamical matrix

     q = (    0.750000000  -0.250000000   0.750000000 ) 

 **************************************************************************
     freq (    1) =       5.392336 [THz] =     179.868957 [cm-1]
     freq (    2) =       6.727093 [THz] =     224.391665 [cm-1]
     freq (    3) =       8.791383 [THz] =     293.248982 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =        179.9  [cm-1]   --> A''                
     freq (  2 -  2) =        224.4  [cm-1]   --> A'                 
     freq (  3 -  3) =        293.2  [cm-1]   --> A'                 
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0085   gamma=    0.32 GHz
     lambda(    2)=  0.0210   gamma=    1.22 GHz
     lambda(    3)=  0.0282   gamma=    2.79 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0619   gamma=    3.23 GHz
     lambda(    2)=  0.1351   gamma=   10.99 GHz
     lambda(    3)=  0.2006   gamma=   27.86 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.0789   gamma=    4.65 GHz
     lambda(    2)=  0.1337   gamma=   12.27 GHz
     lambda(    3)=  0.2248   gamma=   35.23 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.0854   gamma=    5.34 GHz
     lambda(    2)=  0.1171   gamma=   11.39 GHz
     lambda(    3)=  0.2243   gamma=   37.24 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.0863   gamma=    5.58 GHz
     lambda(    2)=  0.1046   gamma=   10.54 GHz
     lambda(    3)=  0.2160   gamma=   37.14 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.0867   gamma=    5.77 GHz
     lambda(    2)=  0.0977   gamma=   10.12 GHz
     lambda(    3)=  0.2082   gamma=   36.81 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.0876   gamma=    5.97 GHz
     lambda(    2)=  0.0948   gamma=   10.05 GHz
     lambda(    3)=  0.2033   gamma=   36.84 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.0889   gamma=    6.19 GHz
     lambda(    2)=  0.0942   gamma=   10.21 GHz
     lambda(    3)=  0.2011   gamma=   37.23 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.0905   gamma=    6.41 GHz
     lambda(    2)=  0.0950   gamma=   10.48 GHz
     lambda(    3)=  0.2009   gamma=   37.85 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.0919   gamma=    6.61 GHz
     lambda(    2)=  0.0963   gamma=   10.78 GHz
     lambda(    3)=  0.2017   gamma=   38.54 GHz


     Number of q in the star =   24
     List of q in the star:
          1   0.750000000  -0.250000000   0.750000000
          2   0.750000000  -0.250000000  -0.750000000
          3  -0.750000000  -0.250000000  -0.750000000
          4  -0.750000000  -0.250000000   0.750000000
          5  -0.750000000   0.250000000  -0.750000000
          6  -0.250000000   0.750000000  -0.750000000
          7  -0.750000000   0.750000000  -0.250000000
          8   0.750000000   0.250000000   0.750000000
          9  -0.750000000   0.250000000   0.750000000
         10   0.750000000   0.250000000  -0.750000000
         11  -0.750000000   0.750000000   0.250000000
         12  -0.250000000   0.750000000   0.750000000
         13   0.250000000   0.750000000  -0.750000000
         14  -0.250000000  -0.750000000  -0.750000000
         15   0.750000000   0.750000000  -0.250000000
         16   0.750000000  -0.750000000   0.250000000
         17  -0.750000000  -0.750000000  -0.250000000
         18   0.250000000  -0.750000000   0.750000000
         19  -0.750000000  -0.750000000   0.250000000
         20   0.250000000   0.750000000   0.750000000
         21  -0.250000000  -0.750000000   0.750000000
         22   0.750000000   0.750000000   0.250000000
         23   0.250000000  -0.750000000  -0.750000000
         24   0.750000000  -0.750000000  -0.250000000

     Calculation of q =    0.5000000   0.0000000   0.5000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     15                  217      217      76
     Max          31      31     16                  218      218      77
     Sum         121     121     61                  869      869     307


     Title: 
     Electron-phonon coefficients for Al                                        


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   328  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.12Mb

     Estimated total allocated dynamical RAM >       0.47Mb

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_6/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.5

     total cpu time spent up to now is        8.9 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.5000000   0.0000000   0.5000000 )

      4 Sym.Ops. (no q -> -q+G )


     G cutoff =   85.4897  (    217 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   328  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  To be done

     Representation     2      1 modes -B_1  D_3  S_3  To be done

     Representation     3      1 modes -B_2  D_4  S_4  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :    22.96s CPU        25.48s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :    25.6 secs   av.it.:   4.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.825E-04

      iter #   2 total cpu time :    25.7 secs   av.it.:   4.9
      thresh= 2.414E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.048E-03

      iter #   3 total cpu time :    25.8 secs   av.it.:   4.1
      thresh= 4.525E-03 alpha_mix =  0.700 |ddv_scf|^2 =  4.208E-08

      iter #   4 total cpu time :    26.0 secs   av.it.:   5.9
      thresh= 2.051E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.527E-09

      iter #   5 total cpu time :    26.1 secs   av.it.:   5.4
      thresh= 3.908E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.511E-11

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :    26.4 secs   av.it.:   3.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.784E-07

      iter #   2 total cpu time :    26.5 secs   av.it.:   5.2
      thresh= 4.224E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.286E-08

      iter #   3 total cpu time :    26.6 secs   av.it.:   5.0
      thresh= 1.134E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.929E-10

      iter #   4 total cpu time :    26.7 secs   av.it.:   5.2
      thresh= 1.389E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.413E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :    27.0 secs   av.it.:   3.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.011E-06

      iter #   2 total cpu time :    27.1 secs   av.it.:   5.4
      thresh= 2.239E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.310E-07

      iter #   3 total cpu time :    27.3 secs   av.it.:   5.3
      thresh= 5.753E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.653E-09

      iter #   4 total cpu time :    27.4 secs   av.it.:   5.2
      thresh= 5.151E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.745E-12

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =   12
     List of q in the star:
          1   0.500000000   0.000000000   0.500000000
          2  -0.500000000   0.000000000   0.500000000
          3  -0.500000000   0.000000000  -0.500000000
          4   0.500000000   0.000000000  -0.500000000
          5   0.000000000   0.500000000  -0.500000000
          6  -0.500000000   0.500000000   0.000000000
          7   0.000000000   0.500000000   0.500000000
          8   0.000000000  -0.500000000  -0.500000000
          9   0.500000000   0.500000000   0.000000000
         10   0.500000000  -0.500000000   0.000000000
         11  -0.500000000  -0.500000000   0.000000000
         12   0.000000000  -0.500000000   0.500000000

     Diagonalizing the dynamical matrix

     q = (    0.500000000   0.000000000   0.500000000 ) 

 **************************************************************************
     freq (    1) =       4.864075 [THz] =     162.248094 [cm-1]
     freq (    2) =       6.528731 [THz] =     217.775011 [cm-1]
     freq (    3) =       8.467305 [THz] =     282.438904 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =        162.2  [cm-1]   --> B_1  D_3  S_3      
     freq (  2 -  2) =        217.8  [cm-1]   --> B_2  D_4  S_4      
     freq (  3 -  3) =        282.4  [cm-1]   --> A_1  D_1  S_1      
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0231   gamma=    0.70 GHz
     lambda(    2)=  0.0561   gamma=    3.06 GHz
     lambda(    3)=  1.3275   gamma=  121.72 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0651   gamma=    2.77 GHz
     lambda(    2)=  0.0805   gamma=    6.17 GHz
     lambda(    3)=  0.8798   gamma=  113.35 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.0533   gamma=    2.56 GHz
     lambda(    2)=  0.1119   gamma=    9.67 GHz
     lambda(    3)=  0.5477   gamma=   79.62 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.0426   gamma=    2.16 GHz
     lambda(    2)=  0.1260   gamma=   11.53 GHz
     lambda(    3)=  0.3883   gamma=   59.81 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.0363   gamma=    1.91 GHz
     lambda(    2)=  0.1256   gamma=   11.91 GHz
     lambda(    3)=  0.3073   gamma=   49.02 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.0334   gamma=    1.81 GHz
     lambda(    2)=  0.1249   gamma=   12.18 GHz
     lambda(    3)=  0.2655   gamma=   43.55 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.0324   gamma=    1.80 GHz
     lambda(    2)=  0.1264   gamma=   12.63 GHz
     lambda(    3)=  0.2437   gamma=   40.96 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.0323   gamma=    1.83 GHz
     lambda(    2)=  0.1290   gamma=   13.16 GHz
     lambda(    3)=  0.2317   gamma=   39.79 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.0325   gamma=    1.87 GHz
     lambda(    2)=  0.1316   gamma=   13.68 GHz
     lambda(    3)=  0.2246   gamma=   39.25 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.0329   gamma=    1.92 GHz
     lambda(    2)=  0.1338   gamma=   14.10 GHz
     lambda(    3)=  0.2196   gamma=   38.94 GHz


     Number of q in the star =   12
     List of q in the star:
          1   0.500000000   0.000000000   0.500000000
          2  -0.500000000   0.000000000   0.500000000
          3  -0.500000000   0.000000000  -0.500000000
          4   0.500000000   0.000000000  -0.500000000
          5   0.000000000   0.500000000  -0.500000000
          6  -0.500000000   0.500000000   0.000000000
          7   0.000000000   0.500000000   0.500000000
          8   0.000000000  -0.500000000  -0.500000000
          9   0.500000000   0.500000000   0.000000000
         10   0.500000000  -0.500000000   0.000000000
         11  -0.500000000  -0.500000000   0.000000000
         12   0.000000000  -0.500000000   0.500000000

     Calculation of q =    0.0000000  -1.0000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     15                  216      216      82
     Max          31      31     16                  218      218      83
     Sum         121     121     61                  869      869     331


     Title: 
     Electron-phonon coefficients for Al                                        


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   118  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.12Mb

     Estimated total allocated dynamical RAM >       0.47Mb

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_7/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.3

     total cpu time spent up to now is        9.9 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000  -1.0000000   0.0000000 )

     17 Sym.Ops. (with q -> -q+G )


     G cutoff =   85.4897  (    218 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   118  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, D_4h(4/mmm) point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_2u X_4' M_4'  To be done

     Representation     2      2 modes -E_u  X_5' M_5'  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :    26.92s CPU        29.93s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :    30.0 secs   av.it.:   3.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.091E-05

      iter #   2 total cpu time :    30.0 secs   av.it.:   5.1
      thresh= 7.804E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.345E-05

      iter #   3 total cpu time :    30.1 secs   av.it.:   4.9
      thresh= 4.843E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.559E-09

      iter #   4 total cpu time :    30.1 secs   av.it.:   5.1
      thresh= 6.752E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.523E-11

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   2  3

     Self-consistent Calculation

      iter #   1 total cpu time :    30.3 secs   av.it.:   3.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.779E-07

      iter #   2 total cpu time :    30.4 secs   av.it.:   5.9
      thresh= 5.271E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.926E-09

      iter #   3 total cpu time :    30.5 secs   av.it.:   5.5
      thresh= 6.266E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.700E-10

      iter #   4 total cpu time :    30.6 secs   av.it.:   5.4
      thresh= 1.923E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.519E-14

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    3
     List of q in the star:
          1   0.000000000  -1.000000000   0.000000000
          2  -1.000000000   0.000000000   0.000000000
          3   0.000000000   0.000000000  -1.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000  -1.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       6.062697 [THz] =     202.229809 [cm-1]
     freq (    2) =       6.062697 [THz] =     202.229809 [cm-1]
     freq (    3) =       9.881070 [THz] =     329.597010 [cm-1]
 **************************************************************************

     Mode symmetry, D_4h(4/mmm) point group:

     freq (  1 -  2) =        202.2  [cm-1]   --> E_u  X_5' M_5'     
     freq (  3 -  3) =        329.6  [cm-1]   --> A_2u X_4' M_4'     
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0244   gamma=    1.15 GHz
     lambda(    2)=  0.0244   gamma=    1.15 GHz
     lambda(    3)=  0.0002   gamma=    0.02 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.1849   gamma=   12.21 GHz
     lambda(    2)=  0.1833   gamma=   12.11 GHz
     lambda(    3)=  0.0909   gamma=   15.95 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.1777   gamma=   13.24 GHz
     lambda(    2)=  0.1681   gamma=   12.53 GHz
     lambda(    3)=  0.1880   gamma=   37.21 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.1597   gamma=   12.61 GHz
     lambda(    2)=  0.1445   gamma=   11.41 GHz
     lambda(    3)=  0.2032   gamma=   42.62 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.1482   gamma=   12.12 GHz
     lambda(    2)=  0.1312   gamma=   10.73 GHz
     lambda(    3)=  0.1871   gamma=   40.65 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.1396   gamma=   11.74 GHz
     lambda(    2)=  0.1229   gamma=   10.34 GHz
     lambda(    3)=  0.1681   gamma=   37.55 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.1327   gamma=   11.44 GHz
     lambda(    2)=  0.1169   gamma=   10.07 GHz
     lambda(    3)=  0.1551   gamma=   35.49 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.1273   gamma=   11.21 GHz
     lambda(    2)=  0.1124   gamma=    9.89 GHz
     lambda(    3)=  0.1484   gamma=   34.70 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.1234   gamma=   11.06 GHz
     lambda(    2)=  0.1092   gamma=    9.78 GHz
     lambda(    3)=  0.1458   gamma=   34.72 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.1207   gamma=   10.97 GHz
     lambda(    2)=  0.1071   gamma=    9.73 GHz
     lambda(    3)=  0.1455   gamma=   35.12 GHz


     Number of q in the star =    3
     List of q in the star:
          1   0.000000000  -1.000000000   0.000000000
          2  -1.000000000   0.000000000   0.000000000
          3   0.000000000   0.000000000  -1.000000000

     Calculation of q =   -0.5000000  -1.0000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          30      30     15                  216      216      82
     Max          31      31     16                  218      218      83
     Sum         121     121     61                  869      869     331


     Title: 
     Electron-phonon coefficients for Al                                        


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   174  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.12Mb

     Estimated total allocated dynamical RAM >       0.47Mb

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_8/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.3

     total cpu time spent up to now is       11.0 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

     Electron-phonon coefficients for Al                                        

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.5000000  -1.0000000   0.0000000 )

      8 Sym.Ops. (no q -> -q+G )


     G cutoff =   85.4897  (    218 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   174  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     ./Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, D_2d (-42m) point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -B_2  X_3  W_2  To be done

     Representation     2      2 modes -E    X_5  W_3  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :    30.07s CPU        33.35s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :    33.4 secs   av.it.:   3.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.363E-06

      iter #   2 total cpu time :    33.5 secs   av.it.:   5.5
      thresh= 2.892E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.288E-06

      iter #   3 total cpu time :    33.6 secs   av.it.:   5.4
      thresh= 1.135E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.421E-09

      iter #   4 total cpu time :    33.6 secs   av.it.:   5.3
      thresh= 5.849E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.815E-12

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   2  3

     Self-consistent Calculation

      iter #   1 total cpu time :    33.9 secs   av.it.:   4.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.840E-06

      iter #   2 total cpu time :    34.0 secs   av.it.:   5.9
      thresh= 2.200E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.031E-06

      iter #   3 total cpu time :    34.1 secs   av.it.:   5.8
      thresh= 1.015E-04 alpha_mix =  0.700 |ddv_scf|^2 =  9.253E-10

      iter #   4 total cpu time :    34.3 secs   av.it.:   5.8
      thresh= 3.042E-06 alpha_mix =  0.700 |ddv_scf|^2 =  8.860E-13

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    6
     List of q in the star:
          1  -0.500000000  -1.000000000   0.000000000
          2   0.000000000   1.000000000   0.500000000
          3   0.000000000  -1.000000000  -0.500000000
          4   0.500000000   1.000000000   0.000000000
          5  -1.000000000  -0.500000000   0.000000000
          6   0.000000000  -0.500000000  -1.000000000

     Diagonalizing the dynamical matrix

     q = (   -0.500000000  -1.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       6.453881 [THz] =     215.278296 [cm-1]
     freq (    2) =       7.605739 [THz] =     253.700152 [cm-1]
     freq (    3) =       7.605739 [THz] =     253.700152 [cm-1]
 **************************************************************************

     Mode symmetry, D_2d (-42m) point group:

     freq (  1 -  1) =        215.3  [cm-1]   --> B_2  X_3  W_2      
     freq (  2 -  3) =        253.7  [cm-1]   --> E    X_5  W_3      
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0002   gamma=    0.01 GHz
     lambda(    2)=  0.0004   gamma=    0.03 GHz
     lambda(    3)=  0.0004   gamma=    0.03 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0330   gamma=    2.47 GHz
     lambda(    2)=  0.0635   gamma=    6.60 GHz
     lambda(    3)=  0.0631   gamma=    6.55 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.0685   gamma=    5.79 GHz
     lambda(    2)=  0.1081   gamma=   12.68 GHz
     lambda(    3)=  0.1074   gamma=   12.60 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.0740   gamma=    6.62 GHz
     lambda(    2)=  0.1144   gamma=   14.21 GHz
     lambda(    3)=  0.1141   gamma=   14.17 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.0718   gamma=    6.65 GHz
     lambda(    2)=  0.1157   gamma=   14.89 GHz
     lambda(    3)=  0.1158   gamma=   14.91 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.0695   gamma=    6.63 GHz
     lambda(    2)=  0.1197   gamma=   15.85 GHz
     lambda(    3)=  0.1201   gamma=   15.90 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.0687   gamma=    6.71 GHz
     lambda(    2)=  0.1268   gamma=   17.20 GHz
     lambda(    3)=  0.1273   gamma=   17.27 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299955 eV
     lambda(    1)=  0.0694   gamma=    6.92 GHz
     lambda(    2)=  0.1358   gamma=   18.82 GHz
     lambda(    3)=  0.1364   gamma=   18.89 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.0712   gamma=    7.23 GHz
     lambda(    2)=  0.1453   gamma=   20.49 GHz
     lambda(    3)=  0.1458   gamma=   20.56 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.0735   gamma=    7.57 GHz
     lambda(    2)=  0.1539   gamma=   22.02 GHz
     lambda(    3)=  0.1544   gamma=   22.09 GHz


     Number of q in the star =    6
     List of q in the star:
          1  -0.500000000  -1.000000000   0.000000000
          2   0.000000000   1.000000000   0.500000000
          3   0.000000000  -1.000000000  -0.500000000
          4   0.500000000   1.000000000   0.000000000
          5  -1.000000000  -0.500000000   0.000000000
          6   0.000000000  -0.500000000  -1.000000000

     init_run     :      0.03s CPU      0.05s WALL (       7 calls)
     electrons    :      8.11s CPU      9.06s WALL (       7 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       7 calls)
     potinit      :      0.00s CPU      0.01s WALL (       7 calls)

     Called by electrons:
     c_bands      :      8.10s CPU      9.04s WALL (       7 calls)
     v_of_rho     :      0.00s CPU      0.00s WALL (       8 calls)

     Called by c_bands:
     init_us_2    :      0.17s CPU      0.19s WALL (   18420 calls)
     cegterg      :      7.62s CPU      8.44s WALL (    1847 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :      7.18s CPU      8.43s WALL (   96516 calls)
     g_psi        :      0.06s CPU      0.04s WALL (   23871 calls)
     cdiaghg      :      3.09s CPU      3.42s WALL (   25637 calls)

     Called by h_psi:
     h_psi:pot    :      7.02s CPU      8.25s WALL (   96516 calls)
     h_psi:calbec :      0.55s CPU      0.69s WALL (   96516 calls)
     vloc_psi     :      6.00s CPU      7.03s WALL (   96516 calls)
     add_vuspsi   :      0.18s CPU      0.23s WALL (   96516 calls)

     General routines
     calbec       :      0.80s CPU      1.00s WALL (  185641 calls)
     fft          :      0.02s CPU      0.02s WALL (     344 calls)
     ffts         :      0.04s CPU      0.04s WALL (    2839 calls)
     fftw         :      6.50s CPU      7.80s WALL (  605934 calls)
     davcio       :      0.18s CPU      0.26s WALL (   68861 calls)

     Parallel routines
     fft_scatter  :      2.67s CPU      3.17s WALL (  609117 calls)

     PHONON       :    32.49s CPU        35.98s WALL

     INITIALIZATION: 
     phq_setup    :      0.02s CPU      0.03s WALL (       8 calls)
     phq_init     :      0.10s CPU      0.14s WALL (       8 calls)

     phq_init     :      0.10s CPU      0.14s WALL (       8 calls)
     init_vloc    :      0.00s CPU      0.00s WALL (       8 calls)
     init_us_1    :      0.02s CPU      0.01s WALL (       8 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.04s CPU      0.08s WALL (       8 calls)
     phqscf       :      9.42s CPU     11.47s WALL (       8 calls)
     dynmatrix    :      0.02s CPU      0.02s WALL (       8 calls)

     phqscf       :      9.42s CPU     11.47s WALL (       8 calls)
     solve_linter :      9.26s CPU     11.28s WALL (      17 calls)
     drhodv       :      0.12s CPU      0.15s WALL (      17 calls)

     dynmat0      :      0.04s CPU      0.08s WALL (       8 calls)
     dynmat_us    :      0.04s CPU      0.07s WALL (       8 calls)
     d2ionq       :      0.00s CPU      0.00s WALL (       8 calls)

     dynmat_us    :      0.04s CPU      0.07s WALL (       8 calls)

     phqscf       :      9.42s CPU     11.47s WALL (       8 calls)
     solve_linter :      9.26s CPU     11.28s WALL (      17 calls)

     solve_linter :      9.26s CPU     11.28s WALL (      17 calls)
     dvqpsi_us    :      0.52s CPU      0.64s WALL (    2736 calls)
     ortho        :      0.15s CPU      0.16s WALL (   12020 calls)
     cgsolve      :      5.24s CPU      6.47s WALL (   12020 calls)
     incdrhoscf   :      0.49s CPU      0.68s WALL (   12020 calls)
     vpsifft      :      0.34s CPU      0.56s WALL (    9284 calls)
     dv_of_drho   :      0.02s CPU      0.02s WALL (      98 calls)
     mix_pot      :      0.01s CPU      0.01s WALL (      74 calls)
     ef_shift     :      0.00s CPU      0.00s WALL (       4 calls)
     localdos     :      0.00s CPU      0.00s WALL (       1 calls)
     psymdvscf    :      0.18s CPU      0.20s WALL (      74 calls)

     dvqpsi_us    :      0.52s CPU      0.64s WALL (    2736 calls)
     dvqpsi_us_on :      0.07s CPU      0.07s WALL (    2736 calls)

     cgsolve      :      5.24s CPU      6.47s WALL (   12020 calls)
     ch_psi       :      4.61s CPU      5.72s WALL (   69032 calls)

     ch_psi       :      4.61s CPU      5.72s WALL (   69032 calls)
     h_psi        :      7.18s CPU      8.43s WALL (   96516 calls)
     last         :      0.75s CPU      0.89s WALL (   69032 calls)

     h_psi        :      7.18s CPU      8.43s WALL (   96516 calls)
     add_vuspsi   :      0.18s CPU      0.23s WALL (   96516 calls)

     incdrhoscf   :      0.49s CPU      0.68s WALL (   12020 calls)


      General routines
     calbec       :      0.80s CPU      1.00s WALL (  185641 calls)
     fft          :      0.02s CPU      0.02s WALL (     344 calls)
     ffts         :      0.04s CPU      0.04s WALL (    2839 calls)
     fftw         :      6.50s CPU      7.80s WALL (  605934 calls)
     davcio       :      0.18s CPU      0.26s WALL (   68861 calls)
     write_rec    :      0.10s CPU      0.14s WALL (      91 calls)


     PHONON       :    32.49s CPU        35.98s WALL


   This run was terminated on:  10:43:43   7Dec2016            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_6 = """
     Program PHONON v.6.2 (svn rev. 13952M) starts on 25Oct2017 at 20:58:55 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
         "P. Giannozzi et al., J. Phys.:Condens. Matter 29 465901 (2017);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors

     MPI processes distributed on     1 nodes
     R & G space division:  proc/nbgrp/npool/nimage =       4

     Reading data from directory:
     /home/paulatto/espresso/tempdir/bn.save/

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want

               file B.pbe-n-kjpaw_psl.0.1.UPF: wavefunction(s)  2P renormalized
               file N.pbe-n-kjpaw_psl.0.1.UPF: wavefunction(s)  2P renormalized
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     27                15432     5472     963
     Max         169      88     28                15455     5478     974
     Sum         673     349    109                61757    21901    3869
 
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     negative rho (up, down):  1.316E-05 0.000E+00


     Dynamical matrices for ( 8, 8, 1)  uniform grid of q-points
     (  10 q-points):
       N         xq(1)         xq(2)         xq(3) 
       1   0.000000000   0.000000000   0.000000000
       2   0.000000000   0.144337567   0.000000000
       3   0.000000000   0.288675135   0.000000000
       4   0.000000000   0.433012702   0.000000000
       5   0.000000000  -0.577350269   0.000000000
       6   0.125000000   0.216506351   0.000000000
       7   0.125000000   0.360843918   0.000000000
       8   0.125000000   0.505181486   0.000000000
       9   0.250000000   0.433012702   0.000000000
      10   0.250000000   0.577350269   0.000000000

     Calculation of q =    0.0000000   0.0000000   0.0000000

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   0.0000000 )
 
     13 Sym.Ops. (with q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=    19

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, D_3h (-62m) point group:


     Electric field:
     Dielectric constant
     Born effective charges in two ways 


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      2 modes -E'  To be done

     Representation     2      2 modes -E'  To be done

     Representation     3      1 modes -A''2  To be done

     Representation     4      1 modes -A''2  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       :     3.18s CPU         3.30s WALL


     Electric Fields Calculation

      iter #   1 total cpu time :     8.9 secs   av.it.:   6.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.021E-05

      iter #   2 total cpu time :    11.3 secs   av.it.:   9.0
      thresh= 4.496E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.389E-05

      iter #   3 total cpu time :    13.2 secs   av.it.:   6.3
      thresh= 4.887E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.428E-08

      iter #   4 total cpu time :    15.6 secs   av.it.:   9.7
      thresh= 1.195E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.527E-10

      iter #   5 total cpu time :    18.2 secs   av.it.:   9.9
      thresh= 1.236E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.167E-12

      iter #   6 total cpu time :    20.9 secs   av.it.:  10.1
      thresh= 1.780E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.960E-14

      iter #   7 total cpu time :    23.7 secs   av.it.:  10.8
      thresh= 1.400E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.920E-17

     End of electric fields calculation

          Dielectric constant in cartesian axis 

          (       2.165990189       0.000000000       0.000000000 )
          (       0.000000000       2.165990189       0.000000000 )
          (       0.000000000       0.000000000       1.192206706 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   B  
      Ex  (        2.71348       -0.00000        0.00000 )
      Ey  (       -0.00000        2.71348        0.00000 )
      Ez  (       -0.00000       -0.00000        0.24303 )
           atom      2   N  
      Ex  (       -2.72028       -0.00000        0.00000 )
      Ey  (        0.00000       -2.72028        0.00000 )
      Ez  (        0.00000       -0.00000       -0.24985 )


     Representation #  1 modes #   1  2

     Self-consistent Calculation

      iter #   1 total cpu time :    26.9 secs   av.it.:   7.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.077E-08

      iter #   2 total cpu time :    29.0 secs   av.it.:  12.2
      thresh= 2.842E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.549E-10

      iter #   3 total cpu time :    31.0 secs   av.it.:  11.5
      thresh= 1.596E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.363E-11

      iter #   4 total cpu time :    33.0 secs   av.it.:  11.5
      thresh= 4.862E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.161E-13

      iter #   5 total cpu time :    35.0 secs   av.it.:  11.1
      thresh= 4.649E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.987E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   3  4

     Self-consistent Calculation

      iter #   1 total cpu time :    36.8 secs   av.it.:   6.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.841E-08

      iter #   2 total cpu time :    38.7 secs   av.it.:  10.7
      thresh= 2.616E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.916E-09

      iter #   3 total cpu time :    40.6 secs   av.it.:  10.2
      thresh= 7.691E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.868E-12

      iter #   4 total cpu time :    42.6 secs   av.it.:  11.0
      thresh= 2.621E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.186E-13

      iter #   5 total cpu time :    44.5 secs   av.it.:  10.8
      thresh= 5.644E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.291E-15

      iter #   6 total cpu time :    46.5 secs   av.it.:  11.9
      thresh= 3.593E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.887E-18

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :    47.6 secs   av.it.:   6.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.516E-05

      iter #   2 total cpu time :    48.6 secs   av.it.:   9.3
      thresh= 5.929E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.433E-05

      iter #   3 total cpu time :    49.5 secs   av.it.:   8.7
      thresh= 6.658E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.465E-08

      iter #   4 total cpu time :    50.4 secs   av.it.:   8.4
      thresh= 1.861E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.318E-10

      iter #   5 total cpu time :    51.4 secs   av.it.:   8.7
      thresh= 1.148E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.634E-12

      iter #   6 total cpu time :    52.3 secs   av.it.:   9.3
      thresh= 1.906E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.356E-13

      iter #   7 total cpu time :    53.3 secs   av.it.:   9.1
      thresh= 3.682E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.795E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :    54.3 secs   av.it.:   7.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.631E-04

      iter #   2 total cpu time :    55.3 secs   av.it.:   9.4
      thresh= 1.277E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.325E-04

      iter #   3 total cpu time :    56.2 secs   av.it.:   8.7
      thresh= 1.525E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.359E-07

      iter #   4 total cpu time :    57.1 secs   av.it.:   8.9
      thresh= 4.857E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.263E-10

      iter #   5 total cpu time :    58.1 secs   av.it.:   9.3
      thresh= 1.806E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.069E-12

      iter #   6 total cpu time :    59.1 secs   av.it.:   9.3
      thresh= 1.439E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.259E-13

      iter #   7 total cpu time :    60.1 secs   av.it.:   9.3
      thresh= 3.548E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.545E-17

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

          Dielectric constant in cartesian axis 

          (       2.165990189       0.000000000       0.000000000 )
          (       0.000000000       2.165990189       0.000000000 )
          (       0.000000000       0.000000000       1.192206706 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   B  
      Ex  (        2.71348       -0.00000        0.00000 )
      Ey  (       -0.00000        2.71348        0.00000 )
      Ez  (       -0.00000       -0.00000        0.24303 )
           atom      2   N  
      Ex  (       -2.72028       -0.00000        0.00000 )
      Ey  (        0.00000       -2.72028        0.00000 )
      Ez  (        0.00000       -0.00000       -0.24985 )

          Effective charges (d P / du) in cartesian axis 

           atom      1   B  
      Px  (        2.71347       -0.00000        0.00000 )
      Py  (        0.00000        2.71347        0.00000 )
      Pz  (       -0.00000       -0.00000        0.24303 )
           atom      2   N  
      Px  (       -2.72028        0.00000        0.00000 )
      Py  (        0.00000       -2.72028        0.00000 )
      Pz  (        0.00000        0.00000       -0.24984 )

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =      -1.430840 [THz] =     -47.727679 [cm-1]
     freq (    2) =      -1.430840 [THz] =     -47.727679 [cm-1]
     freq (    3) =      -0.707884 [THz] =     -23.612465 [cm-1]
     freq (    4) =      24.081625 [THz] =     803.276533 [cm-1]
     freq (    5) =      40.332074 [THz] =    1345.333180 [cm-1]
     freq (    6) =      40.332074 [THz] =    1345.333180 [cm-1]
 **************************************************************************

     Mode symmetry, D_3h (-62m) point group:

     freq (  1 -  2) =        -47.7  [cm-1]   --> E'              I+R
     freq (  3 -  3) =        -23.6  [cm-1]   --> A''2            I  
     freq (  4 -  4) =        803.3  [cm-1]   --> A''2            I  
     freq (  5 -  6) =       1345.3  [cm-1]   --> E'              I+R

     Calculation of q =    0.0000000   0.1443376   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     27                15432     5472    1013
     Max         169      88     28                15455     5478    1030
     Sum         673     349    109                61757    21901    4075
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   156

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.42 MB

     Estimated total dynamical RAM >      65.70 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.8

     total cpu time spent up to now is        8.2 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.1443376   0.0000000 )
 
      4 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   156

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  To be done

     Representation     2      1 modes -A_1  D_1  S_1  To be done

     Representation     3      1 modes -B_1  D_3  S_3  To be done

     Representation     4      1 modes -B_1  D_3  S_3  To be done

     Representation     5      1 modes -B_2  D_4  S_4  To be done

     Representation     6      1 modes -B_2  D_4  S_4  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       :  1m 9.06s CPU     1m11.52s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :    73.8 secs   av.it.:   6.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.284E-06

      iter #   2 total cpu time :    76.7 secs   av.it.:  10.3
      thresh= 2.299E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.440E-06

      iter #   3 total cpu time :    79.4 secs   av.it.:  10.0
      thresh= 2.107E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.841E-08

      iter #   4 total cpu time :    82.3 secs   av.it.:  10.5
      thresh= 1.960E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.857E-09

      iter #   5 total cpu time :    85.2 secs   av.it.:  10.6
      thresh= 7.653E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.034E-10

      iter #   6 total cpu time :    88.0 secs   av.it.:  10.1
      thresh= 1.426E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.109E-13

      iter #   7 total cpu time :    90.8 secs   av.it.:  10.0
      thresh= 7.148E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.941E-14

      iter #   8 total cpu time :    93.7 secs   av.it.:  10.4
      thresh= 2.635E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.306E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :    96.2 secs   av.it.:   7.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.129E-05

      iter #   2 total cpu time :    99.0 secs   av.it.:  10.3
      thresh= 7.162E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.955E-05

      iter #   3 total cpu time :   101.7 secs   av.it.:   9.5
      thresh= 7.039E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.627E-08

      iter #   4 total cpu time :   104.8 secs   av.it.:  11.3
      thresh= 1.905E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.919E-08

      iter #   5 total cpu time :   107.8 secs   av.it.:  10.5
      thresh= 1.708E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.594E-11

      iter #   6 total cpu time :   110.8 secs   av.it.:  10.9
      thresh= 5.093E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.137E-12

      iter #   7 total cpu time :   113.7 secs   av.it.:  10.6
      thresh= 1.066E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.116E-14

      iter #   8 total cpu time :   116.7 secs   av.it.:  11.1
      thresh= 1.056E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.876E-15

      iter #   9 total cpu time :   119.6 secs   av.it.:  10.3
      thresh= 8.292E-09 alpha_mix =  0.700 |ddv_scf|^2 =  4.032E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :   122.0 secs   av.it.:   7.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.511E-07

      iter #   2 total cpu time :   125.1 secs   av.it.:  11.7
      thresh= 5.925E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.153E-09

      iter #   3 total cpu time :   128.1 secs   av.it.:  10.9
      thresh= 3.395E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.959E-11

      iter #   4 total cpu time :   131.0 secs   av.it.:  10.9
      thresh= 9.979E-07 alpha_mix =  0.700 |ddv_scf|^2 =  8.957E-13

      iter #   5 total cpu time :   134.0 secs   av.it.:  10.7
      thresh= 9.464E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.971E-15

      iter #   6 total cpu time :   137.1 secs   av.it.:  11.6
      thresh= 5.451E-09 alpha_mix =  0.700 |ddv_scf|^2 =  9.316E-18

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :   139.3 secs   av.it.:   6.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.510E-07

      iter #   2 total cpu time :   142.2 secs   av.it.:  10.1
      thresh= 5.924E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.038E-08

      iter #   3 total cpu time :   144.9 secs   av.it.:   9.4
      thresh= 1.743E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.180E-11

      iter #   4 total cpu time :   147.9 secs   av.it.:  10.3
      thresh= 5.639E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.415E-12

      iter #   5 total cpu time :   150.9 secs   av.it.:  10.4
      thresh= 1.190E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.898E-15

      iter #   6 total cpu time :   154.1 secs   av.it.:  11.5
      thresh= 6.244E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.124E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :   156.5 secs   av.it.:   6.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.024E-06

      iter #   2 total cpu time :   159.5 secs   av.it.:  10.1
      thresh= 2.006E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.236E-06

      iter #   3 total cpu time :   162.4 secs   av.it.:  10.0
      thresh= 1.112E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.473E-08

      iter #   4 total cpu time :   165.3 secs   av.it.:   9.7
      thresh= 1.864E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.327E-11

      iter #   5 total cpu time :   168.2 secs   av.it.:  10.0
      thresh= 5.768E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.949E-13

      iter #   6 total cpu time :   171.2 secs   av.it.:  10.1
      thresh= 6.284E-08 alpha_mix =  0.700 |ddv_scf|^2 =  5.723E-15

      iter #   7 total cpu time :   174.1 secs   av.it.:  10.0
      thresh= 7.565E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.073E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :   176.6 secs   av.it.:   7.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.839E-05

      iter #   2 total cpu time :   179.6 secs   av.it.:  10.4
      thresh= 4.289E-04 alpha_mix =  0.700 |ddv_scf|^2 =  8.285E-06

      iter #   3 total cpu time :   182.6 secs   av.it.:  10.0
      thresh= 2.878E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.295E-08

      iter #   4 total cpu time :   185.5 secs   av.it.:   9.9
      thresh= 1.138E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.735E-11

      iter #   5 total cpu time :   188.4 secs   av.it.:  10.1
      thresh= 6.111E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.629E-13

      iter #   6 total cpu time :   191.2 secs   av.it.:   9.1
      thresh= 5.127E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.365E-14

      iter #   7 total cpu time :   194.1 secs   av.it.:  10.2
      thresh= 1.168E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.570E-17

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    3
     List of q in the star:
          1   0.000000000   0.144337567   0.000000000
          2  -0.125000000  -0.072168784   0.000000000
          3   0.125000000  -0.072168784   0.000000000
     In addition there is the -q list: 
          1   0.000000000  -0.144337567   0.000000000
          2   0.125000000   0.072168784   0.000000000
          3  -0.125000000   0.072168784   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.144337567   0.000000000 ) 

 **************************************************************************
     freq (    1) =       0.466982 [THz] =      15.576834 [cm-1]
     freq (    2) =       6.603116 [THz] =     220.256244 [cm-1]
     freq (    3) =      11.003220 [THz] =     367.027899 [cm-1]
     freq (    4) =      23.469946 [THz] =     782.873121 [cm-1]
     freq (    5) =      39.820184 [THz] =    1328.258361 [cm-1]
     freq (    6) =      45.411410 [THz] =    1514.761578 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =         15.6  [cm-1]   --> B_2  D_4  S_4      
     freq (  2 -  2) =        220.3  [cm-1]   --> B_1  D_3  S_3      
     freq (  3 -  3) =        367.0  [cm-1]   --> A_1  D_1  S_1      
     freq (  4 -  4) =        782.9  [cm-1]   --> B_2  D_4  S_4      
     freq (  5 -  5) =       1328.3  [cm-1]   --> B_1  D_3  S_3      
     freq (  6 -  6) =       1514.8  [cm-1]   --> A_1  D_1  S_1      

     Calculation of q =    0.0000000   0.2886751   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     29                15432     5472    1087
     Max         169      88     31                15455     5478    1095
     Sum         673     349    121                61757    21901    4365
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   156

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.43 MB

     Estimated total dynamical RAM >      65.73 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.9

     total cpu time spent up to now is       16.6 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.2886751   0.0000000 )
 
      4 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5478 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   156

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  To be done

     Representation     2      1 modes -A_1  D_1  S_1  To be done

     Representation     3      1 modes -B_1  D_3  S_3  To be done

     Representation     4      1 modes -B_1  D_3  S_3  To be done

     Representation     5      1 modes -B_2  D_4  S_4  To be done

     Representation     6      1 modes -B_2  D_4  S_4  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       :  3m18.16s CPU     3m25.72s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :   208.1 secs   av.it.:   7.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.505E-06

      iter #   2 total cpu time :   211.3 secs   av.it.:  11.1
      thresh= 2.123E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.395E-06

      iter #   3 total cpu time :   214.4 secs   av.it.:  10.5
      thresh= 2.323E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.085E-07

      iter #   4 total cpu time :   217.4 secs   av.it.:  10.6
      thresh= 3.294E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.083E-09

      iter #   5 total cpu time :   220.5 secs   av.it.:  11.0
      thresh= 6.390E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.435E-10

      iter #   6 total cpu time :   223.6 secs   av.it.:  10.7
      thresh= 1.198E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.837E-13

      iter #   7 total cpu time :   226.6 secs   av.it.:  10.8
      thresh= 5.326E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.376E-15

      iter #   8 total cpu time :   229.7 secs   av.it.:  11.1
      thresh= 7.985E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.285E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :   232.2 secs   av.it.:   7.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.412E-05

      iter #   2 total cpu time :   235.3 secs   av.it.:  11.0
      thresh= 4.911E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.410E-05

      iter #   3 total cpu time :   238.2 secs   av.it.:  10.2
      thresh= 5.839E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.162E-08

      iter #   4 total cpu time :   241.4 secs   av.it.:  11.9
      thresh= 2.040E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.744E-08

      iter #   5 total cpu time :   244.5 secs   av.it.:  11.0
      thresh= 1.321E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.799E-11

      iter #   6 total cpu time :   247.6 secs   av.it.:  10.6
      thresh= 6.164E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.219E-12

      iter #   7 total cpu time :   250.6 secs   av.it.:  10.7
      thresh= 1.104E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.006E-14

      iter #   8 total cpu time :   253.8 secs   av.it.:  11.2
      thresh= 1.003E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.111E-15

      iter #   9 total cpu time :   256.8 secs   av.it.:  10.8
      thresh= 7.817E-09 alpha_mix =  0.700 |ddv_scf|^2 =  5.063E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :   259.3 secs   av.it.:   7.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.207E-07

      iter #   2 total cpu time :   262.5 secs   av.it.:  11.7
      thresh= 6.486E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.748E-09

      iter #   3 total cpu time :   265.7 secs   av.it.:  11.3
      thresh= 4.181E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.117E-10

      iter #   4 total cpu time :   268.8 secs   av.it.:  10.8
      thresh= 1.057E-06 alpha_mix =  0.700 |ddv_scf|^2 =  8.786E-13

      iter #   5 total cpu time :   271.8 secs   av.it.:  10.8
      thresh= 9.373E-08 alpha_mix =  0.700 |ddv_scf|^2 =  3.033E-15

      iter #   6 total cpu time :   274.9 secs   av.it.:  11.4
      thresh= 5.508E-09 alpha_mix =  0.700 |ddv_scf|^2 =  8.738E-18

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :   277.2 secs   av.it.:   6.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.805E-07

      iter #   2 total cpu time :   280.0 secs   av.it.:  10.0
      thresh= 7.619E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.043E-08

      iter #   3 total cpu time :   282.6 secs   av.it.:   9.1
      thresh= 2.246E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.915E-11

      iter #   4 total cpu time :   285.6 secs   av.it.:  10.6
      thresh= 7.011E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.612E-12

      iter #   5 total cpu time :   288.4 secs   av.it.:  10.2
      thresh= 1.270E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.737E-15

      iter #   6 total cpu time :   291.6 secs   av.it.:  11.8
      thresh= 4.167E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.037E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :   294.0 secs   av.it.:   6.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.345E-06

      iter #   2 total cpu time :   297.0 secs   av.it.:   9.9
      thresh= 1.532E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.197E-07

      iter #   3 total cpu time :   299.6 secs   av.it.:  10.2
      thresh= 4.687E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.643E-08

      iter #   4 total cpu time :   302.4 secs   av.it.:  10.4
      thresh= 1.282E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.168E-12

      iter #   5 total cpu time :   305.4 secs   av.it.:  10.7
      thresh= 2.484E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.866E-14

      iter #   6 total cpu time :   308.3 secs   av.it.:  10.2
      thresh= 2.206E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.500E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :   310.8 secs   av.it.:   7.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.944E-06

      iter #   2 total cpu time :   313.8 secs   av.it.:  10.9
      thresh= 1.986E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.948E-07

      iter #   3 total cpu time :   316.7 secs   av.it.:  10.5
      thresh= 7.034E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.127E-09

      iter #   4 total cpu time :   319.6 secs   av.it.:  10.2
      thresh= 3.357E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.505E-11

      iter #   5 total cpu time :   322.6 secs   av.it.:  10.5
      thresh= 3.879E-07 alpha_mix =  0.700 |ddv_scf|^2 =  8.532E-14

      iter #   6 total cpu time :   325.3 secs   av.it.:   9.9
      thresh= 2.921E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.590E-16

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    3
     List of q in the star:
          1   0.000000000   0.288675135   0.000000000
          2  -0.250000000  -0.144337567   0.000000000
          3   0.250000000  -0.144337567   0.000000000
     In addition there is the -q list: 
          1   0.000000000  -0.288675135   0.000000000
          2   0.250000000   0.144337567   0.000000000
          3  -0.250000000   0.144337567   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.288675135   0.000000000 ) 

 **************************************************************************
     freq (    1) =       3.397752 [THz] =     113.336815 [cm-1]
     freq (    2) =      12.028724 [THz] =     401.235050 [cm-1]
     freq (    3) =      21.103853 [THz] =     703.948760 [cm-1]
     freq (    4) =      22.160229 [THz] =     739.185663 [cm-1]
     freq (    5) =      38.681211 [THz] =    1290.266323 [cm-1]
     freq (    6) =      44.483700 [THz] =    1483.816517 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =        113.3  [cm-1]   --> B_2  D_4  S_4      
     freq (  2 -  2) =        401.2  [cm-1]   --> B_1  D_3  S_3      
     freq (  3 -  3) =        703.9  [cm-1]   --> A_1  D_1  S_1      
     freq (  4 -  4) =        739.2  [cm-1]   --> B_2  D_4  S_4      
     freq (  5 -  5) =       1290.3  [cm-1]   --> B_1  D_3  S_3      
     freq (  6 -  6) =       1483.8  [cm-1]   --> A_1  D_1  S_1      

     Calculation of q =    0.0000000   0.4330127   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     30                15432     5472    1158
     Max         169      88     31                15455     5478    1173
     Sum         673     349    121                61757    21901    4653
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   156

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.44 MB

     Estimated total dynamical RAM >      65.75 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.8

     total cpu time spent up to now is       24.9 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.4330127   0.0000000 )
 
      4 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   156

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  To be done

     Representation     2      1 modes -A_1  D_1  S_1  To be done

     Representation     3      1 modes -B_1  D_3  S_3  To be done

     Representation     4      1 modes -B_1  D_3  S_3  To be done

     Representation     5      1 modes -B_2  D_4  S_4  To be done

     Representation     6      1 modes -B_2  D_4  S_4  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       :  5m24.57s CPU     5m36.87s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :   339.2 secs   av.it.:   7.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.822E-06

      iter #   2 total cpu time :   342.3 secs   av.it.:  11.2
      thresh= 1.955E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.342E-06

      iter #   3 total cpu time :   345.4 secs   av.it.:  11.2
      thresh= 1.530E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.111E-07

      iter #   4 total cpu time :   348.4 secs   av.it.:  10.9
      thresh= 5.578E-05 alpha_mix =  0.700 |ddv_scf|^2 =  7.404E-09

      iter #   5 total cpu time :   351.5 secs   av.it.:  11.3
      thresh= 8.605E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.451E-10

      iter #   6 total cpu time :   354.6 secs   av.it.:  11.2
      thresh= 1.204E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.947E-13

      iter #   7 total cpu time :   357.7 secs   av.it.:  11.2
      thresh= 6.282E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.795E-14

      iter #   8 total cpu time :   360.8 secs   av.it.:  11.2
      thresh= 1.672E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.086E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :   363.3 secs   av.it.:   7.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.100E-05

      iter #   2 total cpu time :   366.5 secs   av.it.:  11.2
      thresh= 3.316E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.035E-05

      iter #   3 total cpu time :   369.5 secs   av.it.:  10.8
      thresh= 3.217E-04 alpha_mix =  0.700 |ddv_scf|^2 =  7.098E-08

      iter #   4 total cpu time :   372.7 secs   av.it.:  11.5
      thresh= 2.664E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.164E-09

      iter #   5 total cpu time :   375.9 secs   av.it.:  11.3
      thresh= 7.186E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.287E-11

      iter #   6 total cpu time :   379.0 secs   av.it.:  11.0
      thresh= 5.733E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.899E-13

      iter #   7 total cpu time :   382.1 secs   av.it.:  11.1
      thresh= 8.888E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.689E-14

      iter #   8 total cpu time :   385.3 secs   av.it.:  11.2
      thresh= 1.300E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.111E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :   387.7 secs   av.it.:   7.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.941E-07

      iter #   2 total cpu time :   390.9 secs   av.it.:  11.6
      thresh= 7.029E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.822E-09

      iter #   3 total cpu time :   394.0 secs   av.it.:  11.3
      thresh= 5.312E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.261E-10

      iter #   4 total cpu time :   397.0 secs   av.it.:  10.7
      thresh= 1.123E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.660E-13

      iter #   5 total cpu time :   400.0 secs   av.it.:  10.6
      thresh= 8.752E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.343E-15

      iter #   6 total cpu time :   403.1 secs   av.it.:  11.6
      thresh= 3.665E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.291E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :   405.4 secs   av.it.:   6.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  9.545E-07

      iter #   2 total cpu time :   408.3 secs   av.it.:   9.9
      thresh= 9.770E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.345E-08

      iter #   3 total cpu time :   410.9 secs   av.it.:   9.0
      thresh= 2.889E-05 alpha_mix =  0.700 |ddv_scf|^2 =  7.831E-11

      iter #   4 total cpu time :   414.0 secs   av.it.:  10.9
      thresh= 8.849E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.730E-12

      iter #   5 total cpu time :   416.9 secs   av.it.:  10.0
      thresh= 1.315E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.574E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :   419.1 secs   av.it.:   6.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.454E-06

      iter #   2 total cpu time :   421.8 secs   av.it.:   9.3
      thresh= 1.858E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.080E-07

      iter #   3 total cpu time :   424.4 secs   av.it.:   9.2
      thresh= 5.550E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.502E-10

      iter #   4 total cpu time :   427.4 secs   av.it.:  11.0
      thresh= 2.346E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.909E-12

      iter #   5 total cpu time :   430.4 secs   av.it.:  10.8
      thresh= 1.382E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.550E-14

      iter #   6 total cpu time :   433.3 secs   av.it.:  10.3
      thresh= 1.884E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.383E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :   435.7 secs   av.it.:   7.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.380E-06

      iter #   2 total cpu time :   438.8 secs   av.it.:  11.1
      thresh= 1.175E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.579E-08

      iter #   3 total cpu time :   441.8 secs   av.it.:  10.9
      thresh= 2.140E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.132E-10

      iter #   4 total cpu time :   444.6 secs   av.it.:  10.1
      thresh= 2.265E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.828E-12

      iter #   5 total cpu time :   447.6 secs   av.it.:  10.8
      thresh= 1.682E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.102E-14

      iter #   6 total cpu time :   450.4 secs   av.it.:  10.4
      thresh= 2.025E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.807E-16

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    3
     List of q in the star:
          1   0.000000000   0.433012702   0.000000000
          2  -0.375000000  -0.216506351   0.000000000
          3   0.375000000  -0.216506351   0.000000000
     In addition there is the -q list: 
          1   0.000000000  -0.433012702   0.000000000
          2   0.375000000   0.216506351   0.000000000
          3  -0.375000000   0.216506351   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.433012702   0.000000000 ) 

 **************************************************************************
     freq (    1) =       6.968350 [THz] =     232.439121 [cm-1]
     freq (    2) =      15.281803 [THz] =     509.746082 [cm-1]
     freq (    3) =      20.176104 [THz] =     673.002388 [cm-1]
     freq (    4) =      29.465276 [THz] =     982.855798 [cm-1]
     freq (    5) =      37.669586 [THz] =    1256.522129 [cm-1]
     freq (    6) =      41.671357 [THz] =    1390.006850 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =        232.4  [cm-1]   --> B_2  D_4  S_4      
     freq (  2 -  2) =        509.7  [cm-1]   --> B_1  D_3  S_3      
     freq (  3 -  3) =        673.0  [cm-1]   --> B_2  D_4  S_4      
     freq (  4 -  4) =        982.9  [cm-1]   --> A_1  D_1  S_1      
     freq (  5 -  5) =       1256.5  [cm-1]   --> B_1  D_3  S_3      
     freq (  6 -  6) =       1390.0  [cm-1]   --> A_1  D_1  S_1      

     Calculation of q =    0.0000000  -0.5773503   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     30                15432     5472    1196
     Max         169      88     31                15455     5478    1213
     Sum         673     349    121                61757    21901    4811
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=    86
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.0000000   0.0000000   0.0000000), wk =   0.0138889
        k(    2) = (   0.0000000  -0.5773503   0.0000000), wk =   0.0000000
        k(    3) = (   0.0000000   0.0962250   0.0000000), wk =   0.0277778
        k(    4) = (   0.0000000  -0.4811252   0.0000000), wk =   0.0000000
        k(    5) = (   0.0000000   0.1924501   0.0000000), wk =   0.0277778
        k(    6) = (   0.0000000  -0.3849002   0.0000000), wk =   0.0000000
        k(    7) = (   0.0000000   0.2886751   0.0000000), wk =   0.0277778
        k(    8) = (   0.0000000  -0.2886751   0.0000000), wk =   0.0000000
        k(    9) = (   0.0000000   0.3849002   0.0000000), wk =   0.0277778
        k(   10) = (   0.0000000  -0.1924501   0.0000000), wk =   0.0000000
        k(   11) = (   0.0000000   0.4811252   0.0000000), wk =   0.0277778
        k(   12) = (   0.0000000  -0.0962250   0.0000000), wk =   0.0000000
        k(   13) = (   0.0000000  -0.5773503   0.0000000), wk =   0.0138889
        k(   14) = (   0.0000000  -1.1547005   0.0000000), wk =   0.0000000
        k(   15) = (   0.0833333   0.1443376   0.0000000), wk =   0.0555556
        k(   16) = (   0.0833333  -0.4330127   0.0000000), wk =   0.0000000
        k(   17) = (   0.0833333   0.2405626   0.0000000), wk =   0.0555556
        k(   18) = (   0.0833333  -0.3367877   0.0000000), wk =   0.0000000
        k(   19) = (   0.0833333   0.3367877   0.0000000), wk =   0.0555556
        k(   20) = (   0.0833333  -0.2405626   0.0000000), wk =   0.0000000
        k(   21) = (   0.0833333   0.4330127   0.0000000), wk =   0.0555556
        k(   22) = (   0.0833333  -0.1443376   0.0000000), wk =   0.0000000
        k(   23) = (   0.0833333   0.5292377   0.0000000), wk =   0.0555556
        k(   24) = (   0.0833333  -0.0481125   0.0000000), wk =   0.0000000
        k(   25) = (   0.1666667   0.2886751   0.0000000), wk =   0.0555556
        k(   26) = (   0.1666667  -0.2886751   0.0000000), wk =   0.0000000
        k(   27) = (   0.1666667   0.3849002   0.0000000), wk =   0.0555556
        k(   28) = (   0.1666667  -0.1924501   0.0000000), wk =   0.0000000
        k(   29) = (   0.1666667   0.4811252   0.0000000), wk =   0.0555556
        k(   30) = (   0.1666667  -0.0962250   0.0000000), wk =   0.0000000
        k(   31) = (   0.1666667   0.5773503   0.0000000), wk =   0.0277778
        k(   32) = (   0.1666667   0.0000000   0.0000000), wk =   0.0000000
        k(   33) = (   0.2500000   0.4330127   0.0000000), wk =   0.0555556
        k(   34) = (   0.2500000  -0.1443376   0.0000000), wk =   0.0000000
        k(   35) = (   0.2500000   0.5292377   0.0000000), wk =   0.0555556
        k(   36) = (   0.2500000  -0.0481125   0.0000000), wk =   0.0000000
        k(   37) = (   0.3333333   0.5773503   0.0000000), wk =   0.0277778
        k(   38) = (   0.3333333   0.0000000   0.0000000), wk =   0.0000000
        k(   39) = (  -0.0833333  -0.0481125   0.0000000), wk =   0.0555556
        k(   40) = (  -0.0833333  -0.6254628   0.0000000), wk =   0.0000000
        k(   41) = (  -0.1666667  -0.0962250   0.0000000), wk =   0.0555556
        k(   42) = (  -0.1666667  -0.6735753   0.0000000), wk =   0.0000000
        k(   43) = (  -0.2500000  -0.1443376   0.0000000), wk =   0.0555556
        k(   44) = (  -0.2500000  -0.7216878   0.0000000), wk =   0.0000000
        k(   45) = (  -0.3333333  -0.1924501   0.0000000), wk =   0.0555556
        k(   46) = (  -0.3333333  -0.7698004   0.0000000), wk =   0.0000000
        k(   47) = (  -0.4166667  -0.2405626   0.0000000), wk =   0.0555556
        k(   48) = (  -0.4166667  -0.8179129   0.0000000), wk =   0.0000000
        k(   49) = (   0.5000000   0.2886751   0.0000000), wk =   0.0277778
        k(   50) = (   0.5000000  -0.2886751   0.0000000), wk =   0.0000000
        k(   51) = (   0.1666667   0.0000000   0.0000000), wk =   0.0277778
        k(   52) = (   0.1666667  -0.5773503   0.0000000), wk =   0.0000000
        k(   53) = (  -0.1666667  -0.1924501   0.0000000), wk =   0.0555556
        k(   54) = (  -0.1666667  -0.7698004   0.0000000), wk =   0.0000000
        k(   55) = (   0.2500000  -0.0481125   0.0000000), wk =   0.0555556
        k(   56) = (   0.2500000  -0.6254628   0.0000000), wk =   0.0000000
        k(   57) = (  -0.2500000  -0.2405626   0.0000000), wk =   0.0555556
        k(   58) = (  -0.2500000  -0.8179129   0.0000000), wk =   0.0000000
        k(   59) = (   0.3333333  -0.0962250   0.0000000), wk =   0.0555556
        k(   60) = (   0.3333333  -0.6735753   0.0000000), wk =   0.0000000
        k(   61) = (  -0.3333333  -0.2886751   0.0000000), wk =   0.0555556
        k(   62) = (  -0.3333333  -0.8660254   0.0000000), wk =   0.0000000
        k(   63) = (   0.4166667  -0.1443376   0.0000000), wk =   0.0555556
        k(   64) = (   0.4166667  -0.7216878   0.0000000), wk =   0.0000000
        k(   65) = (  -0.4166667  -0.3367877   0.0000000), wk =   0.0555556
        k(   66) = (  -0.4166667  -0.9141379   0.0000000), wk =   0.0000000
        k(   67) = (   0.5000000  -0.1924501   0.0000000), wk =   0.0555556
        k(   68) = (   0.5000000  -0.7698004   0.0000000), wk =   0.0000000
        k(   69) = (   0.3333333   0.0000000   0.0000000), wk =   0.0277778
        k(   70) = (   0.3333333  -0.5773503   0.0000000), wk =   0.0000000
        k(   71) = (  -0.2500000  -0.3367877   0.0000000), wk =   0.0555556
        k(   72) = (  -0.2500000  -0.9141379   0.0000000), wk =   0.0000000
        k(   73) = (   0.4166667  -0.0481125   0.0000000), wk =   0.0555556
        k(   74) = (   0.4166667  -0.6254628   0.0000000), wk =   0.0000000
        k(   75) = (  -0.3333333  -0.3849002   0.0000000), wk =   0.0555556
        k(   76) = (  -0.3333333  -0.9622504   0.0000000), wk =   0.0000000
        k(   77) = (   0.5000000  -0.0962250   0.0000000), wk =   0.0555556
        k(   78) = (   0.5000000  -0.6735753   0.0000000), wk =   0.0000000
        k(   79) = (  -0.4166667  -0.4330127   0.0000000), wk =   0.0555556
        k(   80) = (  -0.4166667  -1.0103630   0.0000000), wk =   0.0000000
        k(   81) = (   0.5000000   0.0000000   0.0000000), wk =   0.0277778
        k(   82) = (   0.5000000  -0.5773503   0.0000000), wk =   0.0000000
        k(   83) = (  -0.3333333  -0.4811252   0.0000000), wk =   0.0555556
        k(   84) = (  -0.3333333  -1.0584755   0.0000000), wk =   0.0000000
        k(   85) = (   0.5833333  -0.0481125   0.0000000), wk =   0.0555556
        k(   86) = (   0.5833333  -0.6254628   0.0000000), wk =   0.0000000

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.44 MB

     Estimated total dynamical RAM >      65.77 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.8

     total cpu time spent up to now is       29.8 secs

     End of band structure calculation

          k = 0.0000 0.0000 0.0000 (  2761 PWs)   bands (ev):

   -23.4472 -11.0291  -7.1657  -7.1657

          k = 0.0000-0.5774 0.0000 (  2728 PWs)   bands (ev):

   -20.3803 -14.6791  -9.9422  -6.7138

          k = 0.0000 0.0962 0.0000 (  2765 PWs)   bands (ev):

   -23.3142 -10.8505  -7.7229  -7.4220

          k = 0.0000-0.4811 0.0000 (  2722 PWs)   bands (ev):

   -20.7431 -14.0934  -9.8009  -7.2728

          k = 0.0000 0.1925 0.0000 (  2756 PWs)   bands (ev):

   -22.9227 -10.3229  -9.1278  -8.0464

          k = 0.0000-0.3849 0.0000 (  2740 PWs)   bands (ev):

   -21.5143 -12.6773  -9.3927  -8.3883

          k = 0.0000 0.2887 0.0000 (  2740 PWs)   bands (ev):

   -22.2997 -10.9064  -9.4764  -8.7705

          k = 0.0000-0.2887 0.0000 (  2740 PWs)   bands (ev):

   -22.2997 -10.9064  -9.4764  -8.7705

          k = 0.0000 0.3849 0.0000 (  2740 PWs)   bands (ev):

   -21.5143 -12.6773  -9.3927  -8.3883

          k = 0.0000-0.1925 0.0000 (  2756 PWs)   bands (ev):

   -22.9227 -10.3229  -9.1278  -8.0464

          k = 0.0000 0.4811 0.0000 (  2722 PWs)   bands (ev):

   -20.7431 -14.0934  -9.8009  -7.2728

          k = 0.0000-0.0962 0.0000 (  2765 PWs)   bands (ev):

   -23.3142 -10.8505  -7.7229  -7.4220

          k = 0.0000-0.5774 0.0000 (  2728 PWs)   bands (ev):

   -20.3803 -14.6791  -9.9422  -6.7138

          k = 0.0000-1.1547 0.0000 (  2761 PWs)   bands (ev):

   -23.4472 -11.0291  -7.1657  -7.1657

          k = 0.0833 0.1443 0.0000 (  2755 PWs)   bands (ev):

   -23.0518 -10.4973  -8.6365  -7.9300

          k = 0.0833-0.4330 0.0000 (  2731 PWs)   bands (ev):

   -21.0305 -13.4098 -10.0082  -7.6996

          k = 0.0833 0.2406 0.0000 (  2755 PWs)   bands (ev):

   -22.5433 -10.1178  -9.8087  -8.7163

          k = 0.0833-0.3368 0.0000 (  2741 PWs)   bands (ev):

   -21.8326 -11.8280  -9.4595  -8.8338

          k = 0.0833 0.3368 0.0000 (  2741 PWs)   bands (ev):

   -21.8326 -11.8280  -9.4595  -8.8338

          k = 0.0833-0.2406 0.0000 (  2755 PWs)   bands (ev):

   -22.5433 -10.1178  -9.8087  -8.7163

          k = 0.0833 0.4330 0.0000 (  2731 PWs)   bands (ev):

   -21.0305 -13.4098 -10.0082  -7.6996

          k = 0.0833-0.1443 0.0000 (  2755 PWs)   bands (ev):

   -23.0518 -10.4973  -8.6365  -7.9300

          k = 0.0833 0.5292 0.0000 (  2721 PWs)   bands (ev):

   -20.4115 -14.4484 -10.2958  -6.7702

          k = 0.0833-0.0481 0.0000 (  2765 PWs)   bands (ev):

   -23.3142 -10.8505  -7.7229  -7.4220

          k = 0.1667 0.2887 0.0000 (  2745 PWs)   bands (ev):

   -21.9428 -11.3255  -9.7309  -8.9871

          k = 0.1667-0.2887 0.0000 (  2745 PWs)   bands (ev):

   -21.9428 -11.3255  -9.7309  -8.9871

          k = 0.1667 0.3849 0.0000 (  2725 PWs)   bands (ev):

   -21.1936 -12.6731 -10.5949  -7.9384

          k = 0.1667-0.1925 0.0000 (  2755 PWs)   bands (ev):

   -22.5433 -10.1178  -9.8087  -8.7163

          k = 0.1667 0.4811 0.0000 (  2726 PWs)   bands (ev):

   -20.4713 -13.8567 -11.1142  -6.8761

          k = 0.1667-0.0962 0.0000 (  2756 PWs)   bands (ev):

   -22.9227 -10.3229  -9.1278  -8.0464

          k = 0.1667 0.5774 0.0000 (  2724 PWs)   bands (ev):

   -20.1390 -14.3603 -11.2851  -6.3513

          k = 0.1667 0.0000 0.0000 (  2755 PWs)   bands (ev):

   -23.0518 -10.4973  -8.6365  -7.9300

          k = 0.2500 0.4330 0.0000 (  2726 PWs)   bands (ev):

   -20.5001 -13.3969 -11.6794  -6.9261

          k = 0.2500-0.1443 0.0000 (  2740 PWs)   bands (ev):

   -22.2997 -10.9064  -9.4764  -8.7705

          k = 0.2500 0.5292 0.0000 (  2714 PWs)   bands (ev):

   -19.9961 -13.9087 -12.2874  -6.1334

          k = 0.2500-0.0481 0.0000 (  2755 PWs)   bands (ev):

   -22.5433 -10.1178  -9.8087  -8.7163

          k = 0.3333 0.5774 0.0000 (  2736 PWs)   bands (ev):

   -19.8042 -13.6902 -12.9717  -5.8138

          k = 0.3333 0.0000 0.0000 (  2745 PWs)   bands (ev):

   -21.9428 -11.3255  -9.7309  -8.9871

          k =-0.0833-0.0481 0.0000 (  2765 PWs)   bands (ev):

   -23.3142 -10.8505  -7.7229  -7.4220

          k =-0.0833-0.6255 0.0000 (  2721 PWs)   bands (ev):

   -20.4115 -14.4484 -10.2958  -6.7702

          k =-0.1667-0.0962 0.0000 (  2756 PWs)   bands (ev):

   -22.9227 -10.3229  -9.1278  -8.0464

          k =-0.1667-0.6736 0.0000 (  2726 PWs)   bands (ev):

   -20.4713 -13.8567 -11.1142  -6.8761

          k =-0.2500-0.1443 0.0000 (  2740 PWs)   bands (ev):

   -22.2997 -10.9064  -9.4764  -8.7705

          k =-0.2500-0.7217 0.0000 (  2726 PWs)   bands (ev):

   -20.5001 -13.3969 -11.6794  -6.9261

          k =-0.3333-0.1925 0.0000 (  2740 PWs)   bands (ev):

   -21.5143 -12.6773  -9.3927  -8.3883

          k =-0.3333-0.7698 0.0000 (  2726 PWs)   bands (ev):

   -20.4713 -13.8567 -11.1142  -6.8761

          k =-0.4167-0.2406 0.0000 (  2722 PWs)   bands (ev):

   -20.7431 -14.0934  -9.8009  -7.2728

          k =-0.4167-0.8179 0.0000 (  2721 PWs)   bands (ev):

   -20.4115 -14.4484 -10.2958  -6.7702

          k = 0.5000 0.2887 0.0000 (  2728 PWs)   bands (ev):

   -20.3803 -14.6791  -9.9422  -6.7138

          k = 0.5000-0.2887 0.0000 (  2728 PWs)   bands (ev):

   -20.3803 -14.6791  -9.9422  -6.7138

          k = 0.1667 0.0000 0.0000 (  2755 PWs)   bands (ev):

   -23.0518 -10.4973  -8.6365  -7.9300

          k = 0.1667-0.5774 0.0000 (  2724 PWs)   bands (ev):

   -20.1390 -14.3603 -11.2851  -6.3513

          k =-0.1667-0.1925 0.0000 (  2755 PWs)   bands (ev):

   -22.5433 -10.1178  -9.8087  -8.7163

          k =-0.1667-0.7698 0.0000 (  2725 PWs)   bands (ev):

   -21.1936 -12.6731 -10.5949  -7.9384

          k = 0.2500-0.0481 0.0000 (  2755 PWs)   bands (ev):

   -22.5433 -10.1178  -9.8087  -8.7163

          k = 0.2500-0.6255 0.0000 (  2714 PWs)   bands (ev):

   -19.9961 -13.9087 -12.2874  -6.1334

          k =-0.2500-0.2406 0.0000 (  2741 PWs)   bands (ev):

   -21.8326 -11.8280  -9.4595  -8.8338

          k =-0.2500-0.8179 0.0000 (  2725 PWs)   bands (ev):

   -21.1936 -12.6731 -10.5949  -7.9384

          k = 0.3333-0.0962 0.0000 (  2741 PWs)   bands (ev):

   -21.8326 -11.8280  -9.4595  -8.8338

          k = 0.3333-0.6736 0.0000 (  2714 PWs)   bands (ev):

   -19.9961 -13.9087 -12.2874  -6.1334

          k =-0.3333-0.2887 0.0000 (  2731 PWs)   bands (ev):

   -21.0305 -13.4098 -10.0082  -7.6996

          k =-0.3333-0.8660 0.0000 (  2731 PWs)   bands (ev):

   -21.0305 -13.4098 -10.0082  -7.6996

          k = 0.4167-0.1443 0.0000 (  2731 PWs)   bands (ev):

   -21.0305 -13.4098 -10.0082  -7.6996

          k = 0.4167-0.7217 0.0000 (  2724 PWs)   bands (ev):

   -20.1390 -14.3603 -11.2851  -6.3513

          k =-0.4167-0.3368 0.0000 (  2721 PWs)   bands (ev):

   -20.4115 -14.4484 -10.2958  -6.7702

          k =-0.4167-0.9141 0.0000 (  2722 PWs)   bands (ev):

   -20.7431 -14.0934  -9.8009  -7.2728

          k = 0.5000-0.1925 0.0000 (  2721 PWs)   bands (ev):

   -20.4115 -14.4484 -10.2958  -6.7702

          k = 0.5000-0.7698 0.0000 (  2721 PWs)   bands (ev):

   -20.4115 -14.4484 -10.2958  -6.7702

          k = 0.3333 0.0000 0.0000 (  2745 PWs)   bands (ev):

   -21.9428 -11.3255  -9.7309  -8.9871

          k = 0.3333-0.5774 0.0000 (  2736 PWs)   bands (ev):

   -19.8042 -13.6902 -12.9717  -5.8138

          k =-0.2500-0.3368 0.0000 (  2725 PWs)   bands (ev):

   -21.1936 -12.6731 -10.5949  -7.9384

          k =-0.2500-0.9141 0.0000 (  2741 PWs)   bands (ev):

   -21.8326 -11.8280  -9.4595  -8.8338

          k = 0.4167-0.0481 0.0000 (  2725 PWs)   bands (ev):

   -21.1936 -12.6731 -10.5949  -7.9384

          k = 0.4167-0.6255 0.0000 (  2714 PWs)   bands (ev):

   -19.9961 -13.9087 -12.2874  -6.1334

          k =-0.3333-0.3849 0.0000 (  2726 PWs)   bands (ev):

   -20.4713 -13.8567 -11.1142  -6.8761

          k =-0.3333-0.9623 0.0000 (  2740 PWs)   bands (ev):

   -21.5143 -12.6773  -9.3927  -8.3883

          k = 0.5000-0.0962 0.0000 (  2726 PWs)   bands (ev):

   -20.4713 -13.8567 -11.1142  -6.8761

          k = 0.5000-0.6736 0.0000 (  2726 PWs)   bands (ev):

   -20.4713 -13.8567 -11.1142  -6.8761

          k =-0.4167-0.4330 0.0000 (  2724 PWs)   bands (ev):

   -20.1390 -14.3603 -11.2851  -6.3513

          k =-0.4167-1.0104 0.0000 (  2731 PWs)   bands (ev):

   -21.0305 -13.4098 -10.0082  -7.6996

          k = 0.5000 0.0000 0.0000 (  2726 PWs)   bands (ev):

   -20.5001 -13.3969 -11.6794  -6.9261

          k = 0.5000-0.5774 0.0000 (  2726 PWs)   bands (ev):

   -20.5001 -13.3969 -11.6794  -6.9261

          k =-0.3333-0.4811 0.0000 (  2714 PWs)   bands (ev):

   -19.9961 -13.9087 -12.2874  -6.1334

          k =-0.3333-1.0585 0.0000 (  2741 PWs)   bands (ev):

   -21.8326 -11.8280  -9.4595  -8.8338

          k = 0.5833-0.0481 0.0000 (  2714 PWs)   bands (ev):

   -19.9961 -13.9087 -12.2874  -6.1334

          k = 0.5833-0.6255 0.0000 (  2725 PWs)   bands (ev):

   -21.1936 -12.6731 -10.5949  -7.9384

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000  -0.5773503   0.0000000 )
 
      5 Sym.Ops. (with q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=    86

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  To be done

     Representation     2      1 modes -A_1  D_1  S_1  To be done

     Representation     3      1 modes -B_1  D_3  S_3  To be done

     Representation     4      1 modes -B_1  D_3  S_3  To be done

     Representation     5      1 modes -B_2  D_4  S_4  To be done

     Representation     6      1 modes -B_2  D_4  S_4  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       :  7m21.21s CPU     7m38.00s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :   459.4 secs   av.it.:   7.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.544E-06

      iter #   2 total cpu time :   461.2 secs   av.it.:  11.1
      thresh= 2.132E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.438E-06

      iter #   3 total cpu time :   463.0 secs   av.it.:  11.4
      thresh= 1.199E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.662E-07

      iter #   4 total cpu time :   464.8 secs   av.it.:  11.1
      thresh= 7.525E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.036E-08

      iter #   5 total cpu time :   466.7 secs   av.it.:  11.1
      thresh= 1.018E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.688E-10

      iter #   6 total cpu time :   468.5 secs   av.it.:  11.6
      thresh= 1.299E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.321E-13

      iter #   7 total cpu time :   470.3 secs   av.it.:  11.5
      thresh= 7.295E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.536E-14

      iter #   8 total cpu time :   472.2 secs   av.it.:  11.3
      thresh= 2.556E-08 alpha_mix =  0.700 |ddv_scf|^2 =  7.276E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :   473.7 secs   av.it.:   7.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  7.483E-06

      iter #   2 total cpu time :   475.5 secs   av.it.:  11.3
      thresh= 2.735E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.195E-06

      iter #   3 total cpu time :   477.3 secs   av.it.:  11.0
      thresh= 2.048E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.272E-09

      iter #   4 total cpu time :   479.2 secs   av.it.:  11.9
      thresh= 5.720E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.442E-09

      iter #   5 total cpu time :   481.0 secs   av.it.:  11.1
      thresh= 5.866E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.229E-11

      iter #   6 total cpu time :   482.7 secs   av.it.:  10.8
      thresh= 6.503E-07 alpha_mix =  0.700 |ddv_scf|^2 =  8.400E-13

      iter #   7 total cpu time :   484.6 secs   av.it.:  11.4
      thresh= 9.165E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.076E-14

      iter #   8 total cpu time :   486.4 secs   av.it.:  11.4
      thresh= 1.037E-08 alpha_mix =  0.700 |ddv_scf|^2 =  3.271E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :   487.9 secs   av.it.:   7.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.486E-07

      iter #   2 total cpu time :   489.7 secs   av.it.:  11.4
      thresh= 7.407E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.410E-09

      iter #   3 total cpu time :   491.5 secs   av.it.:  11.1
      thresh= 7.356E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.534E-10

      iter #   4 total cpu time :   493.2 secs   av.it.:   9.9
      thresh= 1.880E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.989E-13

      iter #   5 total cpu time :   494.9 secs   av.it.:  10.3
      thresh= 8.360E-08 alpha_mix =  0.700 |ddv_scf|^2 =  3.516E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :   496.2 secs   av.it.:   6.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.443E-06

      iter #   2 total cpu time :   497.9 secs   av.it.:   9.7
      thresh= 1.201E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.270E-07

      iter #   3 total cpu time :   499.4 secs   av.it.:   8.9
      thresh= 3.564E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.840E-11

      iter #   4 total cpu time :   501.2 secs   av.it.:  11.0
      thresh= 9.402E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.762E-12

      iter #   5 total cpu time :   502.9 secs   av.it.:   9.9
      thresh= 1.327E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.316E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :   504.2 secs   av.it.:   6.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.305E-06

      iter #   2 total cpu time :   505.7 secs   av.it.:   8.7
      thresh= 2.303E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.758E-07

      iter #   3 total cpu time :   507.2 secs   av.it.:   8.5
      thresh= 6.898E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.199E-11

      iter #   4 total cpu time :   509.0 secs   av.it.:  11.1
      thresh= 7.211E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.012E-12

      iter #   5 total cpu time :   510.8 secs   av.it.:  10.8
      thresh= 1.006E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.387E-14

      iter #   6 total cpu time :   512.5 secs   av.it.:  10.3
      thresh= 1.840E-08 alpha_mix =  0.700 |ddv_scf|^2 =  5.934E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :   514.0 secs   av.it.:   7.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  9.568E-07

      iter #   2 total cpu time :   515.8 secs   av.it.:  11.0
      thresh= 9.782E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.932E-08

      iter #   3 total cpu time :   517.6 secs   av.it.:  10.9
      thresh= 1.390E-05 alpha_mix =  0.700 |ddv_scf|^2 =  9.720E-10

      iter #   4 total cpu time :   519.2 secs   av.it.:   9.7
      thresh= 3.118E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.796E-13

      iter #   5 total cpu time :   521.0 secs   av.it.:  10.7
      thresh= 9.898E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.727E-14

      iter #   6 total cpu time :   522.7 secs   av.it.:  10.5
      thresh= 1.651E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.444E-16

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    3
     List of q in the star:
          1   0.000000000  -0.577350269   0.000000000
          2   0.500000000   0.288675135   0.000000000
          3  -0.500000000   0.288675135   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000  -0.577350269   0.000000000 ) 

 **************************************************************************
     freq (    1) =       9.063819 [THz] =     302.336454 [cm-1]
     freq (    2) =      16.324556 [THz] =     544.528565 [cm-1]
     freq (    3) =      18.793443 [THz] =     626.881773 [cm-1]
     freq (    4) =      34.469572 [THz] =    1149.781161 [cm-1]
     freq (    5) =      37.283614 [THz] =    1243.647502 [cm-1]
     freq (    6) =      38.511632 [THz] =    1284.609773 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =        302.3  [cm-1]   --> B_2  D_4  S_4      
     freq (  2 -  2) =        544.5  [cm-1]   --> B_1  D_3  S_3      
     freq (  3 -  3) =        626.9  [cm-1]   --> B_2  D_4  S_4      
     freq (  4 -  4) =       1149.8  [cm-1]   --> A_1  D_1  S_1      
     freq (  5 -  5) =       1243.6  [cm-1]   --> B_1  D_3  S_3      
     freq (  6 -  6) =       1284.6  [cm-1]   --> A_1  D_1  S_1      

     Calculation of q =    0.1250000   0.2165064   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     29                15432     5472    1081
     Max         169      88     31                15455     5478    1089
     Sum         673     349    121                61757    21901    4341
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   288

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.43 MB

     Estimated total dynamical RAM >      65.72 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.9

     total cpu time spent up to now is       44.2 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.1250000   0.2165064   0.0000000 )
 
      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5478 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   288

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A'  To be done

     Representation     4      1 modes -A'  To be done

     Representation     5      1 modes -A''  To be done

     Representation     6      1 modes -A''  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       :  8m41.83s CPU     9m 1.02s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :   545.1 secs   av.it.:   7.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.553E-06

      iter #   2 total cpu time :   550.7 secs   av.it.:  11.7
      thresh= 1.246E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.593E-06

      iter #   3 total cpu time :   555.7 secs   av.it.:  10.4
      thresh= 1.262E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.249E-08

      iter #   4 total cpu time :   561.2 secs   av.it.:  11.4
      thresh= 1.118E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.189E-09

      iter #   5 total cpu time :   566.7 secs   av.it.:  11.5
      thresh= 4.679E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.487E-11

      iter #   6 total cpu time :   572.3 secs   av.it.:  11.5
      thresh= 4.987E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.656E-13

      iter #   7 total cpu time :   577.7 secs   av.it.:  11.3
      thresh= 6.824E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.106E-14

      iter #   8 total cpu time :   583.2 secs   av.it.:  11.5
      thresh= 1.052E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.623E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :   587.2 secs   av.it.:   7.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.480E-06

      iter #   2 total cpu time :   592.6 secs   av.it.:  11.2
      thresh= 1.865E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.151E-06

      iter #   3 total cpu time :   597.7 secs   av.it.:  10.5
      thresh= 2.037E-04 alpha_mix =  0.700 |ddv_scf|^2 =  7.129E-08

      iter #   4 total cpu time :   602.8 secs   av.it.:  10.5
      thresh= 2.670E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.451E-09

      iter #   5 total cpu time :   608.2 secs   av.it.:  11.2
      thresh= 3.810E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.180E-11

      iter #   6 total cpu time :   613.5 secs   av.it.:  11.1
      thresh= 9.581E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.549E-13

      iter #   7 total cpu time :   618.8 secs   av.it.:  11.2
      thresh= 5.049E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.730E-14

      iter #   8 total cpu time :   624.1 secs   av.it.:  11.1
      thresh= 1.315E-08 alpha_mix =  0.700 |ddv_scf|^2 =  3.670E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :   627.9 secs   av.it.:   6.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.898E-06

      iter #   2 total cpu time :   633.2 secs   av.it.:  11.2
      thresh= 1.378E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.951E-06

      iter #   3 total cpu time :   638.3 secs   av.it.:  10.5
      thresh= 1.397E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.899E-08

      iter #   4 total cpu time :   643.5 secs   av.it.:  10.7
      thresh= 2.429E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.435E-09

      iter #   5 total cpu time :   648.8 secs   av.it.:  11.3
      thresh= 4.935E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.353E-11

      iter #   6 total cpu time :   654.0 secs   av.it.:  11.2
      thresh= 4.850E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.585E-13

      iter #   7 total cpu time :   659.2 secs   av.it.:  11.2
      thresh= 5.988E-08 alpha_mix =  0.700 |ddv_scf|^2 =  7.720E-15

      iter #   8 total cpu time :   664.4 secs   av.it.:  11.3
      thresh= 8.786E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.267E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :   668.6 secs   av.it.:   7.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.852E-05

      iter #   2 total cpu time :   673.8 secs   av.it.:  11.1
      thresh= 5.340E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.996E-05

      iter #   3 total cpu time :   678.7 secs   av.it.:  10.3
      thresh= 6.322E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.955E-08

      iter #   4 total cpu time :   684.4 secs   av.it.:  12.1
      thresh= 1.719E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.189E-08

      iter #   5 total cpu time :   689.7 secs   av.it.:  11.2
      thresh= 1.479E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.540E-11

      iter #   6 total cpu time :   695.0 secs   av.it.:  10.8
      thresh= 6.738E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.089E-12

      iter #   7 total cpu time :   700.1 secs   av.it.:  10.6
      thresh= 1.445E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.429E-14

      iter #   8 total cpu time :   705.3 secs   av.it.:  10.9
      thresh= 1.195E-08 alpha_mix =  0.700 |ddv_scf|^2 =  7.168E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :   709.3 secs   av.it.:   6.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.319E-06

      iter #   2 total cpu time :   714.3 secs   av.it.:  10.2
      thresh= 1.523E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.517E-07

      iter #   3 total cpu time :   719.5 secs   av.it.:  10.5
      thresh= 5.017E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.663E-08

      iter #   4 total cpu time :   724.5 secs   av.it.:  10.3
      thresh= 1.632E-05 alpha_mix =  0.700 |ddv_scf|^2 =  9.236E-12

      iter #   5 total cpu time :   729.7 secs   av.it.:  10.7
      thresh= 3.039E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.348E-14

      iter #   6 total cpu time :   734.8 secs   av.it.:  10.3
      thresh= 2.519E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.013E-15

      iter #   7 total cpu time :   739.8 secs   av.it.:  10.1
      thresh= 3.182E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.409E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :   744.1 secs   av.it.:   7.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.665E-06

      iter #   2 total cpu time :   749.4 secs   av.it.:  11.0
      thresh= 2.380E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.014E-06

      iter #   3 total cpu time :   754.5 secs   av.it.:  10.5
      thresh= 1.007E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.858E-09

      iter #   4 total cpu time :   759.6 secs   av.it.:  10.3
      thresh= 4.311E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.070E-11

      iter #   5 total cpu time :   764.8 secs   av.it.:  10.6
      thresh= 4.550E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.123E-13

      iter #   6 total cpu time :   769.6 secs   av.it.:   9.8
      thresh= 3.351E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.580E-15

      iter #   7 total cpu time :   774.7 secs   av.it.:  10.7
      thresh= 3.975E-09 alpha_mix =  0.700 |ddv_scf|^2 =  5.545E-18

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    6
     List of q in the star:
          1   0.125000000   0.216506351   0.000000000
          2  -0.125000000   0.216506351   0.000000000
          3  -0.125000000  -0.216506351   0.000000000
          4   0.250000000   0.000000000   0.000000000
          5  -0.250000000   0.000000000   0.000000000
          6   0.125000000  -0.216506351   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.125000000   0.216506351   0.000000000 ) 

 **************************************************************************
     freq (    1) =       2.476612 [THz] =      82.610876 [cm-1]
     freq (    2) =      11.829128 [THz] =     394.577225 [cm-1]
     freq (    3) =      17.914426 [THz] =     597.560933 [cm-1]
     freq (    4) =      22.582156 [THz] =     753.259645 [cm-1]
     freq (    5) =      38.999947 [THz] =    1300.898196 [cm-1]
     freq (    6) =      44.867534 [THz] =    1496.619835 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =         82.6  [cm-1]   --> A''                
     freq (  2 -  2) =        394.6  [cm-1]   --> A'                 
     freq (  3 -  3) =        597.6  [cm-1]   --> A'                 
     freq (  4 -  4) =        753.3  [cm-1]   --> A''                
     freq (  5 -  5) =       1300.9  [cm-1]   --> A'                 
     freq (  6 -  6) =       1496.6  [cm-1]   --> A'                 

     Calculation of q =    0.1250000   0.3608439   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     30                15432     5472    1158
     Max         169      88     31                15455     5478    1173
     Sum         673     349    121                61757    21901    4653
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   288

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.43 MB

     Estimated total dynamical RAM >      65.74 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.9

     total cpu time spent up to now is       58.6 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.1250000   0.3608439   0.0000000 )
 
      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   288

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A'  To be done

     Representation     4      1 modes -A'  To be done

     Representation     5      1 modes -A''  To be done

     Representation     6      1 modes -A''  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       : 12m44.37s CPU    13m12.96s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :   797.1 secs   av.it.:   7.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.047E-07

      iter #   2 total cpu time :   802.7 secs   av.it.:  12.1
      thresh= 8.970E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.479E-07

      iter #   3 total cpu time :   808.1 secs   av.it.:  11.2
      thresh= 4.979E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.453E-08

      iter #   4 total cpu time :   813.4 secs   av.it.:  11.1
      thresh= 1.858E-05 alpha_mix =  0.700 |ddv_scf|^2 =  9.651E-10

      iter #   5 total cpu time :   819.0 secs   av.it.:  11.9
      thresh= 3.107E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.848E-11

      iter #   6 total cpu time :   824.6 secs   av.it.:  11.6
      thresh= 9.924E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.884E-13

      iter #   7 total cpu time :   829.9 secs   av.it.:  11.1
      thresh= 8.879E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.002E-14

      iter #   8 total cpu time :   835.2 secs   av.it.:  11.1
      thresh= 1.415E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.209E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :   839.3 secs   av.it.:   7.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.440E-06

      iter #   2 total cpu time :   844.7 secs   av.it.:  11.4
      thresh= 1.855E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.785E-06

      iter #   3 total cpu time :   850.0 secs   av.it.:  11.2
      thresh= 1.669E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.053E-07

      iter #   4 total cpu time :   855.3 secs   av.it.:  10.9
      thresh= 4.531E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.766E-09

      iter #   5 total cpu time :   861.0 secs   av.it.:  11.5
      thresh= 6.137E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.065E-10

      iter #   6 total cpu time :   866.4 secs   av.it.:  11.4
      thresh= 1.032E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.820E-13

      iter #   7 total cpu time :   871.8 secs   av.it.:  11.2
      thresh= 6.943E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.986E-14

      iter #   8 total cpu time :   877.1 secs   av.it.:  11.3
      thresh= 1.728E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.752E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :   881.1 secs   av.it.:   6.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.295E-06

      iter #   2 total cpu time :   886.4 secs   av.it.:  11.2
      thresh= 1.138E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.387E-07

      iter #   3 total cpu time :   891.7 secs   av.it.:  11.3
      thresh= 7.339E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.602E-07

      iter #   4 total cpu time :   897.2 secs   av.it.:  11.0
      thresh= 4.002E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.169E-09

      iter #   5 total cpu time :   902.9 secs   av.it.:  11.7
      thresh= 4.658E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.424E-11

      iter #   6 total cpu time :   908.7 secs   av.it.:  11.8
      thresh= 7.365E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.687E-13

      iter #   7 total cpu time :   914.2 secs   av.it.:  11.3
      thresh= 7.541E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.596E-14

      iter #   8 total cpu time :   919.7 secs   av.it.:  11.0
      thresh= 1.611E-08 alpha_mix =  0.700 |ddv_scf|^2 =  3.290E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :   924.1 secs   av.it.:   8.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.364E-05

      iter #   2 total cpu time :   929.7 secs   av.it.:  11.4
      thresh= 3.694E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.535E-05

      iter #   3 total cpu time :   935.0 secs   av.it.:  10.8
      thresh= 3.918E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.397E-08

      iter #   4 total cpu time :   940.6 secs   av.it.:  11.9
      thresh= 2.529E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.450E-09

      iter #   5 total cpu time :   946.1 secs   av.it.:  11.3
      thresh= 9.192E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.125E-11

      iter #   6 total cpu time :   951.7 secs   av.it.:  11.2
      thresh= 7.826E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.521E-12

      iter #   7 total cpu time :   957.2 secs   av.it.:  11.2
      thresh= 1.588E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.830E-14

      iter #   8 total cpu time :   962.7 secs   av.it.:  11.3
      thresh= 1.957E-08 alpha_mix =  0.700 |ddv_scf|^2 =  8.799E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :   966.7 secs   av.it.:   6.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.937E-06

      iter #   2 total cpu time :   971.7 secs   av.it.:   9.9
      thresh= 1.714E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.604E-07

      iter #   3 total cpu time :   976.5 secs   av.it.:   9.8
      thresh= 5.103E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.132E-09

      iter #   4 total cpu time :   981.7 secs   av.it.:  11.0
      thresh= 4.617E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.680E-12

      iter #   5 total cpu time :   987.0 secs   av.it.:  11.0
      thresh= 1.637E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.589E-14

      iter #   6 total cpu time :   992.1 secs   av.it.:  10.6
      thresh= 1.894E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.409E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :   996.4 secs   av.it.:   7.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.883E-06

      iter #   2 total cpu time :  1001.7 secs   av.it.:  11.3
      thresh= 1.372E-04 alpha_mix =  0.700 |ddv_scf|^2 =  9.805E-08

      iter #   3 total cpu time :  1006.8 secs   av.it.:  11.1
      thresh= 3.131E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.913E-10

      iter #   4 total cpu time :  1010.6 secs   av.it.:  10.4
      thresh= 2.432E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.424E-12

      iter #   5 total cpu time :  1015.4 secs   av.it.:  10.8
      thresh= 2.329E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.980E-14

      iter #   6 total cpu time :  1020.3 secs   av.it.:  10.4
      thresh= 2.232E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.059E-16

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    6
     List of q in the star:
          1   0.125000000   0.360843918   0.000000000
          2  -0.125000000   0.360843918   0.000000000
          3  -0.250000000  -0.288675135   0.000000000
          4   0.375000000  -0.072168784   0.000000000
          5  -0.375000000  -0.072168784   0.000000000
          6   0.250000000  -0.288675135   0.000000000
     In addition there is the -q list: 
          1  -0.125000000  -0.360843918   0.000000000
          2   0.125000000  -0.360843918   0.000000000
          3   0.250000000   0.288675135   0.000000000
          4  -0.375000000   0.072168784   0.000000000
          5   0.375000000   0.072168784   0.000000000
          6  -0.250000000   0.288675135   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.125000000   0.360843918   0.000000000 ) 

 **************************************************************************
     freq (    1) =       5.508775 [THz] =     183.752959 [cm-1]
     freq (    2) =      16.546000 [THz] =     551.915141 [cm-1]
     freq (    3) =      20.906780 [THz] =     697.375120 [cm-1]
     freq (    4) =      25.516590 [THz] =     851.141811 [cm-1]
     freq (    5) =      38.017912 [THz] =    1268.141044 [cm-1]
     freq (    6) =      42.621795 [THz] =    1421.710063 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =        183.8  [cm-1]   --> A''                
     freq (  2 -  2) =        551.9  [cm-1]   --> A'                 
     freq (  3 -  3) =        697.4  [cm-1]   --> A''                
     freq (  4 -  4) =        851.1  [cm-1]   --> A'                 
     freq (  5 -  5) =       1268.1  [cm-1]   --> A'                 
     freq (  6 -  6) =       1421.7  [cm-1]   --> A'                 

     Calculation of q =    0.1250000   0.5051815   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     30                15432     5472    1206
     Max         169      88     31                15455     5478    1225
     Sum         673     349    121                61757    21901    4847
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   288

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.44 MB

     Estimated total dynamical RAM >      65.77 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.9

     total cpu time spent up to now is       71.1 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.1250000   0.5051815   0.0000000 )
 
      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   288

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A'  To be done

     Representation     4      1 modes -A'  To be done

     Representation     5      1 modes -A''  To be done

     Representation     6      1 modes -A''  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       : 16m38.78s CPU    17m15.80s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :  1039.8 secs   av.it.:   7.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.260E-06

      iter #   2 total cpu time :  1045.3 secs   av.it.:  12.0
      thresh= 1.123E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.503E-07

      iter #   3 total cpu time :  1050.6 secs   av.it.:  11.4
      thresh= 5.919E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.425E-08

      iter #   4 total cpu time :  1054.9 secs   av.it.:  11.3
      thresh= 1.194E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.810E-10

      iter #   5 total cpu time :  1060.4 secs   av.it.:  12.0
      thresh= 1.952E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.531E-10

      iter #   6 total cpu time :  1065.8 secs   av.it.:  11.6
      thresh= 1.237E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.533E-12

      iter #   7 total cpu time :  1071.2 secs   av.it.:  11.5
      thresh= 1.591E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.453E-14

      iter #   8 total cpu time :  1076.3 secs   av.it.:  11.0
      thresh= 1.205E-08 alpha_mix =  0.700 |ddv_scf|^2 =  7.515E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :  1080.3 secs   av.it.:   7.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.852E-06

      iter #   2 total cpu time :  1085.6 secs   av.it.:  11.4
      thresh= 1.963E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.463E-06

      iter #   3 total cpu time :  1091.0 secs   av.it.:  11.6
      thresh= 1.210E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.545E-07

      iter #   4 total cpu time :  1096.2 secs   av.it.:  11.3
      thresh= 6.742E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.784E-09

      iter #   5 total cpu time :  1100.4 secs   av.it.:  11.4
      thresh= 8.236E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.255E-10

      iter #   6 total cpu time :  1104.6 secs   av.it.:  11.8
      thresh= 1.120E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.990E-13

      iter #   7 total cpu time :  1108.7 secs   av.it.:  11.3
      thresh= 7.739E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.408E-14

      iter #   8 total cpu time :  1112.5 secs   av.it.:  11.3
      thresh= 2.531E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.198E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :  1115.3 secs   av.it.:   7.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.509E-06

      iter #   2 total cpu time :  1119.1 secs   av.it.:  10.9
      thresh= 1.228E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.034E-07

      iter #   3 total cpu time :  1123.2 secs   av.it.:  11.8
      thresh= 4.510E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.618E-07

      iter #   4 total cpu time :  1127.2 secs   av.it.:  11.3
      thresh= 4.023E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.139E-09

      iter #   5 total cpu time :  1131.4 secs   av.it.:  11.8
      thresh= 4.625E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.746E-11

      iter #   6 total cpu time :  1135.5 secs   av.it.:  11.9
      thresh= 9.872E-07 alpha_mix =  0.700 |ddv_scf|^2 =  9.037E-13

      iter #   7 total cpu time :  1139.5 secs   av.it.:  11.4
      thresh= 9.507E-08 alpha_mix =  0.700 |ddv_scf|^2 =  3.571E-14

      iter #   8 total cpu time :  1143.4 secs   av.it.:  10.9
      thresh= 1.890E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.460E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :  1147.6 secs   av.it.:   7.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.825E-06

      iter #   2 total cpu time :  1153.1 secs   av.it.:  11.7
      thresh= 2.613E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.122E-06

      iter #   3 total cpu time :  1158.3 secs   av.it.:  11.3
      thresh= 2.030E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.502E-08

      iter #   4 total cpu time :  1163.7 secs   av.it.:  11.7
      thresh= 1.871E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.314E-09

      iter #   5 total cpu time :  1169.1 secs   av.it.:  11.5
      thresh= 6.568E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.974E-11

      iter #   6 total cpu time :  1174.4 secs   av.it.:  11.5
      thresh= 7.053E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.658E-12

      iter #   7 total cpu time :  1179.7 secs   av.it.:  11.5
      thresh= 1.288E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.014E-14

      iter #   8 total cpu time :  1185.1 secs   av.it.:  11.6
      thresh= 1.419E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.050E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :  1188.8 secs   av.it.:   6.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.547E-06

      iter #   2 total cpu time :  1193.3 secs   av.it.:   9.1
      thresh= 2.132E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.079E-07

      iter #   3 total cpu time :  1197.6 secs   av.it.:   8.9
      thresh= 6.387E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.656E-11

      iter #   4 total cpu time :  1202.9 secs   av.it.:  11.3
      thresh= 8.158E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.123E-12

      iter #   5 total cpu time :  1208.0 secs   av.it.:  10.9
      thresh= 1.060E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.389E-14

      iter #   6 total cpu time :  1213.0 secs   av.it.:  10.5
      thresh= 1.841E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.346E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :  1217.1 secs   av.it.:   7.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  9.453E-07

      iter #   2 total cpu time :  1221.7 secs   av.it.:  11.3
      thresh= 9.723E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.722E-08

      iter #   3 total cpu time :  1225.7 secs   av.it.:  11.2
      thresh= 1.312E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.701E-10

      iter #   4 total cpu time :  1229.4 secs   av.it.:  10.2
      thresh= 2.168E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.189E-12

      iter #   5 total cpu time :  1233.2 secs   av.it.:  10.8
      thresh= 1.090E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.001E-14

      iter #   6 total cpu time :  1237.1 secs   av.it.:  10.8
      thresh= 1.732E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.483E-16

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    6
     List of q in the star:
          1   0.125000000   0.505181486   0.000000000
          2  -0.125000000   0.505181486   0.000000000
          3  -0.375000000  -0.360843918   0.000000000
          4   0.500000000  -0.144337567   0.000000000
          5  -0.500000000  -0.144337567   0.000000000
          6   0.375000000  -0.360843918   0.000000000
     In addition there is the -q list: 
          1  -0.125000000  -0.505181486   0.000000000
          2   0.125000000  -0.505181486   0.000000000
          3   0.375000000   0.360843918   0.000000000
          4  -0.500000000   0.144337567   0.000000000
          5   0.500000000   0.144337567   0.000000000
          6  -0.375000000   0.360843918   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.125000000   0.505181486   0.000000000 ) 

 **************************************************************************
     freq (    1) =       8.481405 [THz] =     282.909228 [cm-1]
     freq (    2) =      18.820390 [THz] =     627.780620 [cm-1]
     freq (    3) =      18.942441 [THz] =     631.851814 [cm-1]
     freq (    4) =      32.081582 [THz] =    1070.126398 [cm-1]
     freq (    5) =      37.460261 [THz] =    1249.539822 [cm-1]
     freq (    6) =      39.000160 [THz] =    1300.905303 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =        282.9  [cm-1]   --> A''                
     freq (  2 -  2) =        627.8  [cm-1]   --> A'                 
     freq (  3 -  3) =        631.9  [cm-1]   --> A''                
     freq (  4 -  4) =       1070.1  [cm-1]   --> A'                 
     freq (  5 -  5) =       1249.5  [cm-1]   --> A'                 
     freq (  6 -  6) =       1300.9  [cm-1]   --> A'                 

     Calculation of q =    0.2500000   0.4330127   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     30                15432     5472    1202
     Max         169      88     31                15455     5478    1221
     Sum         673     349    121                61757    21901    4835
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   288

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.44 MB

     Estimated total dynamical RAM >      65.76 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.9

     total cpu time spent up to now is       81.9 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.2500000   0.4330127   0.0000000 )
 
      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   288

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A'  To be done

     Representation     4      1 modes -A'  To be done

     Representation     5      1 modes -A''  To be done

     Representation     6      1 modes -A''  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       : 20m 5.01s CPU    20m51.20s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :  1254.2 secs   av.it.:   7.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.575E-06

      iter #   2 total cpu time :  1258.4 secs   av.it.:  12.0
      thresh= 1.255E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.396E-07

      iter #   3 total cpu time :  1262.4 secs   av.it.:  11.4
      thresh= 6.630E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.504E-08

      iter #   4 total cpu time :  1266.5 secs   av.it.:  11.5
      thresh= 1.582E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.500E-09

      iter #   5 total cpu time :  1270.6 secs   av.it.:  11.7
      thresh= 3.874E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.787E-11

      iter #   6 total cpu time :  1274.7 secs   av.it.:  11.6
      thresh= 9.893E-07 alpha_mix =  0.700 |ddv_scf|^2 =  9.770E-13

      iter #   7 total cpu time :  1279.5 secs   av.it.:  11.1
      thresh= 9.884E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.023E-14

      iter #   8 total cpu time :  1284.8 secs   av.it.:  11.3
      thresh= 1.422E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.199E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :  1288.8 secs   av.it.:   7.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.152E-06

      iter #   2 total cpu time :  1294.1 secs   av.it.:  11.4
      thresh= 1.776E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.272E-06

      iter #   3 total cpu time :  1299.4 secs   av.it.:  11.5
      thresh= 1.128E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.766E-07

      iter #   4 total cpu time :  1304.6 secs   av.it.:  11.2
      thresh= 6.137E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.068E-09

      iter #   5 total cpu time :  1309.2 secs   av.it.:  11.7
      thresh= 4.548E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.334E-11

      iter #   6 total cpu time :  1313.2 secs   av.it.:  11.7
      thresh= 8.564E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.445E-12

      iter #   7 total cpu time :  1317.4 secs   av.it.:  11.3
      thresh= 1.202E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.776E-14

      iter #   8 total cpu time :  1322.6 secs   av.it.:  11.1
      thresh= 2.603E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.472E-15

      iter #   9 total cpu time :  1327.8 secs   av.it.:  11.4
      thresh= 3.836E-09 alpha_mix =  0.700 |ddv_scf|^2 =  8.117E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :  1331.7 secs   av.it.:   7.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.037E-06

      iter #   2 total cpu time :  1337.0 secs   av.it.:  11.3
      thresh= 1.427E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.398E-07

      iter #   3 total cpu time :  1342.4 secs   av.it.:  11.7
      thresh= 7.347E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.899E-07

      iter #   4 total cpu time :  1347.6 secs   av.it.:  11.3
      thresh= 5.384E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.693E-09

      iter #   5 total cpu time :  1353.0 secs   av.it.:  11.6
      thresh= 6.850E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.824E-11

      iter #   6 total cpu time :  1358.4 secs   av.it.:  11.8
      thresh= 9.912E-07 alpha_mix =  0.700 |ddv_scf|^2 =  9.962E-13

      iter #   7 total cpu time :  1363.6 secs   av.it.:  11.2
      thresh= 9.981E-08 alpha_mix =  0.700 |ddv_scf|^2 =  5.144E-14

      iter #   8 total cpu time :  1368.8 secs   av.it.:  11.0
      thresh= 2.268E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.759E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :  1373.0 secs   av.it.:   8.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.584E-06

      iter #   2 total cpu time :  1378.4 secs   av.it.:  11.7
      thresh= 2.566E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.232E-06

      iter #   3 total cpu time :  1383.6 secs   av.it.:  11.3
      thresh= 2.057E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.338E-08

      iter #   4 total cpu time :  1387.7 secs   av.it.:  11.7
      thresh= 2.311E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.700E-09

      iter #   5 total cpu time :  1391.8 secs   av.it.:  11.5
      thresh= 6.856E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.950E-11

      iter #   6 total cpu time :  1395.8 secs   av.it.:  11.5
      thresh= 8.336E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.366E-12

      iter #   7 total cpu time :  1399.8 secs   av.it.:  11.5
      thresh= 1.538E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.818E-14

      iter #   8 total cpu time :  1403.9 secs   av.it.:  11.7
      thresh= 1.954E-08 alpha_mix =  0.700 |ddv_scf|^2 =  8.267E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :  1406.7 secs   av.it.:   6.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.269E-06

      iter #   2 total cpu time :  1410.5 secs   av.it.:   9.1
      thresh= 2.066E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.828E-07

      iter #   3 total cpu time :  1414.9 secs   av.it.:   8.8
      thresh= 6.187E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.623E-11

      iter #   4 total cpu time :  1420.1 secs   av.it.:  11.2
      thresh= 9.286E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.186E-12

      iter #   5 total cpu time :  1425.3 secs   av.it.:  10.9
      thresh= 1.089E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.242E-14

      iter #   6 total cpu time :  1430.2 secs   av.it.:  10.5
      thresh= 1.801E-08 alpha_mix =  0.700 |ddv_scf|^2 =  6.961E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :  1434.3 secs   av.it.:   7.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  9.816E-07

      iter #   2 total cpu time :  1439.6 secs   av.it.:  11.3
      thresh= 9.907E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.911E-08

      iter #   3 total cpu time :  1444.9 secs   av.it.:  11.2
      thresh= 1.382E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.732E-10

      iter #   4 total cpu time :  1449.8 secs   av.it.:  10.2
      thresh= 2.175E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.303E-12

      iter #   5 total cpu time :  1454.4 secs   av.it.:  10.8
      thresh= 1.142E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.188E-14

      iter #   6 total cpu time :  1457.3 secs   av.it.:  10.8
      thresh= 1.785E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.602E-16

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    6
     List of q in the star:
          1   0.250000000   0.433012702   0.000000000
          2  -0.250000000   0.433012702   0.000000000
          3  -0.250000000  -0.433012702   0.000000000
          4   0.500000000   0.000000000   0.000000000
          5  -0.500000000   0.000000000   0.000000000
          6   0.250000000  -0.433012702   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.250000000   0.433012702   0.000000000 ) 

 **************************************************************************
     freq (    1) =       7.885603 [THz] =     263.035394 [cm-1]
     freq (    2) =      19.074108 [THz] =     636.243745 [cm-1]
     freq (    3) =      22.066814 [THz] =     736.069671 [cm-1]
     freq (    4) =      29.121501 [THz] =     971.388711 [cm-1]
     freq (    5) =      37.592979 [THz] =    1253.966813 [cm-1]
     freq (    6) =      39.195715 [THz] =    1307.428330 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =        263.0  [cm-1]   --> A''                
     freq (  2 -  2) =        636.2  [cm-1]   --> A''                
     freq (  3 -  3) =        736.1  [cm-1]   --> A'                 
     freq (  4 -  4) =        971.4  [cm-1]   --> A'                 
     freq (  5 -  5) =       1254.0  [cm-1]   --> A'                 
     freq (  6 -  6) =       1307.4  [cm-1]   --> A'                 

     Calculation of q =    0.2500000   0.5773503   0.0000000
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         168      86     30                15432     5472    1250
     Max         169      88     31                15455     5478    1271
     Sum         673     349    121                61757    21901    5027
 

     Title: 
     Phonon dispersions for BN                                                  


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     material density (g/cm^3) =       0.7004
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      55.0000  Ry
     charge density cutoff     =     440.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   4.741900  celldm(2)=   0.000000  celldm(3)=   4.300000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     new unit-cell volume =    397.06011 a.u.^3 (    58.83824 Ang^3 )
     density =      0.70038 g/cm^3

CELL_PARAMETERS (alat=  4.74190000)
   1.000000000   0.000000000   0.000000000
  -0.500000000   0.866025404   0.000000000
   0.000000000   0.000000000   4.300000000

ATOMIC_POSITIONS
B        0.000000000   0.288675135   0.000000000
N        0.000000000  -0.288675135   0.000000000


     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   4.300000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.232558 )  


     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        B              3.00    10.81000      B( 1.00)
        N              5.00    14.00674      N( 1.00)

     12 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           B   tau(   1) = (   0.0000000   0.2886751   0.0000000  )
         2           N   tau(   2) = (   0.0000000  -0.2886751   0.0000000  )

     number of k points=   288

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:    61757 G-vectors     FFT dimensions: (  32,  32, 144)

     Smooth grid:    21901 G-vectors     FFT dimensions: (  24,  24, 100)

     Estimated max dynamical RAM per process >      16.45 MB

     Estimated total dynamical RAM >      65.78 MB
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D
  The code is running with the 2D cutoff
  Please refer to:
  Sohier, T., Calandra, M., & Mauri, F. (2017), 
  Density functional perturbation theory for gated two-dimensional heterostructu
 res:
  Theoretical developments and application to flexural phonons in graphene.
  Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
 ----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D----2D

     Check: negative/imaginary core charge=   -0.000002    0.000000

     The potential is recalculated from file :
     /home/paulatto/espresso/tempdir/_ph0/bn.save/charge-density.dat


     negative rho (up, down):  1.316E-05 0.000E+00
     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 17.9

     total cpu time spent up to now is       90.2 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     highest occupied level (ev):    -5.8138

     Writing output data file bn.save

     Phonon dispersions for BN                                                  

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.7419  a.u.
     unit-cell volume          =     397.0601 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      55.0000  Ry
     charge density cut-off    =     440.0000  Ry
     convergence threshold     =      1.0E-15
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    4.74190  celldm(2)=    0.00000  celldm(3)=    4.30000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  4.3000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.2326 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     B   10.8100   tau(    1) = (    0.00000    0.28868    0.00000  )
        2     N   14.0067   tau(    2) = (    0.00000   -0.28868    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.2500000   0.5773503   0.0000000 )
 
      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  250.6096  (  15455 G-vectors)     FFT grid: ( 32, 32,144)
     G cutoff =  125.3048  (   5472 G-vectors)  smooth grid: ( 24, 24,100)
     number of k points=   288

     PseudoPot. # 1 for  B read from file:
     /home/paulatto/espresso/pseudo/B.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 92fadd7fada9510bd105ea601c5e9511
     Pseudo is Projector augmented-wave + core cor, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: BESSEL
     Using radial grid of 1059 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     PseudoPot. # 2 for  N read from file:
     /home/paulatto/espresso/pseudo/N.pbe-n-kjpaw_psl.0.1.UPF
     MD5 check sum: 51b9062c40abfbed4b0465f1412b4c8b
     Pseudo is Projector augmented-wave + core cor, Zval =  5.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.2 svn rev. 9415
     Shape of augmentation charge: PSQ
     Using radial grid of 1085 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A'  To be done

     Representation     4      1 modes -A'  To be done

     Representation     5      1 modes -A''  To be done

     Representation     6      1 modes -A''  To be done



     Alpha used in Ewald sum =   2.8000

     negative rho (up, down):  1.316E-05 0.000E+00
     PHONON       : 23m32.31s CPU    24m27.83s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :  1470.2 secs   av.it.:   7.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.354E-06

      iter #   2 total cpu time :  1473.4 secs   av.it.:  11.9
      thresh= 1.534E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.535E-07

      iter #   3 total cpu time :  1476.5 secs   av.it.:  11.5
      thresh= 8.084E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.273E-08

      iter #   4 total cpu time :  1479.7 secs   av.it.:  11.9
      thresh= 1.128E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.358E-09

      iter #   5 total cpu time :  1482.8 secs   av.it.:  11.7
      thresh= 3.685E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.120E-10

      iter #   6 total cpu time :  1485.9 secs   av.it.:  11.6
      thresh= 1.058E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.309E-12

      iter #   7 total cpu time :  1489.0 secs   av.it.:  11.3
      thresh= 1.144E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.202E-14

      iter #   8 total cpu time :  1492.1 secs   av.it.:  11.4
      thresh= 1.484E-08 alpha_mix =  0.700 |ddv_scf|^2 =  8.432E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :  1494.4 secs   av.it.:   7.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.267E-06

      iter #   2 total cpu time :  1497.4 secs   av.it.:  11.2
      thresh= 2.066E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.124E-06

      iter #   3 total cpu time :  1500.6 secs   av.it.:  11.7
      thresh= 1.060E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.287E-07

      iter #   4 total cpu time :  1503.7 secs   av.it.:  11.3
      thresh= 7.929E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.883E-09

      iter #   5 total cpu time :  1506.8 secs   av.it.:  11.8
      thresh= 5.369E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.209E-11

      iter #   6 total cpu time :  1510.0 secs   av.it.:  11.8
      thresh= 8.491E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.412E-12

      iter #   7 total cpu time :  1513.0 secs   av.it.:  11.2
      thresh= 1.188E-07 alpha_mix =  0.700 |ddv_scf|^2 =  8.461E-14

      iter #   8 total cpu time :  1516.0 secs   av.it.:  11.1
      thresh= 2.909E-08 alpha_mix =  0.700 |ddv_scf|^2 =  8.673E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :  1518.3 secs   av.it.:   7.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.381E-06

      iter #   2 total cpu time :  1521.3 secs   av.it.:  11.1
      thresh= 1.543E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.858E-07

      iter #   3 total cpu time :  1524.4 secs   av.it.:  11.9
      thresh= 6.211E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.356E-07

      iter #   4 total cpu time :  1527.5 secs   av.it.:  11.3
      thresh= 5.793E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.743E-09

      iter #   5 total cpu time :  1530.6 secs   av.it.:  11.8
      thresh= 5.237E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.243E-11

      iter #   6 total cpu time :  1533.8 secs   av.it.:  11.8
      thresh= 7.241E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.391E-12

      iter #   7 total cpu time :  1536.9 secs   av.it.:  11.4
      thresh= 1.179E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.489E-14

      iter #   8 total cpu time :  1539.9 secs   av.it.:  10.9
      thresh= 2.343E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.205E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :  1542.3 secs   av.it.:   7.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.938E-06

      iter #   2 total cpu time :  1545.5 secs   av.it.:  11.9
      thresh= 1.985E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.341E-06

      iter #   3 total cpu time :  1548.5 secs   av.it.:  11.5
      thresh= 1.158E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.454E-09

      iter #   4 total cpu time :  1551.7 secs   av.it.:  11.9
      thresh= 8.034E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.615E-10

      iter #   5 total cpu time :  1554.9 secs   av.it.:  12.0
      thresh= 2.148E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.828E-10

      iter #   6 total cpu time :  1558.1 secs   av.it.:  11.5
      thresh= 1.682E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.413E-13

      iter #   7 total cpu time :  1561.2 secs   av.it.:  11.4
      thresh= 8.610E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.215E-14

      iter #   8 total cpu time :  1564.4 secs   av.it.:  11.6
      thresh= 1.488E-08 alpha_mix =  0.700 |ddv_scf|^2 =  8.805E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :  1566.6 secs   av.it.:   6.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.203E-06

      iter #   2 total cpu time :  1569.1 secs   av.it.:   8.5
      thresh= 2.490E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.572E-07

      iter #   3 total cpu time :  1571.5 secs   av.it.:   8.3
      thresh= 7.465E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.636E-11

      iter #   4 total cpu time :  1574.5 secs   av.it.:  11.3
      thresh= 4.045E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.317E-13

      iter #   5 total cpu time :  1577.5 secs   av.it.:  11.0
      thresh= 7.292E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.761E-14

      iter #   6 total cpu time :  1580.4 secs   av.it.:  10.5
      thresh= 1.662E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.508E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :  1582.8 secs   av.it.:   7.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  7.574E-07

      iter #   2 total cpu time :  1585.9 secs   av.it.:  11.3
      thresh= 8.703E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.009E-08

      iter #   3 total cpu time :  1588.9 secs   av.it.:  11.2
      thresh= 1.004E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.722E-10

      iter #   4 total cpu time :  1591.7 secs   av.it.:  10.1
      thresh= 2.173E-06 alpha_mix =  0.700 |ddv_scf|^2 =  8.785E-13

      iter #   5 total cpu time :  1594.7 secs   av.it.:  10.9
      thresh= 9.373E-08 alpha_mix =  0.700 |ddv_scf|^2 =  9.963E-15

      iter #   6 total cpu time :  1597.6 secs   av.it.:  11.0
      thresh= 9.982E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.163E-16

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    6
     List of q in the star:
          1   0.250000000   0.577350269   0.000000000
          2  -0.250000000   0.577350269   0.000000000
          3  -0.375000000  -0.505181486   0.000000000
          4   0.625000000  -0.072168784   0.000000000
          5  -0.625000000  -0.072168784   0.000000000
          6   0.375000000  -0.505181486   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.250000000   0.577350269   0.000000000 ) 

 **************************************************************************
     freq (    1) =       9.120827 [THz] =     304.238039 [cm-1]
     freq (    2) =      17.958197 [THz] =     599.020971 [cm-1]
     freq (    3) =      24.230218 [THz] =     808.233062 [cm-1]
     freq (    4) =      32.076990 [THz] =    1069.973214 [cm-1]
     freq (    5) =      35.834938 [THz] =    1195.324865 [cm-1]
     freq (    6) =      37.677932 [THz] =    1256.800515 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =        304.2  [cm-1]   --> A''                
     freq (  2 -  2) =        599.0  [cm-1]   --> A''                
     freq (  3 -  3) =        808.2  [cm-1]   --> A'                 
     freq (  4 -  4) =       1070.0  [cm-1]   --> A'                 
     freq (  5 -  5) =       1195.3  [cm-1]   --> A'                 
     freq (  6 -  6) =       1256.8  [cm-1]   --> A'                 
 
     init_run     :      5.36s CPU      5.63s WALL (       9 calls)
     electrons    :     80.98s CPU     84.30s WALL (       9 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       9 calls)
     potinit      :      0.54s CPU      0.81s WALL (       9 calls)

     Called by electrons:
     c_bands      :     80.98s CPU     84.30s WALL (       9 calls)
     v_of_rho     :      0.23s CPU      0.23s WALL (      10 calls)
     newd         :      0.20s CPU      0.20s WALL (      10 calls)
     PAW_pot      :      0.24s CPU      0.24s WALL (      10 calls)

     Called by c_bands:
     init_us_2    :     11.17s CPU     11.88s WALL (   56254 calls)
     cegterg      :     71.53s CPU     74.42s WALL (    1994 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :    996.94s CPU   1037.68s WALL (  595003 calls)
     s_psi        :     47.74s CPU     50.04s WALL ( 1194921 calls)
     g_psi        :      0.56s CPU      0.58s WALL (   35634 calls)
     cdiaghg      :      2.85s CPU      2.98s WALL (   37628 calls)

     Called by h_psi:
     h_psi:pot    :    992.05s CPU   1032.52s WALL (  595003 calls)
     h_psi:calbec :     48.22s CPU     49.97s WALL (  595003 calls)
     vloc_psi     :    915.85s CPU    953.36s WALL (  595003 calls)
     add_vuspsi   :     24.08s CPU     25.22s WALL (  595003 calls)

     General routines
     calbec       :     92.25s CPU     96.12s WALL ( 1289558 calls)
     fft          :     11.75s CPU     11.92s WALL (    8535 calls)
     ffts         :      2.38s CPU      2.52s WALL (    7582 calls)
     fftw         :    939.34s CPU    979.87s WALL ( 4873008 calls)
     interpolate  :      1.16s CPU      1.32s WALL (     964 calls)
     davcio       :      5.83s CPU      6.48s WALL (  235279 calls)
 
     Parallel routines
     fft_scatt_xy :    101.48s CPU    105.77s WALL ( 4889125 calls)
     fft_scatt_yz :    371.82s CPU    386.77s WALL ( 4889125 calls)
 
     PHONON       : 25m34.82s CPU    26m37.73s WALL

     INITIALIZATION: 
     phq_setup    :      0.60s CPU      0.60s WALL (      10 calls)
     phq_init     :     28.73s CPU     29.18s WALL (      10 calls)
 
     phq_init     :     28.73s CPU     29.18s WALL (      10 calls)
     set_drhoc    :      9.97s CPU      9.98s WALL (      30 calls)
     init_vloc    :      0.49s CPU      0.50s WALL (      10 calls)
     init_us_1    :      1.37s CPU      1.37s WALL (      10 calls)
     newd         :      0.20s CPU      0.20s WALL (      10 calls)
     dvanqq       :      1.91s CPU      1.91s WALL (      10 calls)
     drho         :     10.35s CPU     10.70s WALL (      10 calls)
 
     DYNAMICAL MATRIX:
     dynmat0      :      8.92s CPU      8.97s WALL (      10 calls)
     phqscf       :   1396.68s CPU   1454.63s WALL (      10 calls)
     dynmatrix    :      0.02s CPU      0.07s WALL (      10 calls)
 
     phqscf       :   1396.68s CPU   1454.63s WALL (      10 calls)
     solve_linter :   1392.42s CPU   1450.09s WALL (      58 calls)
     drhodv       :      3.65s CPU      3.79s WALL (      58 calls)
 
     dynmat0      :      8.92s CPU      8.97s WALL (      10 calls)
     dynmat_us    :      0.90s CPU      0.93s WALL (      10 calls)
     d2ionq       :      0.06s CPU      0.06s WALL (      10 calls)
     dynmatcc     :      7.42s CPU      7.42s WALL (      10 calls)
 
     dynmat_us    :      0.90s CPU      0.93s WALL (      10 calls)
     addusdynmat  :      0.00s CPU      0.00s WALL (      10 calls)
 
     phqscf       :   1396.68s CPU   1454.63s WALL (      10 calls)
     solve_linter :   1392.42s CPU   1450.09s WALL (      58 calls)
 
     solve_linter :   1392.42s CPU   1450.09s WALL (      58 calls)
     dvqpsi_us    :     14.49s CPU     15.19s WALL (    6210 calls)
     ortho        :      9.30s CPU      9.75s WALL (   44480 calls)
     cgsolve      :   1104.30s CPU   1149.40s WALL (   44480 calls)
     incdrhoscf   :     86.89s CPU     91.97s WALL (   44480 calls)
     addusddens   :     12.25s CPU     12.31s WALL (     467 calls)
     vpsifft      :     64.83s CPU     68.59s WALL (   37928 calls)
     dv_of_drho   :      8.54s CPU      8.58s WALL (     444 calls)
     mix_pot      :      1.52s CPU      2.12s WALL (     416 calls)
     psymdvscf    :     14.80s CPU     14.89s WALL (     409 calls)
     newdq        :     10.95s CPU     10.98s WALL (     416 calls)
     adddvscf     :      2.01s CPU      2.12s WALL (   38270 calls)
     drhodvus     :      0.04s CPU      0.09s WALL (      58 calls)
 
     dvqpsi_us    :     14.49s CPU     15.19s WALL (    6210 calls)
     dvqpsi_us_on :      1.10s CPU      1.14s WALL (    6210 calls)
 
     cgsolve      :   1104.30s CPU   1149.40s WALL (   44480 calls)
     ch_psi       :   1068.22s CPU   1111.66s WALL (  555381 calls)
 
     ch_psi       :   1068.22s CPU   1111.66s WALL (  555381 calls)
     h_psi        :    996.94s CPU   1037.68s WALL (  595003 calls)
     last         :    103.00s CPU    107.19s WALL (  555381 calls)
 
     h_psi        :    996.94s CPU   1037.68s WALL (  595003 calls)
     add_vuspsi   :     24.08s CPU     25.22s WALL (  595003 calls)
 
     incdrhoscf   :     86.89s CPU     91.97s WALL (   44480 calls)
     addusdbec    :      5.74s CPU      6.01s WALL (   50576 calls)
 
     drhodvus     :      0.04s CPU      0.09s WALL (      58 calls)
 
      General routines
     calbec       :     92.25s CPU     96.12s WALL ( 1289558 calls)
     fft          :     11.75s CPU     11.92s WALL (    8535 calls)
     ffts         :      2.38s CPU      2.52s WALL (    7582 calls)
     fftw         :    939.34s CPU    979.87s WALL ( 4873008 calls)
     davcio       :      5.83s CPU      6.48s WALL (  235279 calls)
     write_rec    :      0.71s CPU      1.50s WALL (     474 calls)
 
 
     PHONON       : 25m34.82s CPU    26m37.73s WALL

 
   This run was terminated on:  21:25:33  25Oct2017            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_11 = """
     Program PHONON v.6.0 (svn rev. 13286) starts on  7Feb2017 at 14:35:56 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     path-images division:  nimage    =       2
     R & G space division:  proc/nbgrp/npool/nimage =       2

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/alas.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     229
     Max         121     121     43                 1224     1224     230
     Sum         241     241     85                 2445     2445     459



     Dynamical matrices for ( 4, 4, 4)  uniform grid of q-points
     (   8q-points):
       N         xq(1)         xq(2)         xq(3) 
       1   0.000000000   0.000000000   0.000000000
       2  -0.250000000   0.250000000  -0.250000000
       3   0.500000000  -0.500000000   0.500000000
       4   0.000000000   0.500000000   0.000000000
       5   0.750000000  -0.250000000   0.750000000
       6   0.500000000   0.000000000   0.500000000
       7   0.000000000  -1.000000000   0.000000000
       8  -0.500000000  -1.000000000   0.000000000

      Image parallelization. There are  2 images and    38 representations
      The estimated total work is   336 self-consistent (scf) runs
      I am image number     0 and my work is about  165 scf runs. I calculate: 
      q point number     1, representations:
       0 1 2
      q point number     2, representations:
       0 1 2 3 4
      q point number     3, representations:
       0 1 2 3 4
      q point number     4, representations:
       0 1 2 3 4 5 6
      q point number     5, representations:
       0 1 2 3 4

     Calculation of q =    0.0000000   0.0000000   0.0000000

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   0.0000000 )

     25 Sym.Ops. (with q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=     2

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, T_d (-43m)  point group:


     Electric field:
     Dielectric constant
     Born effective charges in two ways 


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      3 modes -T_2  G_15 P_4  To be done

     Representation     2      3 modes -T_2  G_15 P_4  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     0.25s CPU         0.26s WALL


     Electric Fields Calculation

      iter #   1 total cpu time :     0.3 secs   av.it.:   6.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.326E-06

      iter #   2 total cpu time :     0.4 secs   av.it.:   9.3
      thresh= 1.151E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.508E-08

      iter #   3 total cpu time :     0.4 secs   av.it.:   9.5
      thresh= 2.551E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.401E-10

      iter #   4 total cpu time :     0.5 secs   av.it.:   9.8
      thresh= 2.530E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.108E-12

      iter #   5 total cpu time :     0.5 secs   av.it.:   9.0
      thresh= 1.763E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.543E-14

     End of electric fields calculation

          Dielectric constant in cartesian axis 

          (      13.744216098      -0.000000000      -0.000000000 )
          (       0.000000000      13.744216098       0.000000000 )
          (      -0.000000000       0.000000000      13.744216098 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   Al 
      Ex  (        1.88265       -0.00000       -0.00000 )
      Ey  (       -0.00000        1.88265       -0.00000 )
      Ez  (        0.00000        0.00000        1.88265 )
           atom      2   As 
      Ex  (       -3.23374       -0.00000       -0.00000 )
      Ey  (        0.00000       -3.23374       -0.00000 )
      Ez  (       -0.00000       -0.00000       -3.23374 )


     Representation #  1 modes #   1  2  3

     Self-consistent Calculation

      iter #   1 total cpu time :     0.6 secs   av.it.:   5.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.662E-07

      iter #   2 total cpu time :     0.6 secs   av.it.:   9.7
      thresh= 6.828E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.273E-08

      iter #   3 total cpu time :     0.7 secs   av.it.:   9.7
      thresh= 1.508E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.797E-11

      iter #   4 total cpu time :     0.8 secs   av.it.:   9.5
      thresh= 6.162E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.182E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   4  5  6

     Self-consistent Calculation

      iter #   1 total cpu time :     0.8 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.910E-08

      iter #   2 total cpu time :     0.9 secs   av.it.:   9.8
      thresh= 1.706E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.259E-10

      iter #   3 total cpu time :     0.9 secs   av.it.:   9.5
      thresh= 1.805E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.012E-11

      iter #   4 total cpu time :     1.0 secs   av.it.:   9.5
      thresh= 5.488E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.306E-12

      iter #   5 total cpu time :     1.0 secs   av.it.:   9.5
      thresh= 1.143E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.628E-16

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

          Dielectric constant in cartesian axis 

          (      13.744216098      -0.000000000      -0.000000000 )
          (       0.000000000      13.744216098       0.000000000 )
          (      -0.000000000       0.000000000      13.744216098 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   Al 
      Ex  (        1.88265       -0.00000       -0.00000 )
      Ey  (       -0.00000        1.88265       -0.00000 )
      Ez  (        0.00000        0.00000        1.88265 )
           atom      2   As 
      Ex  (       -3.23374       -0.00000       -0.00000 )
      Ey  (        0.00000       -3.23374       -0.00000 )
      Ez  (       -0.00000       -0.00000       -3.23374 )

          Effective charges (d P / du) in cartesian axis 

           atom      1   Al 
      Px  (        1.88284       -0.00000       -0.00000 )
      Py  (       -0.00000        1.88284        0.00000 )
      Pz  (       -0.00000        0.00000        1.88284 )
           atom      2   As 
      Px  (       -3.23837       -0.00000        0.00000 )
      Py  (       -0.00000       -3.23837       -0.00000 )
      Pz  (        0.00000       -0.00000       -3.23837 )

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    2) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    3) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    4) =      11.258806 [THz] =     375.553348 [cm-1]
     freq (    5) =      11.258806 [THz] =     375.553348 [cm-1]
     freq (    6) =      11.258806 [THz] =     375.553348 [cm-1]
 **************************************************************************

     Mode symmetry, T_d (-43m)  point group:

     freq (  1 -  3) =          5.5  [cm-1]   --> T_2  G_15 P_4   I+R
     freq (  4 -  6) =        375.6  [cm-1]   --> T_2  G_15 P_4   I+R

     Calculation of q =   -0.2500000   0.2500000  -0.2500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     264
     Max         121     121     43                 1224     1224     267
     Sum         241     241     85                 2445     2445     531



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    20
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.1875000
        k(    2) = (   0.0000000   0.5000000   0.0000000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.3750000
        k(    4) = (   0.0000000   0.5000000   0.5000000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0625000
        k(    6) = (  -0.5000000   0.5000000  -0.5000000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.1875000
        k(    8) = (  -0.5000000   0.0000000  -0.5000000), wk =   0.0000000
        k(    9) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0625000
        k(   10) = (   0.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   11) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.1875000
        k(   12) = (  -1.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   13) = (  -0.7500000   0.2500000  -0.2500000), wk =   0.1875000
        k(   14) = (  -1.0000000   0.5000000  -0.5000000), wk =   0.0000000
        k(   15) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.3750000
        k(   16) = (  -0.5000000   0.0000000  -1.0000000), wk =   0.0000000
        k(   17) = (   0.2500000  -0.2500000   0.7500000), wk =   0.1875000
        k(   18) = (   0.0000000   0.0000000   0.5000000), wk =   0.0000000
        k(   19) = (  -0.2500000   0.2500000   0.7500000), wk =   0.1875000
        k(   20) = (  -0.5000000   0.5000000   0.5000000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.0

     total cpu time spent up to now is        0.2 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.0000 0.5000 0.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.0000 0.5000 0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.5000 0.5000-0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.5000 0.0000-0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k = 0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.0000 0.0000 0.0000 (   331 PWs)   bands (ev):

    -6.9797   5.1761   5.1761   5.1761

          k =-0.7500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-1.0000 0.0000 0.0000 (   302 PWs)   bands (ev):

    -4.8217  -0.4470   2.9274   2.9274

          k =-0.7500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-1.0000 0.5000-0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.5000 0.0000-1.0000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k = 0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.0000 0.0000 0.5000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.5000 0.5000 0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (  -0.2500000   0.2500000  -0.2500000 )

      6 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    20

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      1 modes -A_1  L_1  To be done

     Representation     3      2 modes -E    L_3  To be done

     Representation     4      2 modes -E    L_3  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     1.26s CPU         1.32s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     1.4 secs   av.it.:   6.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.084E-03

      iter #   2 total cpu time :     1.4 secs   av.it.:   7.6
      thresh= 5.554E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.965E-02

      iter #   3 total cpu time :     1.4 secs   av.it.:   6.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.334E-06

      iter #   4 total cpu time :     1.5 secs   av.it.:   7.2
      thresh= 2.517E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.380E-07

      iter #   5 total cpu time :     1.5 secs   av.it.:   7.6
      thresh= 3.714E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.454E-09

      iter #   6 total cpu time :     1.6 secs   av.it.:   7.0
      thresh= 6.674E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.717E-10

      iter #   7 total cpu time :     1.6 secs   av.it.:   7.2
      thresh= 2.172E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.731E-11

      iter #   8 total cpu time :     1.7 secs   av.it.:   7.2
      thresh= 6.108E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.132E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     1.7 secs   av.it.:   5.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.514E-04

      iter #   2 total cpu time :     1.8 secs   av.it.:   7.6
      thresh= 2.552E-03 alpha_mix =  0.700 |ddv_scf|^2 =  5.928E-03

      iter #   3 total cpu time :     1.8 secs   av.it.:   6.2
      thresh= 7.699E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.208E-07

      iter #   4 total cpu time :     1.9 secs   av.it.:   8.2
      thresh= 4.699E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.970E-09

      iter #   5 total cpu time :     1.9 secs   av.it.:   8.1
      thresh= 8.348E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.212E-10

      iter #   6 total cpu time :     1.9 secs   av.it.:   7.4
      thresh= 2.283E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.038E-09

      iter #   7 total cpu time :     2.0 secs   av.it.:   6.9
      thresh= 3.223E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.551E-11

      iter #   8 total cpu time :     2.0 secs   av.it.:   7.6
      thresh= 3.938E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.433E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 modes #   3  4

     Self-consistent Calculation

      iter #   1 total cpu time :     2.1 secs   av.it.:   5.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.311E-06

      iter #   2 total cpu time :     2.2 secs   av.it.:   9.2
      thresh= 1.145E-04 alpha_mix =  0.700 |ddv_scf|^2 =  9.088E-08

      iter #   3 total cpu time :     2.3 secs   av.it.:   9.2
      thresh= 3.015E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.906E-11

      iter #   4 total cpu time :     2.4 secs   av.it.:   9.2
      thresh= 9.437E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.651E-12

      iter #   5 total cpu time :     2.5 secs   av.it.:   9.0
      thresh= 1.285E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.874E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 modes #   5  6

     Self-consistent Calculation

      iter #   1 total cpu time :     2.6 secs   av.it.:   5.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.108E-07

      iter #   2 total cpu time :     2.7 secs   av.it.:   9.4
      thresh= 3.328E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.511E-09

      iter #   3 total cpu time :     2.8 secs   av.it.:   9.2
      thresh= 6.717E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.323E-10

      iter #   4 total cpu time :     2.9 secs   av.it.:   9.1
      thresh= 1.150E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.943E-12

      iter #   5 total cpu time :     3.0 secs   av.it.:   8.8
      thresh= 2.635E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.124E-15

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    4
     List of q in the star:
          1  -0.250000000   0.250000000  -0.250000000
          2   0.250000000  -0.250000000  -0.250000000
          3  -0.250000000  -0.250000000   0.250000000
          4   0.250000000   0.250000000   0.250000000
     In addition there is the -q list: 
          1   0.250000000  -0.250000000   0.250000000
          2  -0.250000000   0.250000000   0.250000000
          3   0.250000000   0.250000000  -0.250000000
          4  -0.250000000  -0.250000000  -0.250000000

     Diagonalizing the dynamical matrix

     q = (   -0.250000000   0.250000000  -0.250000000 ) 

 **************************************************************************
     freq (    1) =       1.761214 [THz] =      58.747782 [cm-1]
     freq (    2) =       1.761214 [THz] =      58.747782 [cm-1]
     freq (    3) =       4.534095 [THz] =     151.241127 [cm-1]
     freq (    4) =      11.004844 [THz] =     367.082097 [cm-1]
     freq (    5) =      11.004844 [THz] =     367.082097 [cm-1]
     freq (    6) =      12.136604 [THz] =     404.833529 [cm-1]
 **************************************************************************

     Mode symmetry, C_3v (3m)   point group:

     freq (  1 -  2) =         58.7  [cm-1]   --> E    L_3           
     freq (  3 -  3) =        151.2  [cm-1]   --> A_1  L_1           
     freq (  4 -  5) =        367.1  [cm-1]   --> E    L_3           
     freq (  6 -  6) =        404.8  [cm-1]   --> A_1  L_1           

     Calculation of q =    0.5000000  -0.5000000   0.5000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     267
     Max         121     121     43                 1224     1224     270
     Sum         241     241     85                 2445     2445     537



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    10
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.3750000
        k(    2) = (   0.7500000  -0.2500000   0.7500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.7500000
        k(    4) = (   0.7500000  -0.2500000   1.2500000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(    6) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0000000
        k(    7) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.3750000
        k(    8) = (  -0.2500000  -0.7500000   0.7500000), wk =   0.0000000
        k(    9) = (  -0.7500000   0.2500000  -0.2500000), wk =   0.3750000
        k(   10) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.6

     total cpu time spent up to now is        0.3 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.7500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.7500-0.2500 1.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.7500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.7500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.7500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.5000000  -0.5000000   0.5000000 )

      7 Sym.Ops. (with q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    10

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      1 modes -A_1  L_1  To be done

     Representation     3      2 modes -E    L_3  To be done

     Representation     4      2 modes -E    L_3  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     3.03s CPU         3.22s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     3.2 secs   av.it.:   6.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.569E-04

      iter #   2 total cpu time :     3.3 secs   av.it.:   8.2
      thresh= 1.889E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.022E-03

      iter #   3 total cpu time :     3.3 secs   av.it.:   7.4
      thresh= 3.197E-03 alpha_mix =  0.700 |ddv_scf|^2 =  5.259E-08

      iter #   4 total cpu time :     3.3 secs   av.it.:   8.0
      thresh= 2.293E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.120E-09

      iter #   5 total cpu time :     3.3 secs   av.it.:   7.4
      thresh= 9.011E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.293E-11

      iter #   6 total cpu time :     3.4 secs   av.it.:   8.4
      thresh= 6.552E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.553E-12

      iter #   7 total cpu time :     3.4 secs   av.it.:   8.0
      thresh= 2.134E-07 alpha_mix =  0.700 |ddv_scf|^2 =  8.070E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     3.4 secs   av.it.:   5.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.799E-05

      iter #   2 total cpu time :     3.5 secs   av.it.:   8.2
      thresh= 7.615E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.451E-04

      iter #   3 total cpu time :     3.5 secs   av.it.:   7.4
      thresh= 1.205E-03 alpha_mix =  0.700 |ddv_scf|^2 =  6.734E-07

      iter #   4 total cpu time :     3.5 secs   av.it.:   7.6
      thresh= 8.206E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.035E-09

      iter #   5 total cpu time :     3.6 secs   av.it.:   8.0
      thresh= 6.352E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.159E-11

      iter #   6 total cpu time :     3.6 secs   av.it.:   8.4
      thresh= 8.461E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.226E-12

      iter #   7 total cpu time :     3.7 secs   av.it.:   8.2
      thresh= 1.107E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.358E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 modes #   3  4

     Self-consistent Calculation

      iter #   1 total cpu time :     3.8 secs   av.it.:   6.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.539E-06

      iter #   2 total cpu time :     3.8 secs   av.it.:   9.2
      thresh= 1.241E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.304E-07

      iter #   3 total cpu time :     3.9 secs   av.it.:   9.0
      thresh= 3.611E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.951E-11

      iter #   4 total cpu time :     3.9 secs   av.it.:   9.2
      thresh= 9.461E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.026E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 modes #   5  6

     Self-consistent Calculation

      iter #   1 total cpu time :     4.0 secs   av.it.:   4.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.480E-07

      iter #   2 total cpu time :     4.0 secs   av.it.:   9.0
      thresh= 3.847E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.827E-09

      iter #   3 total cpu time :     4.1 secs   av.it.:   9.0
      thresh= 9.395E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.470E-10

      iter #   4 total cpu time :     4.1 secs   av.it.:   9.1
      thresh= 1.212E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.522E-12

      iter #   5 total cpu time :     4.2 secs   av.it.:   8.3
      thresh= 2.743E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.741E-15

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    4
     List of q in the star:
          1   0.500000000  -0.500000000   0.500000000
          2  -0.500000000   0.500000000   0.500000000
          3   0.500000000   0.500000000  -0.500000000
          4  -0.500000000  -0.500000000  -0.500000000

     Diagonalizing the dynamical matrix

     q = (    0.500000000  -0.500000000   0.500000000 ) 

 **************************************************************************
     freq (    1) =       2.016390 [THz] =      67.259545 [cm-1]
     freq (    2) =       2.016390 [THz] =      67.259545 [cm-1]
     freq (    3) =       6.494357 [THz] =     216.628437 [cm-1]
     freq (    4) =      10.940872 [THz] =     364.948217 [cm-1]
     freq (    5) =      10.940872 [THz] =     364.948217 [cm-1]
     freq (    6) =      11.551694 [THz] =     385.323024 [cm-1]
 **************************************************************************

     Mode symmetry, C_3v (3m)   point group:

     freq (  1 -  2) =         67.3  [cm-1]   --> E    L_3           
     freq (  3 -  3) =        216.6  [cm-1]   --> A_1  L_1           
     freq (  4 -  5) =        364.9  [cm-1]   --> E    L_3           
     freq (  6 -  6) =        385.3  [cm-1]   --> A_1  L_1           

     Calculation of q =    0.0000000   0.5000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     264
     Max         121     121     43                 1224     1224     267
     Sum         241     241     85                 2445     2445     531



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    24
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.1250000
        k(    2) = (   0.2500000   0.7500000   0.2500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.2500000
        k(    4) = (   0.2500000   0.7500000   0.7500000), wk =   0.0000000
        k(    5) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.1250000
        k(    6) = (  -0.2500000   0.2500000   0.2500000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.1250000
        k(    8) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0000000
        k(    9) = (   0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(   10) = (   0.2500000   0.7500000  -0.2500000), wk =   0.0000000
        k(   11) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.2500000
        k(   12) = (  -0.2500000   0.2500000   0.7500000), wk =   0.0000000
        k(   13) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.1250000
        k(   14) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.0000000
        k(   15) = (  -0.2500000   0.7500000  -0.2500000), wk =   0.1250000
        k(   16) = (  -0.2500000   1.2500000  -0.2500000), wk =   0.0000000
        k(   17) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.2500000
        k(   18) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.0000000
        k(   19) = (   0.2500000   0.2500000  -0.7500000), wk =   0.2500000
        k(   20) = (   0.2500000   0.7500000  -0.7500000), wk =   0.0000000
        k(   21) = (  -0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   22) = (  -0.2500000   1.2500000   0.2500000), wk =   0.0000000
        k(   23) = (  -0.2500000  -0.7500000  -0.2500000), wk =   0.1250000
        k(   24) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.7

     total cpu time spent up to now is        0.6 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.7500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500 0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500 0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 1.2500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.7500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 1.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.5000000   0.0000000 )

      4 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    24

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  To be done

     Representation     2      1 modes -A_1  D_1  S_1  To be done

     Representation     3      1 modes -B_1  D_3  S_3  To be done

     Representation     4      1 modes -B_1  D_3  S_3  To be done

     Representation     5      1 modes -B_2  D_4  S_4  To be done

     Representation     6      1 modes -B_2  D_4  S_4  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     4.16s CPU         4.52s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     4.6 secs   av.it.:   6.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.919E-03

      iter #   2 total cpu time :     4.6 secs   av.it.:   8.0
      thresh= 4.381E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.597E-02

      iter #   3 total cpu time :     4.7 secs   av.it.:   7.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.108E-06

      iter #   4 total cpu time :     4.7 secs   av.it.:   8.3
      thresh= 1.452E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.578E-08

      iter #   5 total cpu time :     4.8 secs   av.it.:   8.7
      thresh= 1.606E-05 alpha_mix =  0.700 |ddv_scf|^2 =  7.158E-11

      iter #   6 total cpu time :     4.9 secs   av.it.:   8.2
      thresh= 8.460E-07 alpha_mix =  0.700 |ddv_scf|^2 =  9.928E-11

      iter #   7 total cpu time :     4.9 secs   av.it.:   7.1
      thresh= 9.964E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.226E-11

      iter #   8 total cpu time :     5.0 secs   av.it.:   7.2
      thresh= 5.679E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.114E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     5.0 secs   av.it.:   5.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.751E-04

      iter #   2 total cpu time :     5.1 secs   av.it.:   8.0
      thresh= 1.937E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.908E-03

      iter #   3 total cpu time :     5.1 secs   av.it.:   6.7
      thresh= 5.393E-03 alpha_mix =  0.700 |ddv_scf|^2 =  5.668E-07

      iter #   4 total cpu time :     5.2 secs   av.it.:   7.8
      thresh= 7.528E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.601E-09

      iter #   5 total cpu time :     5.3 secs   av.it.:   8.7
      thresh= 7.484E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.285E-11

      iter #   6 total cpu time :     5.3 secs   av.it.:   8.3
      thresh= 7.270E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.730E-12

      iter #   7 total cpu time :     5.4 secs   av.it.:   7.9
      thresh= 2.780E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.587E-11

      iter #   8 total cpu time :     5.4 secs   av.it.:   6.9
      thresh= 3.983E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.867E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :     5.5 secs   av.it.:   5.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.766E-06

      iter #   2 total cpu time :     5.5 secs   av.it.:   8.4
      thresh= 2.961E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.179E-06

      iter #   3 total cpu time :     5.6 secs   av.it.:   8.2
      thresh= 1.086E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.850E-10

      iter #   4 total cpu time :     5.6 secs   av.it.:   8.0
      thresh= 1.962E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.528E-11

      iter #   5 total cpu time :     5.7 secs   av.it.:   8.2
      thresh= 3.908E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.631E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :     5.8 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.095E-06

      iter #   2 total cpu time :     5.8 secs   av.it.:   8.4
      thresh= 1.046E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.263E-07

      iter #   3 total cpu time :     5.9 secs   av.it.:   8.3
      thresh= 3.554E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.599E-10

      iter #   4 total cpu time :     5.9 secs   av.it.:   7.9
      thresh= 2.569E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.110E-11

      iter #   5 total cpu time :     6.0 secs   av.it.:   7.9
      thresh= 4.594E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.831E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :     6.0 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.095E-06

      iter #   2 total cpu time :     6.1 secs   av.it.:   8.4
      thresh= 1.046E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.261E-07

      iter #   3 total cpu time :     6.2 secs   av.it.:   8.2
      thresh= 3.552E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.596E-10

      iter #   4 total cpu time :     6.3 secs   av.it.:   7.9
      thresh= 2.568E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.115E-11

      iter #   5 total cpu time :     6.4 secs   av.it.:   7.8
      thresh= 4.599E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.785E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :     6.4 secs   av.it.:   5.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.766E-06

      iter #   2 total cpu time :     6.5 secs   av.it.:   8.4
      thresh= 2.961E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.179E-06

      iter #   3 total cpu time :     6.5 secs   av.it.:   8.1
      thresh= 1.086E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.848E-10

      iter #   4 total cpu time :     6.6 secs   av.it.:   8.0
      thresh= 1.962E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.525E-11

      iter #   5 total cpu time :     6.6 secs   av.it.:   8.2
      thresh= 3.905E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.624E-14

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    6
     List of q in the star:
          1   0.000000000   0.500000000   0.000000000
          2  -0.500000000   0.000000000   0.000000000
          3   0.000000000  -0.500000000   0.000000000
          4   0.000000000   0.000000000   0.500000000
          5   0.000000000   0.000000000  -0.500000000
          6   0.500000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.500000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       2.421101 [THz] =      80.759230 [cm-1]
     freq (    2) =       2.421101 [THz] =      80.759230 [cm-1]
     freq (    3) =       4.606324 [THz] =     153.650423 [cm-1]
     freq (    4) =      10.666710 [THz] =     355.803149 [cm-1]
     freq (    5) =      10.666710 [THz] =     355.803149 [cm-1]
     freq (    6) =      12.371391 [THz] =     412.665187 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =         80.8  [cm-1]   --> B_1  D_3  S_3      
     freq (  2 -  2) =         80.8  [cm-1]   --> B_2  D_4  S_4      
     freq (  3 -  3) =        153.7  [cm-1]   --> A_1  D_1  S_1      
     freq (  4 -  4) =        355.8  [cm-1]   --> B_1  D_3  S_3      
     freq (  5 -  5) =        355.8  [cm-1]   --> B_2  D_4  S_4      
     freq (  6 -  6) =        412.7  [cm-1]   --> A_1  D_1  S_1      

     Calculation of q =    0.7500000  -0.2500000   0.7500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     48                 1221     1221     322
     Max         121     121     49                 1224     1224     323
     Sum         241     241     97                 2445     2445     645



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    40
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.0625000
        k(    2) = (   1.0000000   0.0000000   1.0000000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(    4) = (   1.0000000   0.0000000   1.5000000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0625000
        k(    6) = (   0.5000000   0.0000000   0.5000000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.1250000
        k(    8) = (   0.5000000  -0.5000000   1.0000000), wk =   0.0000000
        k(    9) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.0625000
        k(   10) = (   0.5000000  -0.5000000   0.5000000), wk =   0.0000000
        k(   11) = (   0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(   12) = (   1.0000000   0.0000000   0.5000000), wk =   0.0000000
        k(   13) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0625000
        k(   14) = (   1.0000000  -0.5000000   1.0000000), wk =   0.0000000
        k(   15) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   16) = (   0.5000000   0.0000000   0.0000000), wk =   0.0000000
        k(   17) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   18) = (   0.5000000  -0.5000000   1.5000000), wk =   0.0000000
        k(   19) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.1250000
        k(   20) = (   0.5000000  -1.0000000   1.0000000), wk =   0.0000000
        k(   21) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.1250000
        k(   22) = (   0.0000000  -0.5000000   1.0000000), wk =   0.0000000
        k(   23) = (  -0.2500000   0.7500000  -0.2500000), wk =   0.0625000
        k(   24) = (   0.5000000   0.5000000   0.5000000), wk =   0.0000000
        k(   25) = (   0.2500000   0.7500000   0.2500000), wk =   0.0625000
        k(   26) = (   1.0000000   0.5000000   1.0000000), wk =   0.0000000
        k(   27) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.1250000
        k(   28) = (   0.5000000  -0.5000000   0.0000000), wk =   0.0000000
        k(   29) = (   0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   30) = (   1.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   31) = (   0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   32) = (   1.0000000  -0.5000000   1.5000000), wk =   0.0000000
        k(   33) = (  -0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(   34) = (   0.5000000   0.0000000   1.5000000), wk =   0.0000000
        k(   35) = (  -0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   36) = (   0.5000000   0.5000000   1.0000000), wk =   0.0000000
        k(   37) = (  -0.2500000  -0.7500000  -0.2500000), wk =   0.0625000
        k(   38) = (   0.5000000  -1.0000000   0.5000000), wk =   0.0000000
        k(   39) = (   0.2500000  -0.7500000   0.2500000), wk =   0.0625000
        k(   40) = (   1.0000000  -1.0000000   1.0000000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.0

     total cpu time spent up to now is        1.0 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 1.0000 0.0000 1.0000 (   302 PWs)   bands (ev):

    -4.8217  -0.4470   2.9274   2.9274

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000 0.0000 1.5000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.5000 0.0000 0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.5000-0.5000 1.0000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.5000-0.5000 0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k = 0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 1.0000 0.0000 0.5000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k = 0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 1.0000-0.5000 1.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.0000 0.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-0.5000 1.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k =-0.2500-0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-1.0000 1.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.7500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.0000-0.5000 1.0000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k =-0.2500 0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.5000 0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k = 0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000 0.5000 1.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500-0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-0.5000 0.0000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k = 0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000 0.0000 0.0000 (   302 PWs)   bands (ev):

    -4.8217  -0.4470   2.9274   2.9274

          k = 0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000-0.5000 1.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.0000 1.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.5000 1.0000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-1.0000 0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k = 0.2500-0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000-1.0000 1.0000 (   331 PWs)   bands (ev):

    -6.9797   5.1761   5.1761   5.1761

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.7500000  -0.2500000   0.7500000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    40

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A'  To be done

     Representation     4      1 modes -A'  To be done

     Representation     5      1 modes -A''  Not done in this run

     Representation     6      1 modes -A''  Not done in this run

     Compute atoms:     1,    2,



     Alpha used in Ewald sum =   0.7000
     PHONON       :     6.38s CPU         7.11s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     7.2 secs   av.it.:   6.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.089E-04

      iter #   2 total cpu time :     7.3 secs   av.it.:   8.7
      thresh= 1.044E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.311E-04

      iter #   3 total cpu time :     7.4 secs   av.it.:   7.8
      thresh= 1.520E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.033E-06

      iter #   4 total cpu time :     7.5 secs   av.it.:   8.5
      thresh= 1.017E-04 alpha_mix =  0.700 |ddv_scf|^2 =  7.459E-09

      iter #   5 total cpu time :     7.6 secs   av.it.:   8.7
      thresh= 8.636E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.044E-10

      iter #   6 total cpu time :     7.7 secs   av.it.:   8.6
      thresh= 2.459E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.956E-12

      iter #   7 total cpu time :     7.8 secs   av.it.:   8.6
      thresh= 3.155E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.216E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     7.8 secs   av.it.:   5.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.249E-05

      iter #   2 total cpu time :     7.9 secs   av.it.:   8.8
      thresh= 5.700E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.416E-05

      iter #   3 total cpu time :     8.0 secs   av.it.:   7.8
      thresh= 8.010E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.653E-07

      iter #   4 total cpu time :     8.1 secs   av.it.:   8.2
      thresh= 5.151E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.893E-09

      iter #   5 total cpu time :     8.2 secs   av.it.:   8.5
      thresh= 6.240E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.760E-10

      iter #   6 total cpu time :     8.3 secs   av.it.:   8.7
      thresh= 1.661E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.660E-11

      iter #   7 total cpu time :     8.4 secs   av.it.:   8.7
      thresh= 4.074E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.743E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :     8.4 secs   av.it.:   6.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.537E-04

      iter #   2 total cpu time :     8.5 secs   av.it.:   8.7
      thresh= 1.240E-03 alpha_mix =  0.700 |ddv_scf|^2 =  3.325E-04

      iter #   3 total cpu time :     8.6 secs   av.it.:   7.8
      thresh= 1.824E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.125E-06

      iter #   4 total cpu time :     8.7 secs   av.it.:   8.4
      thresh= 1.061E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.383E-09

      iter #   5 total cpu time :     8.8 secs   av.it.:   8.8
      thresh= 7.990E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.601E-10

      iter #   6 total cpu time :     8.9 secs   av.it.:   8.5
      thresh= 2.367E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.256E-11

      iter #   7 total cpu time :     9.0 secs   av.it.:   8.4
      thresh= 3.545E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.821E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :     9.0 secs   av.it.:   5.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  9.536E-06

      iter #   2 total cpu time :     9.1 secs   av.it.:   8.9
      thresh= 3.088E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.109E-05

      iter #   3 total cpu time :     9.2 secs   av.it.:   8.2
      thresh= 3.330E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.916E-07

      iter #   4 total cpu time :     9.3 secs   av.it.:   8.2
      thresh= 6.258E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.967E-09

      iter #   5 total cpu time :     9.4 secs   av.it.:   8.6
      thresh= 5.447E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.110E-10

      iter #   6 total cpu time :     9.5 secs   av.it.:   8.7
      thresh= 1.453E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.069E-11

      iter #   7 total cpu time :     9.6 secs   av.it.:   8.7
      thresh= 3.269E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.710E-13

     End of self-consistent calculation

     Convergence has been achieved 

     Not diagonalizing because representation    5 is not done

     init_run     :      0.08s CPU      0.11s WALL (       4 calls)
     electrons    :      0.61s CPU      0.67s WALL (       4 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       4 calls)
     potinit      :      0.00s CPU      0.01s WALL (       4 calls)

     Called by electrons:
     c_bands      :      0.61s CPU      0.67s WALL (       4 calls)
     v_of_rho     :      0.00s CPU      0.00s WALL (       5 calls)

     Called by c_bands:
     init_us_2    :      0.08s CPU      0.09s WALL (    1809 calls)
     cegterg      :      0.53s CPU      0.57s WALL (      94 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :      4.66s CPU      5.15s WALL (   16197 calls)
     g_psi        :      0.01s CPU      0.01s WALL (    1056 calls)
     cdiaghg      :      0.06s CPU      0.07s WALL (    1150 calls)

     Called by h_psi:
     h_psi:pot    :      4.62s CPU      5.09s WALL (   16197 calls)
     h_psi:calbec :      0.35s CPU      0.34s WALL (   16197 calls)
     vloc_psi     :      3.95s CPU      4.46s WALL (   16197 calls)
     add_vuspsi   :      0.26s CPU      0.25s WALL (   16197 calls)

     General routines
     calbec       :      0.64s CPU      0.62s WALL (   32662 calls)
     fft          :      0.04s CPU      0.04s WALL (     542 calls)
     ffts         :      0.03s CPU      0.02s WALL (     386 calls)
     fftw         :      3.95s CPU      4.52s WALL (  141794 calls)
     davcio       :      0.08s CPU      0.09s WALL (    8711 calls)

     Parallel routines
     fft_scatter  :      1.07s CPU      1.33s WALL (  142722 calls)

     PHONON       :     8.72s CPU         9.58s WALL

     INITIALIZATION: 
     phq_setup    :      0.02s CPU      0.02s WALL (       5 calls)
     phq_init     :      0.16s CPU      0.16s WALL (       5 calls)

     phq_init     :      0.16s CPU      0.16s WALL (       5 calls)
     init_vloc    :      0.00s CPU      0.01s WALL (       5 calls)
     init_us_1    :      0.04s CPU      0.05s WALL (       5 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.06s CPU      0.06s WALL (       5 calls)
     phqscf       :      7.02s CPU      7.76s WALL (       5 calls)
     dynmatrix    :      0.01s CPU      0.02s WALL (       5 calls)

     phqscf       :      7.02s CPU      7.76s WALL (       5 calls)
     solve_linter :      6.94s CPU      7.65s WALL (      20 calls)
     drhodv       :      0.05s CPU      0.05s WALL (      20 calls)

     dynmat0      :      0.06s CPU      0.06s WALL (       5 calls)
     dynmat_us    :      0.04s CPU      0.04s WALL (       5 calls)
     d2ionq       :      0.01s CPU      0.01s WALL (       5 calls)

     dynmat_us    :      0.04s CPU      0.04s WALL (       5 calls)

     phqscf       :      7.02s CPU      7.76s WALL (       5 calls)
     solve_linter :      6.94s CPU      7.65s WALL (      20 calls)

     solve_linter :      6.94s CPU      7.65s WALL (      20 calls)
     dvqpsi_us    :      0.15s CPU      0.12s WALL (     266 calls)
     ortho        :      0.03s CPU      0.04s WALL (    1602 calls)
     cgsolve      :      5.28s CPU      5.78s WALL (    1602 calls)
     incdrhoscf   :      0.37s CPU      0.49s WALL (    1596 calls)
     vpsifft      :      0.39s CPU      0.42s WALL (    1312 calls)
     dv_of_drho   :      0.05s CPU      0.04s WALL (     174 calls)
     mix_pot      :      0.06s CPU      0.06s WALL (     127 calls)
     psymdvscf    :      0.53s CPU      0.52s WALL (     122 calls)

     dvqpsi_us    :      0.15s CPU      0.12s WALL (     266 calls)
     dvqpsi_us_on :      0.02s CPU      0.02s WALL (     266 calls)

     cgsolve      :      5.28s CPU      5.78s WALL (    1602 calls)
     ch_psi       :      4.96s CPU      5.47s WALL (   14953 calls)

     ch_psi       :      4.96s CPU      5.47s WALL (   14953 calls)
     h_psi        :      4.66s CPU      5.15s WALL (   16197 calls)
     last         :      0.65s CPU      0.64s WALL (   14953 calls)

     h_psi        :      4.66s CPU      5.15s WALL (   16197 calls)
     add_vuspsi   :      0.26s CPU      0.25s WALL (   16197 calls)

     incdrhoscf   :      0.37s CPU      0.49s WALL (    1596 calls)


      General routines
     calbec       :      0.64s CPU      0.62s WALL (   32662 calls)
     fft          :      0.04s CPU      0.04s WALL (     542 calls)
     ffts         :      0.03s CPU      0.02s WALL (     386 calls)
     fftw         :      3.95s CPU      4.52s WALL (  141794 calls)
     davcio       :      0.08s CPU      0.09s WALL (    8711 calls)
     write_rec    :      0.18s CPU      0.22s WALL (     147 calls)


     PHONON       :     8.72s CPU         9.58s WALL


   This run was terminated on:  14:36: 6   7Feb2017            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_16 = """
     Program PHONON v.6.1 starts on 16Mar2017 at 21: 4:50 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     K-points division:     npool     =       4

     Reading data from directory:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/aluminum.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    109                 3479     3479     749


     Check: negative/imaginary core charge=   -0.000013    0.000000


     Dynamical matrices for ( 4, 4, 4)  uniform grid of q-points
     With a half shift
     (  10q-points):
       N         xq(1)         xq(2)         xq(3) 
       1  -0.125000000   0.125000000   0.125000000
       2  -0.375000000   0.375000000  -0.125000000
       3   0.375000000  -0.375000000   0.625000000
       4   0.125000000  -0.125000000   0.375000000
       5  -0.125000000   0.625000000   0.125000000
       6   0.625000000  -0.125000000   0.875000000
       7   0.375000000   0.125000000   0.625000000
       8  -0.125000000  -0.875000000   0.125000000
       9  -0.375000000   0.375000000   0.375000000
      10   0.375000000  -0.375000000   1.125000000
     Because shifted q grid is used, q2r will not work !

     Calculation of q =   -0.1250000   0.1250000   0.1250000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    121                 3479     3479     869



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  1632 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_1/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.4

     total cpu time spent up to now is       75.8 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.1250000   0.1250000   0.1250000 )

      6 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  1632

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      2 modes -E    L_3  To be done


     PHONON       :  1m10.19s CPU     1m28.44s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn1

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 11.768747
     lambda(    1)=  0.2178   gamma=    2.25 GHz
     lambda(    2)=  0.2178   gamma=    2.25 GHz
     lambda(    3)=  0.1554   gamma=    5.56 GHz

     Calculation of q =   -0.3750000   0.3750000  -0.1250000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    121                 3479     3479     941



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  4352 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_2/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.1

     total cpu time spent up to now is      273.2 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.3750000   0.3750000  -0.1250000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  4352

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A''  To be done


     PHONON       :  4m 8.58s CPU     5m13.46s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn2

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 15.327453
     lambda(    1)=  0.0224   gamma=    0.94 GHz
     lambda(    2)=  0.1269   gamma=    7.46 GHz
     lambda(    3)=  0.2291   gamma=   35.83 GHz

     Calculation of q =    0.3750000  -0.3750000   0.6250000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    121                 3479     3479     941



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  4352 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_3/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.1

     total cpu time spent up to now is      477.6 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.3750000  -0.3750000   0.6250000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  4352

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A''  To be done


     PHONON       :  7m16.28s CPU     9m15.08s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn3

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 16.261583
     lambda(    1)=  0.0898   gamma=    5.03 GHz
     lambda(    2)=  0.0866   gamma=    6.79 GHz
     lambda(    3)=  0.1819   gamma=   42.57 GHz

     Calculation of q =    0.1250000  -0.1250000   0.3750000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    121                 3479     3479     893



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  4352 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_4/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.1

     total cpu time spent up to now is      627.5 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.1250000  -0.1250000   0.3750000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  4352

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A''  To be done


     PHONON       :  9m53.73s CPU    12m14.32s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn4

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 21.576766
     lambda(    1)=  0.2132   gamma=    6.84 GHz
     lambda(    2)=  0.1629   gamma=    5.61 GHz
     lambda(    3)=  0.3227   gamma=   33.92 GHz

     Calculation of q =   -0.1250000   0.6250000   0.1250000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    121                 3479     3479     941



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  4352 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_5/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.1

     total cpu time spent up to now is      767.4 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.1250000   0.6250000   0.1250000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  4352

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A''  To be done


     PHONON       : 12m23.95s CPU    14m59.14s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn5

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 13.351862
     lambda(    1)=  0.1040   gamma=    6.41 GHz
     lambda(    2)=  0.0721   gamma=    5.56 GHz
     lambda(    3)=  0.1982   gamma=   36.07 GHz

     Calculation of q =    0.6250000  -0.1250000   0.8750000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    139                 3479     3479    1067



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  8192 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_6/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.0

     total cpu time spent up to now is     1024.4 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.6250000  -0.1250000   0.8750000 )

     No symmetry!

     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  8192

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_1 (1)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A  To be done

     Representation     2      1 modes -A  To be done

     Representation     3      1 modes -A  To be done


     PHONON       : 16m51.61s CPU    19m55.90s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn6

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 15.162478
     lambda(    1)=  0.1291   gamma=   11.81 GHz
     lambda(    2)=  0.1171   gamma=   17.83 GHz
     lambda(    3)=  0.1580   gamma=   32.43 GHz

     Calculation of q =    0.3750000   0.1250000   0.6250000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    121                 3479     3479     941



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  8192 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_7/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.0

     total cpu time spent up to now is     1300.2 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.3750000   0.1250000   0.6250000 )

     No symmetry!

     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  8192

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_1 (1)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A  To be done

     Representation     2      1 modes -A  To be done

     Representation     3      1 modes -A  To be done


     PHONON       : 21m47.66s CPU    25m21.71s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn7

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 14.097346
     lambda(    1)=  0.0571   gamma=    3.75 GHz
     lambda(    2)=  0.0765   gamma=    7.98 GHz
     lambda(    3)=  0.1933   gamma=   39.58 GHz

     Calculation of q =   -0.1250000  -0.8750000   0.1250000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    127                 3479     3479     965



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  4352 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_8/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.1

     total cpu time spent up to now is     1458.5 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.1250000  -0.8750000   0.1250000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  4352

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A''  To be done


     PHONON       : 24m44.12s CPU    28m36.04s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn8

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 12.099956
     lambda(    1)=  0.0884   gamma=    7.70 GHz
     lambda(    2)=  0.0961   gamma=   10.22 GHz
     lambda(    3)=  0.1584   gamma=   36.81 GHz

     Calculation of q =   -0.3750000   0.3750000   0.3750000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    121                 3479     3479     941



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  1632 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_9/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.4

     total cpu time spent up to now is     1526.0 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.3750000   0.3750000   0.3750000 )

      6 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  1632

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      2 modes -E    L_3  To be done


     PHONON       : 25m59.64s CPU    30m 0.65s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn9

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 13.741609
     lambda(    1)=  0.0402   gamma=    1.88 GHz
     lambda(    2)=  0.0402   gamma=    1.88 GHz
     lambda(    3)=  0.2124   gamma=   44.40 GHz

     Calculation of q =    0.3750000  -0.3750000   1.1250000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input
     [opt_tetra]  Optimized tetrahedron method is used.

     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum         301     301    139                 3479     3479    1067



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      40.0000  Ry
     charge density cutoff     =     150.0000  Ry
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)

     celldm(1)=   7.628217  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=  4352 (tetrahedron method)

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     3479 G-vectors     FFT dimensions: (  24,  24,  24)

     Estimated max dynamical RAM per process >       1.80MB

     Estimated total allocated dynamical RAM >       7.21MB

     Check: negative/imaginary core charge=   -0.000013    0.000000

     The potential is recalculated from file :
     /mnt/c/Users/kawamuura/program/QE/qe_priv/tempdir/_ph0/aluminum.q_10/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 18.1

     total cpu time spent up to now is     1679.1 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     7.9706 ev

     Writing output data file aluminum.save
 
     [dfpt_tetra]  Dos(E_F)[/Ry] :   0.5441187E+01

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.6282  a.u.
     unit-cell volume          =     110.9709 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      40.0000  Ry
     charge density cut-off    =     150.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)


     celldm(1)=    7.62822  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.3750000  -0.3750000   1.1250000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  221.0943  (   3479 G-vectors)     FFT grid: ( 24, 24, 24)
     number of k points=  4352

     PseudoPot. # 1 for Al read from file:
     /mnt/c/Users/kawamuura/program/QE/qe_priv/pseudo/Al.pbe-n-rrkjus_psl.0.1.UPF
     MD5 check sum: 6479f6627750700cef6db36d98bbe09a
     Pseudo is Ultrasoft + core correction, Zval =  3.0
     Generated using "atomic" code by A. Dal Corso  v.5.0.99 svn rev. 10869
     Using radial grid of 1135 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A''  To be done


     PHONON       : 28m37.78s CPU    32m55.83s WALL

     Reading dVscf from file dv
     Reading dynamics matrix from file al.dyn10

     [elph_tetra]   Lowest band which contains FS :          2
     [elph_tetra]  Highest band which contains FS :          3
     [elph_tetra]    # of bands which contains FS :          2

     Tetrahedron method
     DOS =  2.720593 states/spin/Ry/Unit Cell at Ef=  7.970584 eV
     double delta at Ef = 13.069322
     lambda(    1)=  0.0426   gamma=    3.18 GHz
     lambda(    2)=  0.1310   gamma=   19.77 GHz
     lambda(    3)=  0.1601   gamma=   34.31 GHz

     init_run     :      7.28s CPU      9.55s WALL (      10 calls)
     electrons    :   1353.91s CPU   1571.89s WALL (      10 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.09s WALL (      10 calls)
     potinit      :      0.22s CPU      0.65s WALL (      10 calls)

     Called by electrons:
     c_bands      :   1334.50s CPU   1552.34s WALL (      10 calls)
     v_of_rho     :      0.19s CPU      0.24s WALL (      11 calls)
     newd         :      0.06s CPU      0.05s WALL (      11 calls)

     Called by c_bands:
     init_us_2    :     11.28s CPU     12.54s WALL (   45352 calls)
     cegterg      :   1182.03s CPU   1283.89s WALL (   15050 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :    989.28s CPU   1071.20s WALL (  248965 calls)
     s_psi        :     20.83s CPU     23.47s WALL (  248965 calls)
     g_psi        :      7.06s CPU      7.97s WALL (  222475 calls)
     cdiaghg      :     52.80s CPU     57.16s WALL (  233915 calls)

     Called by h_psi:
     h_psi:pot    :    985.09s CPU   1066.09s WALL (  248965 calls)
     h_psi:calbec :     23.08s CPU     25.44s WALL (  248965 calls)
     vloc_psi     :    939.47s CPU   1014.80s WALL (  248965 calls)
     add_vuspsi   :     21.20s CPU     23.28s WALL (  248965 calls)

     General routines
     calbec       :     29.95s CPU     32.30s WALL (  311885 calls)
     fft          :      0.08s CPU      0.11s WALL (     242 calls)
     ffts         :      8.20s CPU      8.70s WALL (   17190 calls)
     fftw         :   1031.47s CPU   1113.42s WALL ( 2430012 calls)
     davcio       :      0.42s CPU      6.67s WALL (   73604 calls)

     Parallel routines
     fft_scatter  :     54.89s CPU     63.41s WALL ( 2447444 calls)

     PHONON       : 28m54.03s CPU    33m13.30s WALL

     INITIALIZATION: 
     phq_setup    :      4.78s CPU      5.90s WALL (      10 calls)
     phq_init     :     98.33s CPU    111.81s WALL (      10 calls)

     phq_init     :     98.33s CPU    111.81s WALL (      10 calls)
     set_drhoc    :      2.14s CPU      2.29s WALL (      10 calls)
     init_vloc    :      0.03s CPU      0.03s WALL (      11 calls)
     init_us_1    :      5.52s CPU      6.53s WALL (      11 calls)
     newd         :      0.06s CPU      0.05s WALL (      11 calls)
     dvanqq       :      0.80s CPU      0.90s WALL (      20 calls)
     drho         :     89.75s CPU    101.94s WALL (      10 calls)




     dvqpsi_us    :    124.23s CPU    135.10s WALL (   17160 calls)
     addusddens   :      0.34s CPU      0.32s WALL (      28 calls)
     newdq        :      0.31s CPU      0.33s WALL (      28 calls)
     adddvscf     :      1.97s CPU      2.24s WALL (   17160 calls)

     dvqpsi_us    :    124.23s CPU    135.10s WALL (   17160 calls)
     dvqpsi_us_on :     12.05s CPU     13.11s WALL (   17160 calls)


     h_psi        :    989.28s CPU   1071.20s WALL (  248965 calls)

     h_psi        :    989.28s CPU   1071.20s WALL (  248965 calls)
     add_vuspsi   :     21.20s CPU     23.28s WALL (  248965 calls)

     addusdbec    :      2.53s CPU      2.81s WALL (   17160 calls)


      General routines
     calbec       :     29.95s CPU     32.30s WALL (  311885 calls)
     fft          :      0.08s CPU      0.11s WALL (     242 calls)
     ffts         :      8.20s CPU      8.70s WALL (   17190 calls)
     fftw         :   1031.47s CPU   1113.42s WALL ( 2430012 calls)
     davcio       :      0.42s CPU      6.67s WALL (   73604 calls)


     PHONON       : 28m54.03s CPU    33m13.30s WALL


   This run was terminated on:  21:38: 4  16Mar2017            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_21 = """
     Program PHONON v.6.0 (svn rev. 13286) starts on  7Feb2017 at 14:59:58 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     2 processors
     R & G space division:  proc/nbgrp/npool/nimage =       2

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/aluminum.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          60      60     21                  434      434      90
     Max          61      61     22                  435      435      91
     Sum         121     121     43                  869      869     181



     Dynamical matrices for ( 4, 4, 4)  uniform grid of q-points
     (   8q-points):
       N         xq(1)         xq(2)         xq(3) 
       1   0.000000000   0.000000000   0.000000000
       2  -0.250000000   0.250000000  -0.250000000
       3   0.500000000  -0.500000000   0.500000000
       4   0.000000000   0.500000000   0.000000000
       5   0.750000000  -0.250000000   0.750000000
       6   0.500000000   0.000000000   0.500000000
       7   0.000000000  -1.000000000   0.000000000
       8  -0.500000000  -1.000000000   0.000000000

     Calculation of q =    0.0000000   0.0000000   0.0000000

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   0.0000000 )

     49 Sym.Ops. (with q -> -q+G )


     G cutoff =   85.4897  (    435 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=    29  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, O_h (m-3m)  point group:


     Atomic displacements:
     There are   1 irreducible representations

     Representation     1      3 modes -T_1u G_15  G_4-  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     0.18s CPU         0.19s WALL



     Representation #  1 modes #   1  2  3

     Self-consistent Calculation

     Pert. #  1: Fermi energy shift (Ry) =    -8.2718E-24    -2.5077E-37
     Pert. #  2: Fermi energy shift (Ry) =    -1.2959E-23     3.6048E-37
     Pert. #  3: Fermi energy shift (Ry) =    -6.8932E-25     3.1347E-38

      iter #   1 total cpu time :     0.3 secs   av.it.:   3.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.257E-08

     Pert. #  1: Fermi energy shift (Ry) =     6.6174E-24    -1.2245E-40
     Pert. #  2: Fermi energy shift (Ry) =    -1.3786E-24     1.0408E-39
     Pert. #  3: Fermi energy shift (Ry) =     2.8951E-24     0.0000E+00

      iter #   2 total cpu time :     0.3 secs   av.it.:   5.5
      thresh= 1.121E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.314E-09

     Pert. #  1: Fermi energy shift (Ry) =     2.3713E-23    -1.6224E-39
     Pert. #  2: Fermi energy shift (Ry) =    -1.2959E-23     1.1020E-39
     Pert. #  3: Fermi energy shift (Ry) =     4.1359E-24    -4.2857E-40

      iter #   3 total cpu time :     0.4 secs   av.it.:   5.3
      thresh= 3.625E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.570E-13

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       0.171792 [THz] =       5.730349 [cm-1]
     freq (    2) =       0.171792 [THz] =       5.730349 [cm-1]
     freq (    3) =       0.171792 [THz] =       5.730349 [cm-1]
 **************************************************************************

     Mode symmetry, O_h (m-3m)  point group:

     freq (  1 -  3) =          5.7  [cm-1]   --> T_1u G_15  G_4- I  
     electron-phonon interaction  ...

     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  1.339210 states/spin/Ry/Unit Cell at Ef=  8.321793 eV
     lambda(    1)=  0.0000   gamma=    0.00 GHz
     lambda(    2)=  0.0000   gamma=    0.00 GHz
     lambda(    3)=  0.0000   gamma=    0.00 GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     DOS =  1.881761 states/spin/Ry/Unit Cell at Ef=  8.327153 eV
     lambda(    1)=  0.0000   gamma=    0.00 GHz
     lambda(    2)=  0.0000   gamma=    0.00 GHz
     lambda(    3)=  0.0000   gamma=    0.00 GHz
     Gaussian Broadening:   0.015 Ry, ngauss=   0
     DOS =  2.123229 states/spin/Ry/Unit Cell at Ef=  8.328621 eV
     lambda(    1)=  0.0000   gamma=    0.00 GHz
     lambda(    2)=  0.0000   gamma=    0.00 GHz
     lambda(    3)=  0.0000   gamma=    0.00 GHz
     Gaussian Broadening:   0.020 Ry, ngauss=   0
     DOS =  2.249739 states/spin/Ry/Unit Cell at Ef=  8.324319 eV
     lambda(    1)=  0.0000   gamma=    0.02 GHz
     lambda(    2)=  0.0000   gamma=    0.03 GHz
     lambda(    3)=  0.0000   gamma=    0.03 GHz
     Gaussian Broadening:   0.025 Ry, ngauss=   0
     DOS =  2.329803 states/spin/Ry/Unit Cell at Ef=  8.317861 eV
     lambda(    1)=  0.0000   gamma=    0.08 GHz
     lambda(    2)=  0.0000   gamma=    0.10 GHz
     lambda(    3)=  0.0000   gamma=    0.09 GHz
     Gaussian Broadening:   0.030 Ry, ngauss=   0
     DOS =  2.396029 states/spin/Ry/Unit Cell at Ef=  8.311296 eV
     lambda(    1)=  0.0000   gamma=    0.16 GHz
     lambda(    2)=  0.0000   gamma=    0.19 GHz
     lambda(    3)=  0.0000   gamma=    0.17 GHz
     Gaussian Broadening:   0.035 Ry, ngauss=   0
     DOS =  2.455226 states/spin/Ry/Unit Cell at Ef=  8.305262 eV
     lambda(    1)=  0.0000   gamma=    0.25 GHz
     lambda(    2)=  0.0000   gamma=    0.28 GHz
     lambda(    3)=  0.0000   gamma=    0.26 GHz
     Gaussian Broadening:   0.040 Ry, ngauss=   0
     DOS =  2.507873 states/spin/Ry/Unit Cell at Ef=  8.299956 eV
     lambda(    1)=  0.0000   gamma=    0.35 GHz
     lambda(    2)=  0.0000   gamma=    0.39 GHz
     lambda(    3)=  0.0000   gamma=    0.37 GHz
     Gaussian Broadening:   0.045 Ry, ngauss=   0
     DOS =  2.552966 states/spin/Ry/Unit Cell at Ef=  8.295411 eV
     lambda(    1)=  0.0000   gamma=    0.46 GHz
     lambda(    2)=  0.0000   gamma=    0.52 GHz
     lambda(    3)=  0.0000   gamma=    0.50 GHz
     Gaussian Broadening:   0.050 Ry, ngauss=   0
     DOS =  2.589582 states/spin/Ry/Unit Cell at Ef=  8.291553 eV
     lambda(    1)=  0.0000   gamma=    0.59 GHz
     lambda(    2)=  0.0000   gamma=    0.66 GHz
     lambda(    3)=  0.0000   gamma=    0.63 GHz


     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

     Calculation of q =   -0.2500000   0.2500000  -0.2500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          60      60     27                  434      434     129
     Max          61      61     28                  435      435     130
     Sum         121     121     55                  869      869     259



     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =         3.00
     number of Kohn-Sham states=            6
     kinetic-energy cutoff     =      15.0000  Ry
     charge density cutoff     =      60.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   7.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   240  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:      869 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.22MB

     Estimated total allocated dynamical RAM >       0.43MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/aluminum.q_2/aluminum.save/charge-density.dat

     Starting wfc are    4 atomic +    2 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  3.33E-10,  avg # of iterations = 13.7

     total cpu time spent up to now is        0.9 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is     8.1776 ev

     Writing output data file aluminum.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.5000  a.u.
     unit-cell volume          =     105.4688 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      15.0000  Ry
     charge density cut-off    =      60.0000  Ry
     convergence threshold     =      1.0E-10
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    7.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (  -0.2500000   0.2500000  -0.2500000 )

      6 Sym.Ops. (no q -> -q+G )


     G cutoff =   85.4897  (    435 G-vectors)     FFT grid: ( 15, 15, 15)

     number of k points=   240  Marzari-Vanderbilt smearing, width (Ry)=  0.0500

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      2 modes -E    L_3  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     2.83s CPU         2.91s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     3.0 secs   av.it.:   4.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.094E-02

      iter #   2 total cpu time :     3.0 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  9.107E-01

      iter #   3 total cpu time :     3.1 secs   av.it.:   4.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.174E-07

      iter #   4 total cpu time :     3.2 secs   av.it.:   5.2
      thresh= 7.193E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.406E-09

     Maximum CPU time exceeded

     max_seconds     =       3.00
     elapsed seconds =       3.03

     PHONON       :     3.04s CPU         3.18s WALL

     INITIALIZATION: 
     phq_setup    :      0.00s CPU      0.00s WALL (       2 calls)
     phq_init     :      0.02s CPU      0.02s WALL (       2 calls)

     phq_init     :      0.02s CPU      0.02s WALL (       2 calls)
     init_vloc    :      0.00s CPU      0.00s WALL (       2 calls)
     init_us_1    :      0.00s CPU      0.00s WALL (       2 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.01s CPU      0.01s WALL (       2 calls)
     phqscf       :      0.40s CPU      0.48s WALL (       2 calls)
     dynmatrix    :      0.00s CPU      0.00s WALL (       1 calls)

     phqscf       :      0.40s CPU      0.48s WALL (       3 calls)
     solve_linter :      0.39s CPU      0.48s WALL (       2 calls)
     drhodv       :      0.00s CPU      0.00s WALL (       1 calls)

     dynmat0      :      0.01s CPU      0.01s WALL (       2 calls)
     dynmat_us    :      0.01s CPU      0.01s WALL (       2 calls)
     d2ionq       :      0.00s CPU      0.00s WALL (       2 calls)

     dynmat_us    :      0.01s CPU      0.01s WALL (       2 calls)

     phqscf       :      0.40s CPU      0.48s WALL (       4 calls)
     solve_linter :      0.39s CPU      0.48s WALL (       3 calls)

     solve_linter :      0.39s CPU      0.48s WALL (       4 calls)
     dvqpsi_us    :      0.03s CPU      0.04s WALL (     207 calls)
     ortho        :      0.00s CPU      0.01s WALL (     741 calls)
     cgsolve      :      0.23s CPU      0.26s WALL (     741 calls)
     incdrhoscf   :      0.02s CPU      0.03s WALL (     741 calls)
     vpsifft      :      0.01s CPU      0.03s WALL (     534 calls)
     dv_of_drho   :      0.00s CPU      0.00s WALL (      13 calls)
     mix_pot      :      0.00s CPU      0.00s WALL (       7 calls)
     ef_shift     :      0.00s CPU      0.00s WALL (       4 calls)
     localdos     :      0.00s CPU      0.00s WALL (       1 calls)
     psymdvscf    :      0.04s CPU      0.04s WALL (       7 calls)

     dvqpsi_us    :      0.03s CPU      0.04s WALL (     207 calls)
     dvqpsi_us_on :      0.00s CPU      0.00s WALL (     207 calls)

     cgsolve      :      0.23s CPU      0.26s WALL (     741 calls)
     ch_psi       :      0.22s CPU      0.24s WALL (    4042 calls)

     ch_psi       :      0.22s CPU      0.24s WALL (    4042 calls)
     h_psi        :      0.55s CPU      0.64s WALL (    7830 calls)
     last         :      0.02s CPU      0.03s WALL (    4042 calls)

     h_psi        :      0.55s CPU      0.64s WALL (    7830 calls)
     add_vuspsi   :      0.01s CPU      0.02s WALL (    7830 calls)

     incdrhoscf   :      0.02s CPU      0.03s WALL (     741 calls)


      General routines
     calbec       :      0.02s CPU      0.04s WALL (   13739 calls)
     fft          :      0.00s CPU      0.00s WALL (      65 calls)
     ffts         :      0.00s CPU      0.00s WALL (     262 calls)
     fftw         :      0.46s CPU      0.56s WALL (   53334 calls)
     davcio       :      0.02s CPU      0.02s WALL (    4193 calls)
     write_rec    :      0.01s CPU      0.01s WALL (       8 calls)


     PHONON       :     3.04s CPU         3.18s WALL


   This run was terminated on:  15: 0: 1   7Feb2017            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
-------------------------------------------------------
Primary job  terminated normally, but 1 process returned
a non-zero exit code.. Per user-direction, the job has been aborted.
-------------------------------------------------------
--------------------------------------------------------------------------
mpirun detected that one or more processes exited with non-zero status, thus causing
the job to be terminated. The first process to do so was:

  Process name: [[33832,1],1]
  Exit code:    1
--------------------------------------------------------------------------
"""


test_26 = """
     Program PHONON v.6.0 (svn rev. 13188M) starts on  7Dec2016 at 13:16:14 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     R & G space division:  proc/nbgrp/npool/nimage =       4

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/nickel.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PW   PBE  PBE ( 1  4  3  4 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want

               file Ni.pbe-nd-rrkjus.UPF: wavefunction(s)  4S renormalized

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         112      40     15                 1604      351      82
     Max         113      41     16                 1607      354      83
     Sum         451     163     61                 6423     1411     331


     Check: negative/imaginary core charge=   -0.000020    0.000000

     Calculation of q =    0.0000000   0.0000000   1.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         112      40     21                 1604      351     132
     Max         113      41     22                 1607      354     135
     Sum         451     163     85                 6423     1411     531


     Title: 
     phonons of Ni at X                                                         


     bravais-lattice index     =            2
     lattice parameter (alat)  =       6.6500  a.u.
     unit-cell volume          =      73.5199 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =        10.00
     number of Kohn-Sham states=           18
     kinetic-energy cutoff     =      27.0000  Ry
     charge density cutoff     =     300.0000  Ry
     Exchange-correlation      =  SLA  PW   PBE  PBE ( 1  4  3  4 0 0)
     Noncollinear calculation without spin-orbit


     celldm(1)=   6.650000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Ni read from file:
     ./Ni.pbe-nd-rrkjus.UPF
     MD5 check sum: 8081f0a005c9a5470caab1a58e82ecb2
     Pseudo is Ultrasoft + core correction, Zval = 10.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of 1203 points,  6 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
                l(5) =   2
                l(6) =   2
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Ni            10.00    58.69340     Ni( 1.00)

     16 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Ni  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=   216  Marzari-Vanderbilt smearing, width (Ry)=  0.0200

     Number of k-points >= 100: set verbosity='high' to print them.

     Dense  grid:     6423 G-vectors     FFT dimensions: (  27,  27,  27)

     Smooth grid:     1411 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       1.61Mb

     Estimated total allocated dynamical RAM >       6.44Mb

     Check: negative/imaginary core charge=   -0.000020    0.000000

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/nickel.save/charge-density.dat

     Starting wfc are   12 atomic +    6 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.00E-10,  avg # of iterations = 14.2

     total cpu time spent up to now is       18.4 secs

     End of band structure calculation

     Number of k-points >= 100: set verbosity='high' to print the bands.

     the Fermi energy is    14.2603 ev

     Writing output data file nickel.save

     Fixed quantization axis for GGA:    -0.000000    1.000000   -0.000000

     phonons of Ni at X                                                         

     bravais-lattice index     =            2
     lattice parameter (alat)  =       6.6500  a.u.
     unit-cell volume          =      73.5199 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      27.0000  Ry
     charge density cut-off    =     300.0000  Ry
     convergence threshold     =      1.0E-16
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PW   PBE  PBE ( 1  4  3  4 0 0)
     Noncollinear calculation without spin-orbit

     celldm(1)=    6.65000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Ni  58.6934   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   1.0000000 )

      8 Sym.Ops. (no q -> -q+G )


     G cutoff =  336.0507  (   1607 G-vectors)     FFT grid: ( 27, 27, 27)
     G cutoff =  120.9783  (    354 G-vectors)  smooth grid: ( 15, 15, 15)

     number of k points=   216  Marzari-Vanderbilt smearing, width (Ry)=  0.0200

     PseudoPot. # 1 for Ni read from file:
     ./Ni.pbe-nd-rrkjus.UPF
     MD5 check sum: 8081f0a005c9a5470caab1a58e82ecb2
     Pseudo is Ultrasoft + core correction, Zval = 10.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of 1203 points,  6 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
                l(5) =   2
                l(6) =   2
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, D_2h (mmm)  point group:


     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      1 modes -B_1u  To be done

     Representation     2      1 modes -B_2u  To be done

     Representation     3      1 modes -B_3u  To be done



     Alpha used in Ewald sum =   2.8000
     PHONON       :    18.15s CPU        23.25s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :    25.7 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.476E-05

      iter #   2 total cpu time :    28.1 secs   av.it.:   7.7
      thresh= 7.400E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.185E-05

      iter #   3 total cpu time :    30.7 secs   av.it.:   7.2
      thresh= 5.643E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.082E-09

      iter #   4 total cpu time :    33.3 secs   av.it.:   7.3
      thresh= 6.389E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.666E-11

      iter #   5 total cpu time :    35.5 secs   av.it.:   6.6
      thresh= 6.055E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.176E-12

      iter #   6 total cpu time :    37.8 secs   av.it.:   7.1
      thresh= 1.085E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.548E-15

      iter #   7 total cpu time :    40.2 secs   av.it.:   7.3
      thresh= 7.448E-09 alpha_mix =  0.700 |ddv_scf|^2 =  7.637E-17

      iter #   8 total cpu time :    42.6 secs   av.it.:   7.2
      thresh= 8.739E-10 alpha_mix =  0.700 |ddv_scf|^2 =  1.260E-17

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :    44.9 secs   av.it.:   4.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.625E-06

      iter #   2 total cpu time :    47.4 secs   av.it.:   7.9
      thresh= 2.574E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.584E-07

      iter #   3 total cpu time :    49.9 secs   av.it.:   7.7
      thresh= 5.083E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.503E-09

      iter #   4 total cpu time :    52.2 secs   av.it.:   7.1
      thresh= 3.877E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.603E-12

      iter #   5 total cpu time :    54.6 secs   av.it.:   7.2
      thresh= 1.613E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.865E-14

      iter #   6 total cpu time :    57.0 secs   av.it.:   7.4
      thresh= 2.422E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.505E-15

      iter #   7 total cpu time :    59.5 secs   av.it.:   8.0
      thresh= 5.005E-09 alpha_mix =  0.700 |ddv_scf|^2 =  5.891E-17

      iter #   8 total cpu time :    62.3 secs   av.it.:   7.9
      thresh= 7.676E-10 alpha_mix =  0.700 |ddv_scf|^2 =  2.594E-18

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :    64.8 secs   av.it.:   4.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.625E-06

      iter #   2 total cpu time :    67.3 secs   av.it.:   7.9
      thresh= 2.574E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.584E-07

      iter #   3 total cpu time :    70.1 secs   av.it.:   7.7
      thresh= 5.083E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.503E-09

      iter #   4 total cpu time :    72.5 secs   av.it.:   7.1
      thresh= 3.877E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.603E-12

      iter #   5 total cpu time :    74.9 secs   av.it.:   7.2
      thresh= 1.613E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.865E-14

      iter #   6 total cpu time :    77.3 secs   av.it.:   7.4
      thresh= 2.422E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.506E-15

      iter #   7 total cpu time :    79.8 secs   av.it.:   8.0
      thresh= 5.006E-09 alpha_mix =  0.700 |ddv_scf|^2 =  5.893E-17

      iter #   8 total cpu time :    82.3 secs   av.it.:   7.9
      thresh= 7.677E-10 alpha_mix =  0.700 |ddv_scf|^2 =  2.679E-18

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    2
     List of q in the star:
          1   0.000000000   0.000000000   1.000000000
          2   1.000000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   1.000000000 ) 

 **************************************************************************
     freq (    1) =       6.382312 [THz] =     212.891027 [cm-1]
     freq (    2) =       6.382357 [THz] =     212.892523 [cm-1]
     freq (    3) =       8.906622 [THz] =     297.092916 [cm-1]
 **************************************************************************

     Mode symmetry, D_2h (mmm)  [C_2h (2/m) ] magnetic point group:

     freq (  1 -  1) =        212.9  [cm-1]   --> B_2u               
     freq (  2 -  2) =        212.9  [cm-1]   --> B_3u               
     freq (  3 -  3) =        297.1  [cm-1]   --> B_1u               

     init_run     :      0.24s CPU      0.35s WALL (       1 calls)
     electrons    :     14.23s CPU     18.03s WALL (       1 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       1 calls)
     potinit      :      0.03s CPU      0.03s WALL (       1 calls)

     Called by electrons:
     c_bands      :     14.22s CPU     18.02s WALL (       1 calls)
     v_of_rho     :      0.04s CPU      0.05s WALL (       2 calls)
     newd         :      0.02s CPU      0.03s WALL (       2 calls)

     Called by c_bands:
     init_us_2    :      0.12s CPU      0.16s WALL (    3564 calls)
     cegterg      :     13.58s CPU     17.15s WALL (     223 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :     26.25s CPU     35.84s WALL (   26146 calls)
     s_psi        :      3.88s CPU      4.71s WALL (   51384 calls)
     g_psi        :      0.02s CPU      0.03s WALL (    3061 calls)
     cdiaghg      :      6.18s CPU      7.42s WALL (    3277 calls)

     Called by h_psi:
     h_psi:pot    :     26.06s CPU     35.67s WALL (   26146 calls)
     h_psi:calbec :      2.20s CPU      3.06s WALL (   26146 calls)
     vloc_psi     :     21.81s CPU     29.98s WALL (   26146 calls)
     add_vuspsi   :      1.96s CPU      2.54s WALL (   26146 calls)

     General routines
     calbec       :      5.23s CPU      7.47s WALL (   57108 calls)
     fft          :      0.37s CPU      0.71s WALL (    1264 calls)
     ffts         :      0.01s CPU      0.03s WALL (     548 calls)
     fftw         :     19.91s CPU     28.91s WALL ( 1249472 calls)
     interpolate  :      0.06s CPU      0.07s WALL (     212 calls)
     davcio       :      0.18s CPU      0.31s WALL (   14046 calls)

     Parallel routines
     fft_scatter  :     10.70s CPU     18.24s WALL ( 1251284 calls)

     PHONON       :  1m 1.50s CPU     1m22.48s WALL

     INITIALIZATION: 
     phq_setup    :      0.10s CPU      0.12s WALL (       1 calls)
     phq_init     :      2.80s CPU      3.75s WALL (       1 calls)

     phq_init     :      2.80s CPU      3.75s WALL (       1 calls)
     set_drhoc    :      0.16s CPU      0.20s WALL (       3 calls)
     init_vloc    :      0.01s CPU      0.01s WALL (       2 calls)
     init_us_1    :      0.37s CPU      0.53s WALL (       2 calls)
     newd         :      0.02s CPU      0.03s WALL (       2 calls)
     dvanqq       :      0.11s CPU      0.16s WALL (       1 calls)
     drho         :      1.97s CPU      2.64s WALL (       1 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.55s CPU      0.72s WALL (       1 calls)
     phqscf       :     43.34s CPU     59.22s WALL (       1 calls)
     dynmatrix    :      0.00s CPU      0.01s WALL (       1 calls)

     phqscf       :     43.34s CPU     59.22s WALL (       1 calls)
     solve_linter :     42.96s CPU     58.74s WALL (       3 calls)
     drhodv       :      0.38s CPU      0.47s WALL (       3 calls)

     dynmat0      :      0.55s CPU      0.72s WALL (       1 calls)
     dynmat_us    :      0.42s CPU      0.58s WALL (       1 calls)
     d2ionq       :      0.00s CPU      0.00s WALL (       1 calls)
     dynmatcc     :      0.12s CPU      0.14s WALL (       1 calls)

     dynmat_us    :      0.42s CPU      0.58s WALL (       1 calls)
     addusdynmat  :      0.00s CPU      0.00s WALL (       1 calls)

     phqscf       :     43.34s CPU     59.22s WALL (       1 calls)
     solve_linter :     42.96s CPU     58.74s WALL (       3 calls)

     solve_linter :     42.96s CPU     58.74s WALL (       3 calls)
     dvqpsi_us    :      1.36s CPU      1.68s WALL (     324 calls)
     ortho        :      1.49s CPU      2.26s WALL (    2592 calls)
     cgsolve      :     31.52s CPU     43.36s WALL (    2592 calls)
     incdrhoscf   :      3.14s CPU      4.07s WALL (    2592 calls)
     addusddens   :      0.58s CPU      0.66s WALL (      27 calls)
     vpsifft      :      1.85s CPU      2.83s WALL (    2268 calls)
     dv_of_drho   :      0.32s CPU      0.43s WALL (      24 calls)
     mix_pot      :      0.04s CPU      0.09s WALL (      24 calls)
     psymdvscf    :      0.31s CPU      0.36s WALL (      24 calls)
     newdq        :      0.46s CPU      0.55s WALL (      24 calls)
     adddvscf     :      0.54s CPU      0.70s WALL (    2268 calls)
     drhodvus     :      0.00s CPU      0.00s WALL (       3 calls)

     dvqpsi_us    :      1.36s CPU      1.68s WALL (     324 calls)
     dvqpsi_us_on :      1.02s CPU      1.07s WALL (     324 calls)

     cgsolve      :     31.52s CPU     43.36s WALL (    2592 calls)
     ch_psi       :     30.22s CPU     41.31s WALL (   22646 calls)

     ch_psi       :     30.22s CPU     41.31s WALL (   22646 calls)
     h_psi        :     26.25s CPU     35.84s WALL (   26146 calls)
     last         :      6.29s CPU      8.84s WALL (   22646 calls)

     h_psi        :     26.25s CPU     35.84s WALL (   26146 calls)
     add_vuspsi   :      1.96s CPU      2.54s WALL (   26146 calls)

     incdrhoscf   :      3.14s CPU      4.07s WALL (    2592 calls)

     drhodvus     :      0.00s CPU      0.00s WALL (       3 calls)

      General routines
     calbec       :      5.23s CPU      7.47s WALL (   57108 calls)
     fft          :      0.37s CPU      0.71s WALL (    1264 calls)
     ffts         :      0.01s CPU      0.03s WALL (     548 calls)
     fftw         :     19.91s CPU     28.91s WALL ( 1249472 calls)
     davcio       :      0.18s CPU      0.31s WALL (   14046 calls)
     write_rec    :      0.04s CPU      0.06s WALL (      27 calls)


     PHONON       :  1m 1.50s CPU     1m22.49s WALL


   This run was terminated on:  13:17:36   7Dec2016            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_31 = """
     Program PHONON v.6.4.1 starts on 20Sep2019 at 16:13:29 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
         "P. Giannozzi et al., J. Phys.:Condens. Matter 29 465901 (2017);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors

     MPI processes distributed on     1 nodes
     R & G space division:  proc/nbgrp/npool/nimage =       4

     Reading data from directory:
     /benchmarks/tempdir/graphite.save/

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation= RVV10
                           (   1   4  13   4   3   0   0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want

               file C.pbe-rrkjus.UPF: wavefunction(s)  2S 2P renormalized
 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          65      46     15                 2276     1235     259
     Max          67      47     16                 2279     1237     264
     Sum         265     187     61                 9111     4943    1045
 


---------------------------------------------------------------------------------
Carrying out rVV10 run using the following parameters:
Nqs =      20    Nr_points =    1024   r_max =   100.000
b_value =  6.30000    beta =  0.00901
     q_mesh =  0.00010000  0.00030000  0.00058939  0.00100810
               0.00161396  0.00249058  0.00375900  0.00559430
               0.00824984  0.01209221  0.01765183  0.02569619
               0.03733578  0.05417739  0.07854596  0.11380545
               0.16482331  0.23864234  0.34545298  0.50000000

Gradients computed in Reciprocal space

---------------------------------------------------------------------------------



     Calculation of q =    0.3333333   0.5773503   0.0000000

     Subspace diagonalization in iterative solution of the eigenvalue problem:
     a serial algorithm will be used

 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          65      46     21                 2276     1235     385
     Max          67      47     22                 2279     1237     392
     Sum         265     187     85                 9111     4943    1553
 

     Title: 
     phonons of graphite                                                        


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.6463  a.u.
     unit-cell volume          =     224.3370 (a.u.)^3
     number of atoms/cell      =            4
     number of atomic types    =            1
     number of electrons       =        16.00
     number of Kohn-Sham states=            8
     kinetic-energy cutoff     =      30.0000  Ry
     charge density cutoff     =     180.0000  Ry
     Exchange-correlation= RVV10
                           (   1   4  13   4   3   0   0)

     celldm(1)=   4.646303  celldm(2)=   0.000000  celldm(3)=   2.582543
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   2.582543 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.387215 )  


     PseudoPot. # 1 for C  read from file:
     /benchmarks/pseudo/C.pbe-rrkjus.UPF
     MD5 check sum: c9ac5a99bc85b198593446162950cd17
     Pseudo is Ultrasoft, Zval =  4.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  627 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        C              4.00    12.01070     C ( 1.00)

     24 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           C   tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           C   tau(   2) = (  -0.0000000   0.5773503   0.0000000  )
         3           C   tau(   3) = (   0.0000000   0.0000000   1.2912714  )
         4           C   tau(   4) = (   0.5000000   0.2886751   1.2912714  )

     number of k points=    20
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.1250000   0.2165064   0.0968038), wk =   0.1250000
        k(    2) = (   0.4583333   0.7938566   0.0968038), wk =   0.0000000
        k(    3) = (   0.1250000   0.5051815   0.0968038), wk =   0.2500000
        k(    4) = (   0.4583333   1.0825318   0.0968038), wk =   0.0000000
        k(    5) = (   0.1250000  -0.3608439   0.0968038), wk =   0.2500000
        k(    6) = (   0.4583333   0.2165064   0.0968038), wk =   0.0000000
        k(    7) = (   0.1250000  -0.0721688   0.0968038), wk =   0.2500000
        k(    8) = (   0.4583333   0.5051815   0.0968038), wk =   0.0000000
        k(    9) = (   0.3750000   0.6495191   0.0968038), wk =   0.1250000
        k(   10) = (   0.7083333   1.2268693   0.0968038), wk =   0.0000000
        k(   11) = (   0.3750000  -0.2165064   0.0968038), wk =   0.2500000
        k(   12) = (   0.7083333   0.3608439   0.0968038), wk =   0.0000000
        k(   13) = (  -0.1250000  -0.2165064  -0.0968038), wk =   0.1250000
        k(   14) = (   0.2083333   0.3608439  -0.0968038), wk =   0.0000000
        k(   15) = (  -0.1250000  -0.5051815  -0.0968038), wk =   0.2500000
        k(   16) = (   0.2083333   0.0721688  -0.0968038), wk =   0.0000000
        k(   17) = (  -0.1250000   0.3608439  -0.0968038), wk =   0.2500000
        k(   18) = (   0.2083333   0.9381942  -0.0968038), wk =   0.0000000
        k(   19) = (  -0.3750000  -0.6495191  -0.0968038), wk =   0.1250000
        k(   20) = (  -0.0416667  -0.0721688  -0.0968038), wk =   0.0000000

     Dense  grid:     9111 G-vectors     FFT dimensions: (  20,  20,  54)

     Smooth grid:     4943 G-vectors     FFT dimensions: (  18,  18,  48)

     Estimated max dynamical RAM per process >       2.66 MB

     Estimated total dynamical RAM >      10.65 MB

     The potential is recalculated from file :
     /benchmarks/tempdir/_ph0/graphite.save/charge-density

     Starting wfcs are   16 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  6.25E-11,  avg # of iterations = 10.5

     total cpu time spent up to now is        0.4 secs

     End of band structure calculation

          k = 0.1250 0.2165 0.0968 (   619 PWs)   bands (ev):

   -10.8738 -10.5682   0.5104   1.4221   1.4730   2.2699   2.5058   2.5418

          k = 0.4583 0.7939 0.0968 (   625 PWs)   bands (ev):

    -6.5480  -6.4244  -5.4692  -5.3773   0.6875   0.7674   5.0640   5.7350

          k = 0.1250 0.5052 0.0968 (   628 PWs)   bands (ev):

    -7.3284  -7.0965  -4.4034  -4.2265   0.2017   0.2943   4.5259   5.6756

          k = 0.4583 1.0825 0.0968 (   629 PWs)   bands (ev):

    -8.2354  -7.9803  -2.8330  -2.6889  -0.1756  -0.0803   3.6011   5.0485

          k = 0.1250-0.3608 0.0968 (   627 PWs)   bands (ev):

    -9.4169  -9.1387  -1.3886  -1.2645   0.9751   1.0527   2.2401   3.8534

          k = 0.4583 0.2165 0.0968 (   626 PWs)   bands (ev):

    -7.5528  -7.3162  -4.5881  -4.4043   1.1686   1.2373   4.2400   5.3927

          k = 0.1250-0.0722 0.0968 (   606 PWs)   bands (ev):

   -11.6140 -11.2945  -0.3831   1.4333   3.2129   3.2687   3.8907   3.9321

          k = 0.4583 0.5052 0.0968 (   630 PWs)   bands (ev):

    -6.8845  -6.6621  -4.0845  -3.9286  -1.2384  -1.1273   5.0799   6.2210

          k = 0.3750 0.6495 0.0968 (   630 PWs)   bands (ev):

    -5.7021  -5.5877  -4.7720  -4.7236  -2.0978  -1.9840   6.2474   6.8892

          k = 0.7083 1.2269 0.0968 (   630 PWs)   bands (ev):

    -6.2421  -6.0363  -4.0214  -3.9251  -2.3286  -2.2667   5.8289   6.9067

          k = 0.3750-0.2165 0.0968 (   623 PWs)   bands (ev):

    -8.7020  -8.4377  -3.0948  -2.9278   1.5762   1.6396   3.0391   4.5238

          k = 0.7083 0.3608 0.0968 (   623 PWs)   bands (ev):

    -9.6566  -9.3738  -1.4534  -1.3139   1.8791   1.9430   1.9532   3.5881

          k =-0.1250-0.2165-0.0968 (   619 PWs)   bands (ev):

   -10.8738 -10.5682   0.5104   1.4221   1.4730   2.2699   2.5058   2.5418

          k = 0.2083 0.3608-0.0968 (   625 PWs)   bands (ev):

    -8.9405  -8.6714  -1.5766  -1.5022  -0.0608  -0.0052   2.8016   4.3630

          k =-0.1250-0.5052-0.0968 (   628 PWs)   bands (ev):

    -7.3284  -7.0965  -4.4034  -4.2265   0.2017   0.2943   4.5259   5.6756

          k = 0.2083 0.0722-0.0968 (   613 PWs)   bands (ev):

   -11.1197 -10.8095   0.2142   1.7777   1.8581   1.9942   3.1268   3.1771

          k =-0.1250 0.3608-0.0968 (   627 PWs)   bands (ev):

    -9.4169  -9.1387  -1.3886  -1.2645   0.9751   1.0527   2.2401   3.8534

          k = 0.2083 0.9382-0.0968 (   621 PWs)   bands (ev):

   -10.3843 -10.0879   0.1808   0.2828   1.0967   2.1148   2.1766   2.8114

          k =-0.3750-0.6495-0.0968 (   630 PWs)   bands (ev):

    -5.7021  -5.5877  -4.7720  -4.7236  -2.0978  -1.9840   6.2474   6.8892

          k =-0.0417-0.0722-0.0968 (   594 PWs)   bands (ev):

   -11.8622 -11.5380  -0.6844   1.1498   4.1138   4.1409   4.3315   4.3552

     highest occupied level (ev):     6.8892

     Writing output data file /benchmarks/tempdir/_ph0/graphite.save/

     phonons of graphite                                                        

     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.6463  a.u.
     unit-cell volume          =     224.3370 (a.u.)^3
     number of atoms/cell      =            4
     number of atomic types    =            1
     kinetic-energy cut-off    =      30.0000  Ry
     charge density cut-off    =     180.0000  Ry
     convergence threshold     =      1.0E-18
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation= RVV10
                           (   1   4  13   4   3   0   0)


     celldm(1)=    4.64630  celldm(2)=    0.00000  celldm(3)=    2.58254
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  1.0000  0.0000  0.0000 )  
               a(2) = ( -0.5000  0.8660  0.0000 )  
               a(3) = (  0.0000  0.0000  2.5825 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.0000  0.5774  0.0000 )  
               b(2) = (  0.0000  1.1547  0.0000 )  
               b(3) = (  0.0000  0.0000  0.3872 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     C   12.0107   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     C   12.0107   tau(    2) = (   -0.00000    0.57735    0.00000  )
        3     C   12.0107   tau(    3) = (    0.00000    0.00000    1.29127  )
        4     C   12.0107   tau(    4) = (    0.50000    0.28868    1.29127  )

     Computing dynamical matrix for 
                    q = (   0.3333333   0.5773503   0.0000000 )
 
     12 Sym.Ops. (no q -> -q+G )


     G cutoff =   98.4301  (   2276 G-vectors)     FFT grid: ( 20, 20, 54)
     G cutoff =   65.6201  (   1236 G-vectors)  smooth grid: ( 18, 18, 48)
     number of k points=    20
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.1250000   0.2165064   0.0968038), wk =   0.1250000
        k(    2) = (   0.4583333   0.7938566   0.0968038), wk =   0.0000000
        k(    3) = (   0.1250000   0.5051815   0.0968038), wk =   0.2500000
        k(    4) = (   0.4583333   1.0825318   0.0968038), wk =   0.0000000
        k(    5) = (   0.1250000  -0.3608439   0.0968038), wk =   0.2500000
        k(    6) = (   0.4583333   0.2165064   0.0968038), wk =   0.0000000
        k(    7) = (   0.1250000  -0.0721688   0.0968038), wk =   0.2500000
        k(    8) = (   0.4583333   0.5051815   0.0968038), wk =   0.0000000
        k(    9) = (   0.3750000   0.6495191   0.0968038), wk =   0.1250000
        k(   10) = (   0.7083333   1.2268693   0.0968038), wk =   0.0000000
        k(   11) = (   0.3750000  -0.2165064   0.0968038), wk =   0.2500000
        k(   12) = (   0.7083333   0.3608439   0.0968038), wk =   0.0000000
        k(   13) = (  -0.1250000  -0.2165064  -0.0968038), wk =   0.1250000
        k(   14) = (   0.2083333   0.3608439  -0.0968038), wk =   0.0000000
        k(   15) = (  -0.1250000  -0.5051815  -0.0968038), wk =   0.2500000
        k(   16) = (   0.2083333   0.0721688  -0.0968038), wk =   0.0000000
        k(   17) = (  -0.1250000   0.3608439  -0.0968038), wk =   0.2500000
        k(   18) = (   0.2083333   0.9381942  -0.0968038), wk =   0.0000000
        k(   19) = (  -0.3750000  -0.6495191  -0.0968038), wk =   0.1250000
        k(   20) = (  -0.0416667  -0.0721688  -0.0968038), wk =   0.0000000

     PseudoPot. # 1 for C  read from file:
     /benchmarks/pseudo/C.pbe-rrkjus.UPF
     MD5 check sum: c9ac5a99bc85b198593446162950cd17
     Pseudo is Ultrasoft, Zval =  4.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  627 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, D_3h (-62m) point group:


     Atomic displacements:
     There are   8 irreducible representations

     Representation     1      1 modes -A'_1  To be done

     Representation     2      1 modes -A'_2  To be done

     Representation     3      2 modes -E'  To be done

     Representation     4      2 modes -E'  To be done

     Representation     5      2 modes -E'  To be done

     Representation     6      1 modes -A''1  To be done

     Representation     7      1 modes -A''2  To be done

     Representation     8      2 modes -E''  To be done



     Alpha used in Ewald sum =   1.8000
     PHONON       :      2.32s CPU      2.54s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     2.7 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.226E-05

      iter #   2 total cpu time :     3.0 secs   av.it.:   8.1
      thresh= 5.680E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.838E-05

      iter #   3 total cpu time :     3.2 secs   av.it.:   7.5
      thresh= 6.195E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.760E-08

      iter #   4 total cpu time :     3.4 secs   av.it.:   7.5
      thresh= 1.939E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.839E-10

      iter #   5 total cpu time :     3.6 secs   av.it.:   7.1
      thresh= 1.685E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.508E-13

      iter #   6 total cpu time :     3.9 secs   av.it.:   8.2
      thresh= 5.922E-08 alpha_mix =  0.700 |ddv_scf|^2 =  5.382E-14

      iter #   7 total cpu time :     4.1 secs   av.it.:   7.8
      thresh= 2.320E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.883E-14

      iter #   8 total cpu time :     4.3 secs   av.it.:   7.6
      thresh= 1.698E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.513E-16

      iter #   9 total cpu time :     4.5 secs   av.it.:   7.7
      thresh= 1.230E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.921E-17

      iter #  10 total cpu time :     4.7 secs   av.it.:   7.6
      thresh= 5.405E-10 alpha_mix =  0.700 |ddv_scf|^2 =  3.097E-19

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     5.0 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.265E-05

      iter #   2 total cpu time :     5.2 secs   av.it.:   8.2
      thresh= 4.759E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.174E-05

      iter #   3 total cpu time :     5.4 secs   av.it.:   7.6
      thresh= 4.663E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.360E-08

      iter #   4 total cpu time :     5.6 secs   av.it.:   7.7
      thresh= 1.833E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.497E-10

      iter #   5 total cpu time :     5.8 secs   av.it.:   7.2
      thresh= 1.580E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.790E-13

      iter #   6 total cpu time :     6.0 secs   av.it.:   7.6
      thresh= 6.156E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.440E-14

      iter #   7 total cpu time :     6.3 secs   av.it.:   8.2
      thresh= 1.562E-08 alpha_mix =  0.700 |ddv_scf|^2 =  1.663E-14

      iter #   8 total cpu time :     6.5 secs   av.it.:   7.6
      thresh= 1.289E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.800E-17

      iter #   9 total cpu time :     6.7 secs   av.it.:   7.4
      thresh= 6.928E-10 alpha_mix =  0.700 |ddv_scf|^2 =  4.215E-18

      iter #  10 total cpu time :     7.0 secs   av.it.:   7.9
      thresh= 2.053E-10 alpha_mix =  0.700 |ddv_scf|^2 =  8.588E-20

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 modes #   3  4

     Self-consistent Calculation

      iter #   1 total cpu time :     7.3 secs   av.it.:   5.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.570E-06

      iter #   2 total cpu time :     7.8 secs   av.it.:   8.9
      thresh= 1.890E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.858E-06

      iter #   3 total cpu time :     8.2 secs   av.it.:   8.7
      thresh= 1.363E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.212E-07

      iter #   4 total cpu time :     8.6 secs   av.it.:   8.4
      thresh= 3.482E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.478E-10

      iter #   5 total cpu time :     9.0 secs   av.it.:   8.6
      thresh= 2.116E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.247E-11

      iter #   6 total cpu time :     9.4 secs   av.it.:   8.4
      thresh= 3.532E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.797E-14

      iter #   7 total cpu time :     9.8 secs   av.it.:   8.5
      thresh= 1.340E-08 alpha_mix =  0.700 |ddv_scf|^2 =  2.169E-15

      iter #   8 total cpu time :    10.3 secs   av.it.:   8.8
      thresh= 4.657E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.077E-15

      iter #   9 total cpu time :    10.7 secs   av.it.:   8.5
      thresh= 3.282E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.092E-17

      iter #  10 total cpu time :    11.2 secs   av.it.:   8.7
      thresh= 3.305E-10 alpha_mix =  0.700 |ddv_scf|^2 =  3.862E-19

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 modes #   5  6

     Self-consistent Calculation

      iter #   1 total cpu time :    11.5 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.197E-06

      iter #   2 total cpu time :    12.1 secs   av.it.:   8.8
      thresh= 1.788E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.775E-07

      iter #   3 total cpu time :    12.5 secs   av.it.:   8.5
      thresh= 6.910E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.355E-08

      iter #   4 total cpu time :    13.0 secs   av.it.:   8.6
      thresh= 1.164E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.823E-10

      iter #   5 total cpu time :    13.4 secs   av.it.:   8.5
      thresh= 2.413E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.228E-12

      iter #   6 total cpu time :    13.8 secs   av.it.:   8.3
      thresh= 1.493E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.610E-15

      iter #   7 total cpu time :    14.2 secs   av.it.:   8.3
      thresh= 8.723E-09 alpha_mix =  0.700 |ddv_scf|^2 =  4.968E-16

      iter #   8 total cpu time :    14.6 secs   av.it.:   8.5
      thresh= 2.229E-09 alpha_mix =  0.700 |ddv_scf|^2 =  4.293E-17

      iter #   9 total cpu time :    15.1 secs   av.it.:   8.8
      thresh= 6.552E-10 alpha_mix =  0.700 |ddv_scf|^2 =  1.068E-17

      iter #  10 total cpu time :    15.5 secs   av.it.:   8.6
      thresh= 3.268E-10 alpha_mix =  0.700 |ddv_scf|^2 =  3.651E-19

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 modes #   7  8

     Self-consistent Calculation

      iter #   1 total cpu time :    15.9 secs   av.it.:   5.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.477E-06

      iter #   2 total cpu time :    16.3 secs   av.it.:   9.1
      thresh= 2.545E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.225E-06

      iter #   3 total cpu time :    16.8 secs   av.it.:   8.4
      thresh= 2.495E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.876E-08

      iter #   4 total cpu time :    17.2 secs   av.it.:   8.6
      thresh= 1.969E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.026E-10

      iter #   5 total cpu time :    17.7 secs   av.it.:   8.6
      thresh= 2.242E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.465E-11

      iter #   6 total cpu time :    18.1 secs   av.it.:   8.4
      thresh= 3.828E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.141E-14

      iter #   7 total cpu time :    18.5 secs   av.it.:   8.7
      thresh= 1.463E-08 alpha_mix =  0.700 |ddv_scf|^2 =  3.200E-15

      iter #   8 total cpu time :    19.0 secs   av.it.:   8.9
      thresh= 5.657E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.970E-15

      iter #   9 total cpu time :    19.4 secs   av.it.:   8.4
      thresh= 4.438E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.098E-17

      iter #  10 total cpu time :    19.9 secs   av.it.:   8.8
      thresh= 3.314E-10 alpha_mix =  0.700 |ddv_scf|^2 =  1.466E-18

      iter #  11 total cpu time :    20.3 secs   av.it.:   8.8
      thresh= 1.211E-10 alpha_mix =  0.700 |ddv_scf|^2 =  1.150E-20

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   9

     Self-consistent Calculation

      iter #   1 total cpu time :    20.6 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.621E-06

      iter #   2 total cpu time :    20.8 secs   av.it.:   8.1
      thresh= 1.903E-04 alpha_mix =  0.700 |ddv_scf|^2 =  4.244E-08

      iter #   3 total cpu time :    21.0 secs   av.it.:   7.8
      thresh= 2.060E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.458E-09

      iter #   4 total cpu time :    21.2 secs   av.it.:   6.5
      thresh= 3.818E-06 alpha_mix =  0.700 |ddv_scf|^2 =  9.844E-13

      iter #   5 total cpu time :    21.4 secs   av.it.:   6.6
      thresh= 9.922E-08 alpha_mix =  0.700 |ddv_scf|^2 =  8.631E-15

      iter #   6 total cpu time :    21.6 secs   av.it.:   6.2
      thresh= 9.290E-09 alpha_mix =  0.700 |ddv_scf|^2 =  6.627E-17

      iter #   7 total cpu time :    21.8 secs   av.it.:   7.0
      thresh= 8.141E-10 alpha_mix =  0.700 |ddv_scf|^2 =  1.041E-19

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  7 mode #  10

     Self-consistent Calculation

      iter #   1 total cpu time :    22.0 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.793E-06

      iter #   2 total cpu time :    22.2 secs   av.it.:   8.1
      thresh= 1.339E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.482E-08

      iter #   3 total cpu time :    22.4 secs   av.it.:   7.0
      thresh= 1.218E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.552E-09

      iter #   4 total cpu time :    22.6 secs   av.it.:   6.2
      thresh= 3.940E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.195E-12

      iter #   5 total cpu time :    22.8 secs   av.it.:   6.1
      thresh= 1.787E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.168E-15

      iter #   6 total cpu time :    23.0 secs   av.it.:   6.8
      thresh= 5.629E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.161E-16

      iter #   7 total cpu time :    23.2 secs   av.it.:   6.9
      thresh= 1.470E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.572E-19

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  8 modes #  11 12

     Self-consistent Calculation

      iter #   1 total cpu time :    23.6 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.643E-07

      iter #   2 total cpu time :    24.0 secs   av.it.:   8.8
      thresh= 8.150E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.223E-09

      iter #   3 total cpu time :    24.7 secs   av.it.:   8.6
      thresh= 7.227E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.981E-10

      iter #   4 total cpu time :    25.1 secs   av.it.:   7.3
      thresh= 1.727E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.636E-13

      iter #   5 total cpu time :    25.5 secs   av.it.:   8.3
      thresh= 5.134E-08 alpha_mix =  0.700 |ddv_scf|^2 =  5.411E-15

      iter #   6 total cpu time :    25.9 secs   av.it.:   7.8
      thresh= 7.356E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.215E-17

      iter #   7 total cpu time :    26.3 secs   av.it.:   8.1
      thresh= 4.706E-10 alpha_mix =  0.700 |ddv_scf|^2 =  6.675E-19

     End of self-consistent calculation

     Convergence has been achieved 
 
     Number of q in the star =    2
     List of q in the star:
          1   0.333333330   0.577350270   0.000000000
          2  -0.333333330  -0.577350270   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.333333330   0.577350270   0.000000000 ) 

 **************************************************************************
     freq (    1) =      14.906571 [THz] =     497.229703 [cm-1]
     freq (    2) =      15.257136 [THz] =     508.923280 [cm-1]
     freq (    3) =      15.279082 [THz] =     509.655325 [cm-1]
     freq (    4) =      15.279082 [THz] =     509.655325 [cm-1]
     freq (    5) =      30.274803 [THz] =    1009.858720 [cm-1]
     freq (    6) =      30.274803 [THz] =    1009.858720 [cm-1]
     freq (    7) =      36.689639 [THz] =    1223.834634 [cm-1]
     freq (    8) =      36.740092 [THz] =    1225.517545 [cm-1]
     freq (    9) =      36.850791 [THz] =    1229.210068 [cm-1]
     freq (   10) =      36.850791 [THz] =    1229.210068 [cm-1]
     freq (   11) =      39.415836 [THz] =    1314.770769 [cm-1]
     freq (   12) =      39.415836 [THz] =    1314.770769 [cm-1]
 **************************************************************************

     Mode symmetry, D_3h (-62m) point group:

     freq (  1 -  1) =        497.2  [cm-1]   --> A''2               
     freq (  2 -  2) =        508.9  [cm-1]   --> A''1               
     freq (  3 -  4) =        509.7  [cm-1]   --> E''                
     freq (  5 -  6) =       1009.9  [cm-1]   --> E'                 
     freq (  7 -  7) =       1223.8  [cm-1]   --> A'_1               
     freq (  8 -  8) =       1225.5  [cm-1]   --> A'_2               
     freq (  9 - 10) =       1229.2  [cm-1]   --> E'                 
     freq ( 11 - 12) =       1314.8  [cm-1]   --> E'                 
 
     init_run     :      0.03s CPU      0.05s WALL (       1 calls)
     electrons    :      0.29s CPU      0.34s WALL (       1 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       1 calls)
     potinit      :      0.02s CPU      0.02s WALL (       1 calls)
     hinit0       :      0.01s CPU      0.01s WALL (       1 calls)

     Called by electrons:
     c_bands      :      0.28s CPU      0.34s WALL (       1 calls)
     v_of_rho     :      1.65s CPU      1.69s WALL (       2 calls)
     newd         :      0.00s CPU      0.00s WALL (       2 calls)
     vdW_ffts     :      0.93s CPU      1.07s WALL (     220 calls)

     Called by c_bands:
     init_us_2    :      0.03s CPU      0.03s WALL (     860 calls)
     cegterg      :      0.24s CPU      0.28s WALL (      20 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :      8.26s CPU      9.68s WALL (   11176 calls)
     s_psi        :      0.59s CPU      0.69s WALL (   23202 calls)
     g_psi        :      0.00s CPU      0.00s WALL (     210 calls)
     cdiaghg      :      0.04s CPU      0.04s WALL (     230 calls)

     Called by h_psi:
     h_psi:calbec :      0.46s CPU      0.55s WALL (   11176 calls)
     vloc_psi     :      7.45s CPU      8.72s WALL (   11176 calls)
     add_vuspsi   :      0.28s CPU      0.34s WALL (   11176 calls)

     General routines
     calbec       :      0.99s CPU      1.17s WALL (   25042 calls)
     fft          :      1.34s CPU      1.56s WALL (   12143 calls)
     ffts         :      0.03s CPU      0.04s WALL (     462 calls)
     fftw         :      8.58s CPU     10.06s WALL (  185946 calls)
     interpolate  :      0.04s CPU      0.06s WALL (     234 calls)
     davcio       :      0.08s CPU      0.24s WALL (    5636 calls)
 
     Parallel routines
     fft_scatt_xy :      1.22s CPU      1.43s WALL (  198551 calls)
     fft_scatt_yz :      4.48s CPU      5.26s WALL (  198551 calls)
 
     PHONON       :     21.11s CPU     26.36s WALL

     INITIALIZATION: 
     phq_setup    :      0.01s CPU      0.02s WALL (       1 calls)
     phq_init     :      0.18s CPU      0.20s WALL (       1 calls)
 
     phq_init     :      0.18s CPU      0.20s WALL (       1 calls)
     init_vloc    :      0.00s CPU      0.00s WALL (       2 calls)
     init_us_1    :      0.01s CPU      0.02s WALL (       2 calls)
     newd         :      0.00s CPU      0.00s WALL (       2 calls)
     dvanqq       :      0.04s CPU      0.04s WALL (       1 calls)
     drho         :      0.12s CPU      0.13s WALL (       1 calls)
 
     DYNAMICAL MATRIX:
     dynmat0      :      0.01s CPU      0.02s WALL (       1 calls)
     phqscf       :     18.79s CPU     23.76s WALL (       1 calls)
     dynmatrix    :      0.00s CPU      0.01s WALL (       1 calls)
 
     phqscf       :     18.79s CPU     23.76s WALL (       1 calls)
     solve_linter :     18.74s CPU     23.52s WALL (       8 calls)
     drhodv       :      0.04s CPU      0.05s WALL (       8 calls)
 
     dynmat0      :      0.01s CPU      0.02s WALL (       1 calls)
     dynmat_us    :      0.01s CPU      0.01s WALL (       1 calls)
     d2ionq       :      0.00s CPU      0.00s WALL (       1 calls)
 
     dynmat_us    :      0.01s CPU      0.01s WALL (       1 calls)
     addusdynmat  :      0.00s CPU      0.00s WALL (       1 calls)
 
     phqscf       :     18.79s CPU     23.76s WALL (       1 calls)
     solve_linter :     18.74s CPU     23.52s WALL (       8 calls)
 
     solve_linter :     18.74s CPU     23.52s WALL (       8 calls)
     dvqpsi_us    :      0.14s CPU      0.16s WALL (     120 calls)
     ortho        :      0.13s CPU      0.15s WALL (    1100 calls)
     cgsolve      :      9.96s CPU     11.67s WALL (    1100 calls)
     incdrhoscf   :      1.06s CPU      1.25s WALL (    1100 calls)
     addusddens   :      0.30s CPU      0.31s WALL (      80 calls)
     vpsifft      :      0.87s CPU      1.02s WALL (     980 calls)
     dv_of_drho   :      4.54s CPU      5.45s WALL (     110 calls)
     mix_pot      :      0.15s CPU      0.81s WALL (      72 calls)
     psymdvscf    :      0.97s CPU      0.99s WALL (      72 calls)
     newdq        :      0.33s CPU      0.35s WALL (      72 calls)
     adddvscf     :      0.04s CPU      0.05s WALL (     980 calls)
     drhodvus     :      0.00s CPU      0.00s WALL (       8 calls)
 
     dvqpsi_us    :      0.14s CPU      0.16s WALL (     120 calls)
     dvqpsi_us_on :      0.03s CPU      0.03s WALL (     120 calls)
 
     cgsolve      :      9.96s CPU     11.67s WALL (    1100 calls)
     ch_psi       :      9.65s CPU     11.30s WALL (   10926 calls)
 
     ch_psi       :      9.65s CPU     11.30s WALL (   10926 calls)
     h_psi        :      8.26s CPU      9.68s WALL (   11176 calls)
     last         :      1.19s CPU      1.39s WALL (   10926 calls)
 
     h_psi        :      8.26s CPU      9.68s WALL (   11176 calls)
     add_vuspsi   :      0.28s CPU      0.34s WALL (   11176 calls)
 
     incdrhoscf   :      1.06s CPU      1.25s WALL (    1100 calls)
     addusdbec    :      0.07s CPU      0.08s WALL (    1220 calls)
 
     drhodvus     :      0.00s CPU      0.00s WALL (       8 calls)
 
      General routines
     calbec       :      0.99s CPU      1.17s WALL (   25042 calls)
     fft          :      1.34s CPU      1.56s WALL (   12143 calls)
     ffts         :      0.03s CPU      0.04s WALL (     462 calls)
     fftw         :      8.58s CPU     10.06s WALL (  185946 calls)
     davcio       :      0.08s CPU      0.24s WALL (    5636 calls)
     write_rec    :      0.08s CPU      1.26s WALL (      80 calls)
 
 
     PHONON       :     21.11s CPU     26.36s WALL

 
   This run was terminated on:  16:13:55  20Sep2019            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_36 = """
     Program PWSCF v.6.4.1 starts on 20Sep2019 at 16: 8:32 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
         "P. Giannozzi et al., J. Phys.:Condens. Matter 29 465901 (2017);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors

     MPI processes distributed on     1 nodes
     R & G space division:  proc/nbgrp/npool/nimage =       4
     Waiting for input...
     Reading input from standard input

     Current dimensions of program PWSCF are:
     Max number of different atomic species (ntypx) = 10
     Max number of k-points (npk) =  40000
     Max angular momentum in pseudopotentials (lmaxx) =  3
               file C.pbe-rrkjus.UPF: wavefunction(s)  2S 2P renormalized

     -------------------------------------------------
     Parameters for Dispersion (Grimme-D2) Correction:
     -------------------------------------------------
       atom      VdW radius       C_6     

        C          2.744         60.710

     Subspace diagonalization in iterative solution of the eigenvalue problem:
     a serial algorithm will be used

 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          65      46     15                 2276     1235     259
     Max          67      47     16                 2279     1237     264
     Sum         265     187     61                 9111     4943    1045
 


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.6463  a.u.
     unit-cell volume          =     224.3370 (a.u.)^3
     number of atoms/cell      =            4
     number of atomic types    =            1
     number of electrons       =        16.00
     number of Kohn-Sham states=            8
     kinetic-energy cutoff     =      30.0000  Ry
     charge density cutoff     =     180.0000  Ry
     convergence threshold     =      1.0E-12
     mixing beta               =       0.7000
     number of iterations used =            8  plain     mixing
     Exchange-correlation= SLA PW PBE PBE
                           (   1   4   3   4   0   0   0)

     celldm(1)=   4.646303  celldm(2)=   0.000000  celldm(3)=   2.582543
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   2.582543 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.387215 )  


     PseudoPot. # 1 for C  read from file:
     /benchmarks/pseudo/C.pbe-rrkjus.UPF
     MD5 check sum: c9ac5a99bc85b198593446162950cd17
     Pseudo is Ultrasoft, Zval =  4.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  627 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        C              4.00    12.00000     C ( 1.00)

     24 Sym. Ops., with inversion, found (12 have fractional translation)



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           C   tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           C   tau(   2) = (  -0.0000000   0.5773503   0.0000000  )
         3           C   tau(   3) = (   0.0000000   0.0000000   1.2912714  )
         4           C   tau(   4) = (   0.5000000   0.2886751   1.2912714  )

     number of k points=     6
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.1250000   0.2165064   0.0968038), wk =   0.2500000
        k(    2) = (   0.1250000   0.5051815   0.0968038), wk =   0.5000000
        k(    3) = (   0.1250000  -0.3608439   0.0968038), wk =   0.5000000
        k(    4) = (   0.1250000  -0.0721688   0.0968038), wk =   0.2500000
        k(    5) = (   0.3750000   0.6495191   0.0968038), wk =   0.2500000
        k(    6) = (   0.3750000  -0.2165064   0.0968038), wk =   0.2500000

     Dense  grid:     9111 G-vectors     FFT dimensions: (  20,  20,  54)

     Smooth grid:     4943 G-vectors     FFT dimensions: (  18,  18,  48)

     Estimated max dynamical RAM per process >       3.39 MB

     Estimated total dynamical RAM >      13.58 MB

     Initial potential from superposition of free atoms

     starting charge   15.99979, renormalised to   16.00000
     Starting wfcs are   16 randomized atomic wfcs

     total cpu time spent up to now is        0.2 secs

     Self-consistent Calculation

     iteration #  1     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-02,  avg # of iterations =  2.0

     total cpu time spent up to now is        0.4 secs

     total energy              =     -45.55818426 Ry
     Harris-Foulkes estimate   =     -45.78070131 Ry
     estimated scf accuracy    <       0.38648672 Ry

     iteration #  2     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.42E-03,  avg # of iterations =  2.0

     total cpu time spent up to now is        0.4 secs

     total energy              =     -45.61210746 Ry
     Harris-Foulkes estimate   =     -45.61204330 Ry
     estimated scf accuracy    <       0.00421417 Ry

     iteration #  3     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.63E-05,  avg # of iterations =  2.5

     total cpu time spent up to now is        0.4 secs

     total energy              =     -45.61275456 Ry
     Harris-Foulkes estimate   =     -45.61270647 Ry
     estimated scf accuracy    <       0.00026869 Ry

     iteration #  4     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.68E-06,  avg # of iterations =  1.8

     total cpu time spent up to now is        0.4 secs

     total energy              =     -45.61278001 Ry
     Harris-Foulkes estimate   =     -45.61278036 Ry
     estimated scf accuracy    <       0.00000370 Ry

     iteration #  5     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.31E-08,  avg # of iterations =  3.3

     total cpu time spent up to now is        0.5 secs

     total energy              =     -45.61278275 Ry
     Harris-Foulkes estimate   =     -45.61278260 Ry
     estimated scf accuracy    <       0.00000010 Ry

     iteration #  6     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  6.36E-10,  avg # of iterations =  3.0

     total cpu time spent up to now is        0.5 secs

     total energy              =     -45.61278278 Ry
     Harris-Foulkes estimate   =     -45.61278279 Ry
     estimated scf accuracy    <       0.00000003 Ry

     iteration #  7     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.79E-10,  avg # of iterations =  2.8

     total cpu time spent up to now is        0.5 secs

     total energy              =     -45.61278279 Ry
     Harris-Foulkes estimate   =     -45.61278279 Ry
     estimated scf accuracy    <          6.0E-09 Ry

     iteration #  8     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  3.72E-11,  avg # of iterations =  2.3

     total cpu time spent up to now is        0.6 secs

     total energy              =     -45.61278279 Ry
     Harris-Foulkes estimate   =     -45.61278279 Ry
     estimated scf accuracy    <          5.3E-10 Ry

     iteration #  9     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  3.34E-12,  avg # of iterations =  3.2

     total cpu time spent up to now is        0.6 secs

     total energy              =     -45.61278279 Ry
     Harris-Foulkes estimate   =     -45.61278279 Ry
     estimated scf accuracy    <          1.0E-10 Ry

     iteration # 10     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  6.29E-13,  avg # of iterations =  2.5

     total cpu time spent up to now is        0.6 secs

     total energy              =     -45.61278279 Ry
     Harris-Foulkes estimate   =     -45.61278279 Ry
     estimated scf accuracy    <          1.7E-12 Ry

     iteration # 11     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-13,  avg # of iterations =  1.0

     total cpu time spent up to now is        0.7 secs

     total energy              =     -45.61278279 Ry
     Harris-Foulkes estimate   =     -45.61278279 Ry
     estimated scf accuracy    <          1.6E-12 Ry

     iteration # 12     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-13,  avg # of iterations =  1.0

     total cpu time spent up to now is        0.7 secs

     End of self-consistent calculation

          k = 0.1250 0.2165 0.0968 (   619 PWs)   bands (ev):

   -10.8945 -10.5952   0.5344   1.4501   1.5006   2.2444   2.5206   2.5563

          k = 0.1250 0.5052 0.0968 (   628 PWs)   bands (ev):

    -7.3330  -7.1046  -4.4133  -4.2383   0.2015   0.2933   4.5392   5.6590

          k = 0.1250-0.3608 0.0968 (   627 PWs)   bands (ev):

    -9.4327  -9.1599  -1.3776  -1.2545   0.9774   1.0547   2.2612   3.8298

          k = 0.1250-0.0722 0.0968 (   606 PWs)   bands (ev):

   -11.6367 -11.3241  -0.3581   1.4069   3.2479   3.3034   3.9269   3.9682

          k = 0.3750 0.6495 0.0968 (   630 PWs)   bands (ev):

    -5.6886  -5.5759  -4.7821  -4.7335  -2.1298  -2.0155   6.2547   6.8861

          k = 0.3750-0.2165 0.0968 (   623 PWs)   bands (ev):

    -8.7147  -8.4553  -3.0976  -2.9327   1.5908   1.6535   3.0578   4.5017

     highest occupied level (ev):     6.8861

!    total energy              =     -45.61278279 Ry
     Harris-Foulkes estimate   =     -45.61278279 Ry
     estimated scf accuracy    <          4.5E-14 Ry

     The total energy is the sum of the following terms:

     one-electron contribution =      -7.06559982 Ry
     hartree contribution      =      11.83976083 Ry
     xc contribution           =     -14.07677349 Ry
     ewald contribution        =     -36.26917740 Ry
     Dispersion Correction     =      -0.04099291 Ry

     convergence has been achieved in  12 iterations

     Writing output data file /benchmarks/tempdir/graphite.save/
 
     init_run     :      0.05s CPU      0.06s WALL (       1 calls)
     electrons    :      0.43s CPU      0.48s WALL (       1 calls)

     Called by init_run:
     wfcinit      :      0.02s CPU      0.02s WALL (       1 calls)
     potinit      :      0.01s CPU      0.01s WALL (       1 calls)
     hinit0       :      0.02s CPU      0.02s WALL (       1 calls)

     Called by electrons:
     c_bands      :      0.22s CPU      0.24s WALL (      12 calls)
     sum_band     :      0.05s CPU      0.06s WALL (      12 calls)
     v_of_rho     :      0.03s CPU      0.03s WALL (      13 calls)
     newd         :      0.01s CPU      0.02s WALL (      13 calls)
     mix_rho      :      0.01s CPU      0.01s WALL (      12 calls)

     Called by c_bands:
     init_us_2    :      0.00s CPU      0.01s WALL (     150 calls)
     cegterg      :      0.20s CPU      0.22s WALL (      72 calls)

     Called by sum_band:
     sum_band:bec :      0.00s CPU      0.00s WALL (      72 calls)
     addusdens    :      0.01s CPU      0.01s WALL (      12 calls)

     Called by *egterg:
     h_psi        :      0.16s CPU      0.17s WALL (     243 calls)
     s_psi        :      0.01s CPU      0.01s WALL (     243 calls)
     g_psi        :      0.00s CPU      0.00s WALL (     165 calls)
     cdiaghg      :      0.03s CPU      0.03s WALL (     237 calls)

     Called by h_psi:
     h_psi:calbec :      0.01s CPU      0.01s WALL (     243 calls)
     vloc_psi     :      0.14s CPU      0.15s WALL (     243 calls)
     add_vuspsi   :      0.01s CPU      0.01s WALL (     243 calls)

     General routines
     calbec       :      0.01s CPU      0.02s WALL (     315 calls)
     fft          :      0.02s CPU      0.02s WALL (     166 calls)
     ffts         :      0.00s CPU      0.00s WALL (      25 calls)
     fftw         :      0.15s CPU      0.16s WALL (    3858 calls)
     interpolate  :      0.00s CPU      0.00s WALL (      13 calls)
 
     Parallel routines
     fft_scatt_xy :      0.02s CPU      0.03s WALL (    4049 calls)
     fft_scatt_yz :      0.06s CPU      0.07s WALL (    4049 calls)
 
     PWSCF        :      0.62s CPU      0.76s WALL

 
   This run was terminated on:  16: 8:33  20Sep2019            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_41 = """
     Program PWSCF v.6.4.1 starts on 20Sep2019 at 16: 9: 7 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
         "P. Giannozzi et al., J. Phys.:Condens. Matter 29 465901 (2017);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors

     MPI processes distributed on     1 nodes
     R & G space division:  proc/nbgrp/npool/nimage =       4
     Waiting for input...
     Reading input from standard input

     Current dimensions of program PWSCF are:
     Max number of different atomic species (ntypx) = 10
     Max number of k-points (npk) =  40000
     Max angular momentum in pseudopotentials (lmaxx) =  3

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation= VDW-DF
                           (   1   4   4   0   1   0   0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want

               file C.pbe-rrkjus.UPF: wavefunction(s)  2S 2P renormalized

     Subspace diagonalization in iterative solution of the eigenvalue problem:
     a serial algorithm will be used

 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min          65      46     15                 2276     1235     259
     Max          67      47     16                 2279     1237     264
     Sum         265     187     61                 9111     4943    1045
 


     bravais-lattice index     =            4
     lattice parameter (alat)  =       4.6463  a.u.
     unit-cell volume          =     224.3370 (a.u.)^3
     number of atoms/cell      =            4
     number of atomic types    =            1
     number of electrons       =        16.00
     number of Kohn-Sham states=            8
     kinetic-energy cutoff     =      30.0000  Ry
     charge density cutoff     =     180.0000  Ry
     convergence threshold     =      1.0E-12
     mixing beta               =       0.7000
     number of iterations used =            8  plain     mixing
     Exchange-correlation= VDW-DF
                           (   1   4   4   0   1   0   0)

     celldm(1)=   4.646303  celldm(2)=   0.000000  celldm(3)=   2.582543
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (   1.000000   0.000000   0.000000 )  
               a(2) = (  -0.500000   0.866025   0.000000 )  
               a(3) = (   0.000000   0.000000   2.582543 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = (  1.000000  0.577350  0.000000 )  
               b(2) = (  0.000000  1.154701  0.000000 )  
               b(3) = (  0.000000  0.000000  0.387215 )  


     PseudoPot. # 1 for C  read from file:
     /benchmarks/pseudo/C.pbe-rrkjus.UPF
     MD5 check sum: c9ac5a99bc85b198593446162950cd17
     Pseudo is Ultrasoft, Zval =  4.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  627 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        C              4.00    12.00000     C ( 1.00)

     24 Sym. Ops., with inversion, found (12 have fractional translation)



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           C   tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           C   tau(   2) = (  -0.0000000   0.5773503   0.0000000  )
         3           C   tau(   3) = (   0.0000000   0.0000000   1.2912714  )
         4           C   tau(   4) = (   0.5000000   0.2886751   1.2912714  )

     number of k points=     6
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.1250000   0.2165064   0.0968038), wk =   0.2500000
        k(    2) = (   0.1250000   0.5051815   0.0968038), wk =   0.5000000
        k(    3) = (   0.1250000  -0.3608439   0.0968038), wk =   0.5000000
        k(    4) = (   0.1250000  -0.0721688   0.0968038), wk =   0.2500000
        k(    5) = (   0.3750000   0.6495191   0.0968038), wk =   0.2500000
        k(    6) = (   0.3750000  -0.2165064   0.0968038), wk =   0.2500000

     Dense  grid:     9111 G-vectors     FFT dimensions: (  20,  20,  54)

     Smooth grid:     4943 G-vectors     FFT dimensions: (  18,  18,  48)

     Estimated max dynamical RAM per process >       3.39 MB

     Estimated total dynamical RAM >      13.58 MB

     Initial potential from superposition of free atoms

     starting charge   15.99979, renormalised to   16.00000


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %                                                                      %
     % You are using vdW-DF, which was implemented by the Thonhauser group. %
     % Please cite the following two papers that made this development      %
     % possible and the two reviews that describe the various versions:     %
     %                                                                      %
     %   T. Thonhauser et al., PRL 115, 136402 (2015).                      %
     %   T. Thonhauser et al., PRB 76, 125112 (2007).                       %
     %   K. Berland et al., Rep. Prog. Phys. 78, 066501 (2015).             %
     %   D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009). %
     %                                                                      %
     %                                                                      %
     % If you are calculating the stress with vdW-DF, please also cite:     %
     %                                                                      %
     %   R. Sabatini et al., J. Phys.: Condens. Matter 24, 424209 (2012).   %
     %                                                                      %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     Starting wfcs are   16 randomized atomic wfcs

     total cpu time spent up to now is       25.4 secs

     Self-consistent Calculation

     iteration #  1     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-02,  avg # of iterations =  2.0

     total cpu time spent up to now is       25.4 secs

     total energy              =     -45.81537064 Ry
     Harris-Foulkes estimate   =     -46.05259572 Ry
     estimated scf accuracy    <       0.42906803 Ry

     iteration #  2     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.68E-03,  avg # of iterations =  2.0

     total cpu time spent up to now is       25.5 secs

     total energy              =     -45.87729814 Ry
     Harris-Foulkes estimate   =     -45.87595856 Ry
     estimated scf accuracy    <       0.00525176 Ry

     iteration #  3     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  3.28E-05,  avg # of iterations =  2.5

     total cpu time spent up to now is       25.5 secs

     total energy              =     -45.87812971 Ry
     Harris-Foulkes estimate   =     -45.87794212 Ry
     estimated scf accuracy    <       0.00039868 Ry

     iteration #  4     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.49E-06,  avg # of iterations =  1.5

     total cpu time spent up to now is       25.5 secs

     total energy              =     -45.87816598 Ry
     Harris-Foulkes estimate   =     -45.87816437 Ry
     estimated scf accuracy    <       0.00000423 Ry

     iteration #  5     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.64E-08,  avg # of iterations =  3.5

     total cpu time spent up to now is       25.6 secs

     total energy              =     -45.87816949 Ry
     Harris-Foulkes estimate   =     -45.87816934 Ry
     estimated scf accuracy    <       0.00000034 Ry

     iteration #  6     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  2.14E-09,  avg # of iterations =  3.3

     total cpu time spent up to now is       25.6 secs

     total energy              =     -45.87816953 Ry
     Harris-Foulkes estimate   =     -45.87816960 Ry
     estimated scf accuracy    <       0.00000013 Ry

     iteration #  7     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  8.30E-10,  avg # of iterations =  2.7

     total cpu time spent up to now is       25.7 secs

     total energy              =     -45.87816956 Ry
     Harris-Foulkes estimate   =     -45.87816957 Ry
     estimated scf accuracy    <       0.00000003 Ry

     iteration #  8     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.76E-10,  avg # of iterations =  2.7

     total cpu time spent up to now is       25.7 secs

     total energy              =     -45.87816957 Ry
     Harris-Foulkes estimate   =     -45.87816957 Ry
     estimated scf accuracy    <          1.8E-10 Ry

     iteration #  9     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.10E-12,  avg # of iterations =  3.8

     total cpu time spent up to now is       25.8 secs

     total energy              =     -45.87816957 Ry
     Harris-Foulkes estimate   =     -45.87816957 Ry
     estimated scf accuracy    <          1.1E-10 Ry

     iteration # 10     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  6.67E-13,  avg # of iterations =  2.7

     total cpu time spent up to now is       25.8 secs

     total energy              =     -45.87816957 Ry
     Harris-Foulkes estimate   =     -45.87816957 Ry
     estimated scf accuracy    <          3.1E-12 Ry

     iteration # 11     ecut=    30.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-13,  avg # of iterations =  1.8

     total cpu time spent up to now is       25.9 secs

     End of self-consistent calculation

          k = 0.1250 0.2165 0.0968 (   619 PWs)   bands (ev):

   -10.9292 -10.6363   0.5597   1.2917   1.3409   2.2425   2.3780   2.4123

          k = 0.1250 0.5052 0.0968 (   628 PWs)   bands (ev):

    -7.4066  -7.1849  -4.4847  -4.3166   0.0918   0.1800   4.5348   5.6329

          k = 0.1250-0.3608 0.0968 (   627 PWs)   bands (ev):

    -9.4801  -9.2135  -1.4960  -1.3776   0.8616   0.9353   2.2769   3.8186

          k = 0.1250-0.0722 0.0968 (   606 PWs)   bands (ev):

   -11.6656 -11.3594  -0.3285   1.4093   3.0718   3.1253   3.7433   3.7832

          k = 0.3750 0.6495 0.0968 (   630 PWs)   bands (ev):

    -5.8014  -5.6911  -4.8560  -4.8135  -2.1716  -2.0649   6.2356   6.8601

          k = 0.3750-0.2165 0.0968 (   623 PWs)   bands (ev):

    -8.7696  -8.5164  -3.1851  -3.0257   1.4468   1.5077   3.0669   4.4844

     highest occupied level (ev):     6.8601

!    total energy              =     -45.87816957 Ry
     Harris-Foulkes estimate   =     -45.87816957 Ry
     estimated scf accuracy    <          8.5E-13 Ry

     The total energy is the sum of the following terms:

     one-electron contribution =      -7.16846936 Ry
     hartree contribution      =      11.99235762 Ry
     xc contribution           =     -14.43288043 Ry
     ewald contribution        =     -36.26917740 Ry

     convergence has been achieved in  11 iterations

     Writing output data file /benchmarks/tempdir/graphite.save/
 
     init_run     :     25.23s CPU     25.25s WALL (       1 calls)
     electrons    :      0.42s CPU      0.49s WALL (       1 calls)

     Called by init_run:
     wfcinit      :      0.01s CPU      0.02s WALL (       1 calls)
     potinit      :     25.20s CPU     25.21s WALL (       1 calls)
     hinit0       :      0.01s CPU      0.01s WALL (       1 calls)

     Called by electrons:
     c_bands      :      0.21s CPU      0.24s WALL (      11 calls)
     sum_band     :      0.05s CPU      0.06s WALL (      11 calls)
     v_of_rho     :     25.33s CPU     25.37s WALL (      12 calls)
     newd         :      0.01s CPU      0.02s WALL (      12 calls)
     mix_rho      :      0.01s CPU      0.01s WALL (      11 calls)
     vdW_kernel   :     25.18s CPU     25.19s WALL (       1 calls)

     Called by c_bands:
     init_us_2    :      0.00s CPU      0.01s WALL (     138 calls)
     cegterg      :      0.20s CPU      0.22s WALL (      66 calls)

     Called by sum_band:
     sum_band:bec :      0.00s CPU      0.00s WALL (      66 calls)
     addusdens    :      0.01s CPU      0.01s WALL (      11 calls)

     Called by *egterg:
     h_psi        :      0.15s CPU      0.17s WALL (     243 calls)
     s_psi        :      0.01s CPU      0.01s WALL (     243 calls)
     g_psi        :      0.00s CPU      0.00s WALL (     171 calls)
     cdiaghg      :      0.03s CPU      0.03s WALL (     237 calls)

     Called by h_psi:
     h_psi:calbec :      0.01s CPU      0.01s WALL (     243 calls)
     vloc_psi     :      0.13s CPU      0.14s WALL (     243 calls)
     add_vuspsi   :      0.01s CPU      0.01s WALL (     243 calls)

     General routines
     calbec       :      0.01s CPU      0.01s WALL (     309 calls)
     fft          :      0.08s CPU      0.08s WALL (     753 calls)
     ffts         :      0.00s CPU      0.00s WALL (      23 calls)
     fftw         :      0.14s CPU      0.15s WALL (    3744 calls)
     interpolate  :      0.00s CPU      0.00s WALL (      12 calls)
 
     Parallel routines
     fft_scatt_xy :      0.03s CPU      0.03s WALL (    4520 calls)
     fft_scatt_yz :      0.08s CPU      0.09s WALL (    4520 calls)
 
     PWSCF        :     25.76s CPU     25.93s WALL

 
   This run was terminated on:  16: 9:33  20Sep2019            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_46 = """
     Program PHONON v.6.0 (svn rev. 13188M) starts on  7Dec2016 at 13:43:37 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     R & G space division:  proc/nbgrp/npool/nimage =       4

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/carbon.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         115      40     15                 1641      351      82
     Max         117      41     16                 1642      354      83
     Sum         463     163     61                 6567     1411     331


     Check: negative/imaginary core charge=   -0.000005    0.000000

     Calculation of q =    1.0000000   0.0000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         115      40     21                 1641      351     132
     Max         117      41     22                 1642      354     135
     Sum         463     163     85                 6567     1411     531


     Title: 
     phonons of C at X                                                          


     bravais-lattice index     =            2
     lattice parameter (alat)  =       6.6800  a.u.
     unit-cell volume          =      74.5194 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            1
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      27.0000  Ry
     charge density cutoff     =     300.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=   6.680000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for  C read from file:
     ./C.pz-kjpaw.UPF
     MD5 check sum: 414e6e825ae75add557e798061b49a04
     Pseudo is Projector augmented-wave + core cor, Zval =  4.0
     Generated using "atomic" code by A. Dal Corso (espresso distribution)
     Shape of augmentation charge: BESSEL
     Using radial grid of 1073 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        C              4.00    12.01070      C( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           C   tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           C   tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    40
                       cart. coord. in units 2pi/alat
        k(    1) = (  -0.1250000   0.1250000   0.1250000), wk =   0.0625000
        k(    2) = (   0.8750000   0.1250000   0.1250000), wk =   0.0000000
        k(    3) = (  -0.3750000   0.3750000  -0.1250000), wk =   0.1250000
        k(    4) = (   0.6250000   0.3750000  -0.1250000), wk =   0.0000000
        k(    5) = (   0.3750000  -0.3750000   0.6250000), wk =   0.1250000
        k(    6) = (   1.3750000  -0.3750000   0.6250000), wk =   0.0000000
        k(    7) = (   0.1250000  -0.1250000   0.3750000), wk =   0.1250000
        k(    8) = (   1.1250000  -0.1250000   0.3750000), wk =   0.0000000
        k(    9) = (  -0.1250000   0.6250000   0.1250000), wk =   0.1250000
        k(   10) = (   0.8750000   0.6250000   0.1250000), wk =   0.0000000
        k(   11) = (   0.6250000  -0.1250000   0.8750000), wk =   0.1250000
        k(   12) = (   1.6250000  -0.1250000   0.8750000), wk =   0.0000000
        k(   13) = (   0.3750000   0.1250000   0.6250000), wk =   0.1250000
        k(   14) = (   1.3750000   0.1250000   0.6250000), wk =   0.0000000
        k(   15) = (  -0.1250000  -0.8750000   0.1250000), wk =   0.1250000
        k(   16) = (   0.8750000  -0.8750000   0.1250000), wk =   0.0000000
        k(   17) = (  -0.3750000   0.3750000   0.3750000), wk =   0.0625000
        k(   18) = (   0.6250000   0.3750000   0.3750000), wk =   0.0000000
        k(   19) = (   0.3750000  -0.3750000   1.1250000), wk =   0.1250000
        k(   20) = (   1.3750000  -0.3750000   1.1250000), wk =   0.0000000
        k(   21) = (  -0.1250000  -0.3750000   0.3750000), wk =   0.0625000
        k(   22) = (   0.8750000  -0.3750000   0.3750000), wk =   0.0000000
        k(   23) = (   0.6250000   0.3750000  -0.3750000), wk =   0.0625000
        k(   24) = (   1.6250000   0.3750000  -0.3750000), wk =   0.0000000
        k(   25) = (   0.3750000   0.1250000  -0.1250000), wk =   0.0625000
        k(   26) = (   1.3750000   0.1250000  -0.1250000), wk =   0.0000000
        k(   27) = (   0.6250000   0.1250000  -0.1250000), wk =   0.0625000
        k(   28) = (   1.6250000   0.1250000  -0.1250000), wk =   0.0000000
        k(   29) = (  -0.1250000   0.8750000   0.6250000), wk =   0.1250000
        k(   30) = (   0.8750000   0.8750000   0.6250000), wk =   0.0000000
        k(   31) = (   0.8750000   0.6250000  -0.1250000), wk =   0.1250000
        k(   32) = (   1.8750000   0.6250000  -0.1250000), wk =   0.0000000
        k(   33) = (   0.1250000   0.6250000   0.3750000), wk =   0.1250000
        k(   34) = (   1.1250000   0.6250000   0.3750000), wk =   0.0000000
        k(   35) = (   0.6250000   0.3750000   0.1250000), wk =   0.1250000
        k(   36) = (   1.6250000   0.3750000   0.1250000), wk =   0.0000000
        k(   37) = (  -0.8750000   0.1250000  -0.1250000), wk =   0.0625000
        k(   38) = (   0.1250000   0.1250000  -0.1250000), wk =   0.0000000
        k(   39) = (   1.1250000   0.3750000  -0.3750000), wk =   0.0625000
        k(   40) = (   2.1250000   0.3750000  -0.3750000), wk =   0.0000000

     Dense  grid:     6567 G-vectors     FFT dimensions: (  32,  32,  32)

     Smooth grid:     1411 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.77Mb

     Estimated total allocated dynamical RAM >       3.06Mb

     Check: negative/imaginary core charge=   -0.000005    0.000000

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/carbon.save/charge-density.dat

     Starting wfc are    8 atomic wfcs
     Checking if some PAW data can be deallocated... 

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.7

     total cpu time spent up to now is        0.4 secs

     End of band structure calculation

          k =-0.1250 0.1250 0.1250 (   172 PWs)   bands (ev):

    -7.6538  11.2130  13.0027  13.0027

          k = 0.8750 0.1250 0.1250 (   178 PWs)   bands (ev):

    -0.9281   2.9486   6.1002   7.5711

          k =-0.3750 0.3750-0.1250 (   173 PWs)   bands (ev):

    -5.3443   5.7928   9.2182  11.8188

          k = 0.6250 0.3750-0.1250 (   175 PWs)   bands (ev):

    -3.1174   3.6701   7.2810   9.4825

          k = 0.3750-0.3750 0.6250 (   172 PWs)   bands (ev):

    -2.1992   0.8512   9.4650  10.0574

          k = 1.3750-0.3750 0.6250 (   172 PWs)   bands (ev):

    -2.1992   0.8512   9.4650  10.0574

          k = 0.1250-0.1250 0.3750 (   172 PWs)   bands (ev):

    -6.4897   8.5250  10.7067  11.6521

          k = 1.1250-0.1250 0.3750 (   180 PWs)   bands (ev):

    -0.1270   2.3914   5.0210   7.1187

          k =-0.1250 0.6250 0.1250 (   180 PWs)   bands (ev):

    -4.2071   6.0838   8.3417   8.6814

          k = 0.8750 0.6250 0.1250 (   180 PWs)   bands (ev):

    -0.1270   2.3914   5.0210   7.1187

          k = 0.6250-0.1250 0.8750 (   180 PWs)   bands (ev):

    -0.1270   2.3914   5.0210   7.1187

          k = 1.6250-0.1250 0.8750 (   180 PWs)   bands (ev):

    -0.1270   2.3914   5.0210   7.1187

          k = 0.3750 0.1250 0.6250 (   175 PWs)   bands (ev):

    -3.1174   3.6701   7.2810   9.4825

          k = 1.3750 0.1250 0.6250 (   178 PWs)   bands (ev):

    -1.1767   1.8657   5.7621   9.3468

          k =-0.1250-0.8750 0.1250 (   178 PWs)   bands (ev):

    -0.9281   2.9486   6.1002   7.5711

          k = 0.8750-0.8750 0.1250 (   178 PWs)   bands (ev):

    -0.9281   2.9486   6.1002   7.5711

          k =-0.3750 0.3750 0.3750 (   177 PWs)   bands (ev):

    -4.2308   2.7971  11.0868  11.0868

          k = 0.6250 0.3750 0.3750 (   172 PWs)   bands (ev):

    -2.1992   0.8512   9.4650  10.0574

          k = 0.3750-0.3750 1.1250 (   178 PWs)   bands (ev):

    -1.1767   1.8657   5.7621   9.3468

          k = 1.3750-0.3750 1.1250 (   175 PWs)   bands (ev):

    -3.1174   3.6701   7.2810   9.4825

          k =-0.1250-0.3750 0.3750 (   173 PWs)   bands (ev):

    -5.3443   5.7928   9.2182  11.8188

          k = 0.8750-0.3750 0.3750 (   178 PWs)   bands (ev):

    -1.1767   1.8657   5.7621   9.3468

          k = 0.6250 0.3750-0.3750 (   172 PWs)   bands (ev):

    -2.1992   0.8512   9.4650  10.0574

          k = 1.6250 0.3750-0.3750 (   177 PWs)   bands (ev):

    -4.2308   2.7971  11.0868  11.0868

          k = 0.3750 0.1250-0.1250 (   172 PWs)   bands (ev):

    -6.4897   8.5250  10.7067  11.6521

          k = 1.3750 0.1250-0.1250 (   180 PWs)   bands (ev):

    -4.2071   6.0838   8.3417   8.6814

          k = 0.6250 0.1250-0.1250 (   180 PWs)   bands (ev):

    -4.2071   6.0838   8.3417   8.6814

          k = 1.6250 0.1250-0.1250 (   172 PWs)   bands (ev):

    -6.4897   8.5250  10.7067  11.6521

          k =-0.1250 0.8750 0.6250 (   180 PWs)   bands (ev):

    -0.1270   2.3914   5.0210   7.1187

          k = 0.8750 0.8750 0.6250 (   172 PWs)   bands (ev):

    -6.4897   8.5250  10.7067  11.6521

          k = 0.8750 0.6250-0.1250 (   180 PWs)   bands (ev):

    -0.1270   2.3914   5.0210   7.1187

          k = 1.8750 0.6250-0.1250 (   180 PWs)   bands (ev):

    -4.2071   6.0838   8.3417   8.6814

          k = 0.1250 0.6250 0.3750 (   175 PWs)   bands (ev):

    -3.1174   3.6701   7.2810   9.4825

          k = 1.1250 0.6250 0.3750 (   175 PWs)   bands (ev):

    -3.1174   3.6701   7.2810   9.4825

          k = 0.6250 0.3750 0.1250 (   175 PWs)   bands (ev):

    -3.1174   3.6701   7.2810   9.4825

          k = 1.6250 0.3750 0.1250 (   173 PWs)   bands (ev):

    -5.3443   5.7928   9.2182  11.8188

          k =-0.8750 0.1250-0.1250 (   178 PWs)   bands (ev):

    -0.9281   2.9486   6.1002   7.5711

          k = 0.1250 0.1250-0.1250 (   172 PWs)   bands (ev):

    -7.6538  11.2130  13.0027  13.0027

          k = 1.1250 0.3750-0.3750 (   178 PWs)   bands (ev):

    -1.1767   1.8657   5.7621   9.3468

          k = 2.1250 0.3750-0.3750 (   173 PWs)   bands (ev):

    -5.3443   5.7928   9.2182  11.8188

     highest occupied level (ev):    13.0027

     Writing output data file carbon.save

     phonons of C at X                                                          

     bravais-lattice index     =            2
     lattice parameter (alat)  =       6.6800  a.u.
     unit-cell volume          =      74.5194 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            1
     kinetic-energy cut-off    =      27.0000  Ry
     charge density cut-off    =     300.0000  Ry
     convergence threshold     =      1.0E-14
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=    6.68000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     C   12.0107   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     C   12.0107   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   1.0000000   0.0000000   0.0000000 )

     17 Sym.Ops. (with q -> -q+G )


     G cutoff =  339.0896  (   1642 G-vectors)     FFT grid: ( 32, 32, 32)
     G cutoff =  122.0722  (    354 G-vectors)  smooth grid: ( 15, 15, 15)
     number of k points=    40
                       cart. coord. in units 2pi/alat
        k(    1) = (  -0.1250000   0.1250000   0.1250000), wk =   0.0625000
        k(    2) = (   0.8750000   0.1250000   0.1250000), wk =   0.0000000
        k(    3) = (  -0.3750000   0.3750000  -0.1250000), wk =   0.1250000
        k(    4) = (   0.6250000   0.3750000  -0.1250000), wk =   0.0000000
        k(    5) = (   0.3750000  -0.3750000   0.6250000), wk =   0.1250000
        k(    6) = (   1.3750000  -0.3750000   0.6250000), wk =   0.0000000
        k(    7) = (   0.1250000  -0.1250000   0.3750000), wk =   0.1250000
        k(    8) = (   1.1250000  -0.1250000   0.3750000), wk =   0.0000000
        k(    9) = (  -0.1250000   0.6250000   0.1250000), wk =   0.1250000
        k(   10) = (   0.8750000   0.6250000   0.1250000), wk =   0.0000000
        k(   11) = (   0.6250000  -0.1250000   0.8750000), wk =   0.1250000
        k(   12) = (   1.6250000  -0.1250000   0.8750000), wk =   0.0000000
        k(   13) = (   0.3750000   0.1250000   0.6250000), wk =   0.1250000
        k(   14) = (   1.3750000   0.1250000   0.6250000), wk =   0.0000000
        k(   15) = (  -0.1250000  -0.8750000   0.1250000), wk =   0.1250000
        k(   16) = (   0.8750000  -0.8750000   0.1250000), wk =   0.0000000
        k(   17) = (  -0.3750000   0.3750000   0.3750000), wk =   0.0625000
        k(   18) = (   0.6250000   0.3750000   0.3750000), wk =   0.0000000
        k(   19) = (   0.3750000  -0.3750000   1.1250000), wk =   0.1250000
        k(   20) = (   1.3750000  -0.3750000   1.1250000), wk =   0.0000000
        k(   21) = (  -0.1250000  -0.3750000   0.3750000), wk =   0.0625000
        k(   22) = (   0.8750000  -0.3750000   0.3750000), wk =   0.0000000
        k(   23) = (   0.6250000   0.3750000  -0.3750000), wk =   0.0625000
        k(   24) = (   1.6250000   0.3750000  -0.3750000), wk =   0.0000000
        k(   25) = (   0.3750000   0.1250000  -0.1250000), wk =   0.0625000
        k(   26) = (   1.3750000   0.1250000  -0.1250000), wk =   0.0000000
        k(   27) = (   0.6250000   0.1250000  -0.1250000), wk =   0.0625000
        k(   28) = (   1.6250000   0.1250000  -0.1250000), wk =   0.0000000
        k(   29) = (  -0.1250000   0.8750000   0.6250000), wk =   0.1250000
        k(   30) = (   0.8750000   0.8750000   0.6250000), wk =   0.0000000
        k(   31) = (   0.8750000   0.6250000  -0.1250000), wk =   0.1250000
        k(   32) = (   1.8750000   0.6250000  -0.1250000), wk =   0.0000000
        k(   33) = (   0.1250000   0.6250000   0.3750000), wk =   0.1250000
        k(   34) = (   1.1250000   0.6250000   0.3750000), wk =   0.0000000
        k(   35) = (   0.6250000   0.3750000   0.1250000), wk =   0.1250000
        k(   36) = (   1.6250000   0.3750000   0.1250000), wk =   0.0000000
        k(   37) = (  -0.8750000   0.1250000  -0.1250000), wk =   0.0625000
        k(   38) = (   0.1250000   0.1250000  -0.1250000), wk =   0.0000000
        k(   39) = (   1.1250000   0.3750000  -0.3750000), wk =   0.0625000
        k(   40) = (   2.1250000   0.3750000  -0.3750000), wk =   0.0000000

     PseudoPot. # 1 for  C read from file:
     ./C.pz-kjpaw.UPF
     MD5 check sum: 414e6e825ae75add557e798061b49a04
     Pseudo is Projector augmented-wave + core cor, Zval =  4.0
     Generated using "atomic" code by A. Dal Corso (espresso distribution)
     Shape of augmentation charge: BESSEL
     Using radial grid of 1073 points,  4 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
     Q(r) pseudized with 0 coefficients 



     Atomic displacements:
     There are   3 irreducible representations

     Representation     1      2 modes -  To be done

     Representation     2      2 modes -  To be done

     Representation     3      2 modes -  To be done



     Alpha used in Ewald sum =   2.8000
     PHONON       :     1.12s CPU         1.43s WALL



     Representation #  1 modes #   1  2

     Self-consistent Calculation

      iter #   1 total cpu time :     1.8 secs   av.it.:   6.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.331E-08

      iter #   2 total cpu time :     2.1 secs   av.it.:  11.2
      thresh= 1.825E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.995E-09

      iter #   3 total cpu time :     2.4 secs   av.it.:  11.0
      thresh= 5.472E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.689E-10

      iter #   4 total cpu time :     2.8 secs   av.it.:  10.2
      thresh= 1.640E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.900E-13

      iter #   5 total cpu time :     3.1 secs   av.it.:  11.5
      thresh= 6.245E-08 alpha_mix =  0.700 |ddv_scf|^2 =  7.037E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   3  4

     Self-consistent Calculation

      iter #   1 total cpu time :     3.4 secs   av.it.:   7.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.887E-06

      iter #   2 total cpu time :     3.7 secs   av.it.:  11.3
      thresh= 1.699E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.048E-08

      iter #   3 total cpu time :     4.0 secs   av.it.:  11.0
      thresh= 2.459E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.374E-10

      iter #   4 total cpu time :     4.4 secs   av.it.:  10.9
      thresh= 2.894E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.209E-11

      iter #   5 total cpu time :     4.6 secs   av.it.:   9.1
      thresh= 7.217E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.822E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 modes #   5  6

     Self-consistent Calculation

      iter #   1 total cpu time :     4.9 secs   av.it.:   7.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.552E-05

      iter #   2 total cpu time :     5.3 secs   av.it.:  11.1
      thresh= 3.940E-04 alpha_mix =  0.700 |ddv_scf|^2 =  9.905E-06

      iter #   3 total cpu time :     5.6 secs   av.it.:  10.5
      thresh= 3.147E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.810E-08

      iter #   4 total cpu time :     6.1 secs   av.it.:  10.5
      thresh= 1.952E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.355E-10

      iter #   5 total cpu time :     6.4 secs   av.it.:  10.2
      thresh= 2.087E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.463E-11

      iter #   6 total cpu time :     6.7 secs   av.it.:   9.8
      thresh= 4.963E-07 alpha_mix =  0.700 |ddv_scf|^2 =  9.768E-14

      iter #   7 total cpu time :     7.1 secs   av.it.:  10.4
      thresh= 3.125E-08 alpha_mix =  0.700 |ddv_scf|^2 =  4.784E-15

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    3
     List of q in the star:
          1   1.000000000   0.000000000   0.000000000
          2   0.000000000   0.000000000   1.000000000
          3   0.000000000   1.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    1.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =      23.718499 [THz] =     791.163955 [cm-1]
     freq (    2) =      23.718499 [THz] =     791.163955 [cm-1]
     freq (    3) =      31.976134 [THz] =    1066.609020 [cm-1]
     freq (    4) =      31.976134 [THz] =    1066.609020 [cm-1]
     freq (    5) =      36.055300 [THz] =    1202.675339 [cm-1]
     freq (    6) =      36.055300 [THz] =    1202.675339 [cm-1]
 **************************************************************************

     init_run     :      0.11s CPU      0.16s WALL (       1 calls)
     electrons    :      0.19s CPU      0.27s WALL (       1 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       1 calls)
     potinit      :      0.02s CPU      0.04s WALL (       1 calls)

     Called by electrons:
     c_bands      :      0.19s CPU      0.27s WALL (       1 calls)
     v_of_rho     :      0.01s CPU      0.00s WALL (       2 calls)
     newd         :      0.00s CPU      0.01s WALL (       2 calls)
     PAW_pot      :      0.02s CPU      0.04s WALL (       2 calls)

     Called by c_bands:
     init_us_2    :      0.02s CPU      0.02s WALL (     520 calls)
     cegterg      :      0.17s CPU      0.23s WALL (      40 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :      1.34s CPU      1.95s WALL (    8641 calls)
     s_psi        :      0.18s CPU      0.24s WALL (   17414 calls)
     g_psi        :      0.00s CPU      0.00s WALL (     468 calls)
     cdiaghg      :      0.04s CPU      0.06s WALL (     508 calls)

     Called by h_psi:
     h_psi:pot    :      1.33s CPU      1.94s WALL (    8641 calls)
     h_psi:calbec :      0.13s CPU      0.23s WALL (    8641 calls)
     vloc_psi     :      1.06s CPU      1.55s WALL (    8641 calls)
     add_vuspsi   :      0.12s CPU      0.13s WALL (    8641 calls)

     General routines
     calbec       :      0.31s CPU      0.48s WALL (   18974 calls)
     fft          :      0.06s CPU      0.15s WALL (     341 calls)
     ffts         :      0.01s CPU      0.01s WALL (     220 calls)
     fftw         :      0.94s CPU      1.42s WALL (   72498 calls)
     interpolate  :      0.02s CPU      0.04s WALL (      76 calls)
     davcio       :      0.03s CPU      0.04s WALL (    3040 calls)

     Parallel routines
     fft_scatter  :      0.50s CPU      0.87s WALL (   73059 calls)

     PHONON       :     5.17s CPU         7.09s WALL

     INITIALIZATION: 
     phq_setup    :      0.01s CPU      0.01s WALL (       1 calls)
     phq_init     :      0.31s CPU      0.41s WALL (       1 calls)

     phq_init     :      0.31s CPU      0.41s WALL (       1 calls)
     set_drhoc    :      0.14s CPU      0.14s WALL (       3 calls)
     init_vloc    :      0.01s CPU      0.01s WALL (       2 calls)
     init_us_1    :      0.14s CPU      0.17s WALL (       2 calls)
     newd         :      0.00s CPU      0.01s WALL (       2 calls)
     dvanqq       :      0.04s CPU      0.05s WALL (       1 calls)
     drho         :      0.04s CPU      0.07s WALL (       1 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.11s CPU      0.15s WALL (       1 calls)
     phqscf       :      4.05s CPU      5.64s WALL (       1 calls)
     dynmatrix    :      0.00s CPU      0.02s WALL (       1 calls)

     phqscf       :      4.05s CPU      5.64s WALL (       1 calls)
     solve_linter :      4.04s CPU      5.61s WALL (       3 calls)
     drhodv       :      0.01s CPU      0.02s WALL (       3 calls)

     dynmat0      :      0.11s CPU      0.15s WALL (       1 calls)
     dynmat_us    :      0.01s CPU      0.02s WALL (       1 calls)
     d2ionq       :      0.00s CPU      0.00s WALL (       1 calls)
     dynmatcc     :      0.09s CPU      0.10s WALL (       1 calls)

     dynmat_us    :      0.01s CPU      0.02s WALL (       1 calls)
     addusdynmat  :      0.00s CPU      0.00s WALL (       1 calls)

     phqscf       :      4.05s CPU      5.64s WALL (       1 calls)
     solve_linter :      4.04s CPU      5.61s WALL (       3 calls)

     solve_linter :      4.04s CPU      5.61s WALL (       3 calls)
     dvqpsi_us    :      0.04s CPU      0.05s WALL (     120 calls)
     ortho        :      0.05s CPU      0.04s WALL (     680 calls)
     cgsolve      :      1.84s CPU      2.65s WALL (     680 calls)
     incdrhoscf   :      0.07s CPU      0.14s WALL (     680 calls)
     addusddens   :      0.09s CPU      0.16s WALL (      20 calls)
     vpsifft      :      0.08s CPU      0.09s WALL (     560 calls)
     dv_of_drho   :      0.02s CPU      0.05s WALL (      34 calls)
     mix_pot      :      0.03s CPU      0.04s WALL (      17 calls)
     psymdvscf    :      0.97s CPU      1.16s WALL (      17 calls)
     newdq        :      0.12s CPU      0.15s WALL (      17 calls)
     adddvscf     :      0.01s CPU      0.01s WALL (     560 calls)
     drhodvus     :      0.00s CPU      0.00s WALL (       3 calls)

     dvqpsi_us    :      0.04s CPU      0.05s WALL (     120 calls)
     dvqpsi_us_on :      0.02s CPU      0.02s WALL (     120 calls)

     cgsolve      :      1.84s CPU      2.65s WALL (     680 calls)
     ch_psi       :      1.70s CPU      2.47s WALL (    8093 calls)

     ch_psi       :      1.70s CPU      2.47s WALL (    8093 calls)
     h_psi        :      1.34s CPU      1.95s WALL (    8641 calls)
     last         :      0.33s CPU      0.49s WALL (    8093 calls)

     h_psi        :      1.34s CPU      1.95s WALL (    8641 calls)
     add_vuspsi   :      0.12s CPU      0.13s WALL (    8641 calls)

     incdrhoscf   :      0.07s CPU      0.14s WALL (     680 calls)
     addusdbec    :      0.03s CPU      0.03s WALL (     800 calls)

     drhodvus     :      0.00s CPU      0.00s WALL (       3 calls)

      General routines
     calbec       :      0.31s CPU      0.48s WALL (   18974 calls)
     fft          :      0.06s CPU      0.15s WALL (     341 calls)
     ffts         :      0.01s CPU      0.01s WALL (     220 calls)
     fftw         :      0.94s CPU      1.42s WALL (   72498 calls)
     davcio       :      0.03s CPU      0.04s WALL (    3040 calls)
     write_rec    :      0.02s CPU      0.04s WALL (      20 calls)


     PHONON       :     5.17s CPU         7.09s WALL


   This run was terminated on:  13:43:44   7Dec2016            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_51 = """
     Program PHONON v.6.0 (svn rev. 13188M) starts on  7Dec2016 at 13:12: 3 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     R & G space division:  proc/nbgrp/npool/nimage =       4

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/platinum.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         118      55     18                 1712      556     102
     Max         119      56     19                 1715      558     104
     Sum         475     223     73                 6855     2229     411


     Check: negative/imaginary core charge=   -0.000004    0.000000

     Calculation of q =    1.0000000   0.0000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         118      55     22                 1712      556     152
     Max         119      56     23                 1715      558     153
     Sum         475     223     91                 6855     2229     609


     Title: 
     phonons of Pt at X                                                         


     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.4200  a.u.
     unit-cell volume          =     102.1296 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =        10.00
     number of Kohn-Sham states=           18
     kinetic-energy cutoff     =      30.0000  Ry
     charge density cutoff     =     250.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Non magnetic calculation with spin-orbit


     celldm(1)=   7.420000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Pt read from file:
     ./Pt.rel-pz-n-rrkjus.UPF
     MD5 check sum: 4baafe8ec1942611396c7a5466f52249
     Pseudo is Ultrasoft + core correction, Zval = 10.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of 1277 points,  6 beta functions with: 
                l(1) =   2
                l(2) =   2
                l(3) =   2
                l(4) =   2
                l(5) =   1
                l(6) =   1
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Pt            10.00   195.07800     Pt( 1.00)

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Pt  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=     6  Marzari-Vanderbilt smearing, width (Ry)=  0.0200
                       cart. coord. in units 2pi/alat
        k(    1) = (  -0.2500000   0.2500000   0.2500000), wk =   0.2500000
        k(    2) = (   0.7500000   0.2500000   0.2500000), wk =   0.0000000
        k(    3) = (   0.2500000  -0.2500000   0.7500000), wk =   0.5000000
        k(    4) = (   1.2500000  -0.2500000   0.7500000), wk =   0.0000000
        k(    5) = (   0.7500000   0.2500000  -0.2500000), wk =   0.2500000
        k(    6) = (   1.7500000   0.2500000  -0.2500000), wk =   0.0000000

     Dense  grid:     6855 G-vectors     FFT dimensions: (  27,  27,  27)

     Smooth grid:     2229 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       1.79Mb

     Estimated total allocated dynamical RAM >       7.17Mb

     Check: negative/imaginary core charge=   -0.000004    0.000000

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/platinum.save/charge-density.dat

     Starting wfc are   12 atomic +    6 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.00E-10,  avg # of iterations = 15.2

     total cpu time spent up to now is        0.9 secs

     End of band structure calculation

          k =-0.2500 0.2500 0.2500 (   289 PWs)   bands (ev):

     9.3170   9.3170  13.3107  13.3107  13.5800  13.5800  14.7744  14.7744
    16.0692  16.0692  16.6624  16.6624  31.1506  31.1506  35.9701  35.9701
    39.8081  39.8081

          k = 0.7500 0.2500 0.2500 (   283 PWs)   bands (ev):

    11.2910  11.2910  12.4161  12.4161  13.9359  13.9359  15.5889  15.5889
    17.8747  17.8747  20.6641  20.6641  25.0087  25.0087  31.6343  31.6343
    33.8373  33.8373

          k = 0.2500-0.2500 0.7500 (   283 PWs)   bands (ev):

    11.2910  11.2910  12.4161  12.4161  13.9359  13.9359  15.5889  15.5889
    17.8747  17.8747  20.6641  20.6641  25.0087  25.0087  31.6343  31.6343
    33.8373  33.8373

          k = 1.2500-0.2500 0.7500 (   283 PWs)   bands (ev):

    11.2910  11.2910  12.4161  12.4161  13.9359  13.9359  15.5889  15.5889
    17.8747  17.8747  20.6641  20.6641  25.0087  25.0087  31.6343  31.6343
    33.8373  33.8373

          k = 0.7500 0.2500-0.2500 (   283 PWs)   bands (ev):

    11.2910  11.2910  12.4161  12.4161  13.9359  13.9359  15.5889  15.5889
    17.8747  17.8747  20.6641  20.6641  25.0087  25.0087  31.6343  31.6343
    33.8373  33.8373

          k = 1.7500 0.2500-0.2500 (   289 PWs)   bands (ev):

     9.3170   9.3170  13.3107  13.3107  13.5800  13.5800  14.7744  14.7744
    16.0692  16.0692  16.6624  16.6624  31.1506  31.1506  35.9701  35.9701
    39.8081  39.8081

     the Fermi energy is    17.9731 ev

     Writing output data file platinum.save

     phonons of Pt at X                                                         

     bravais-lattice index     =            2
     lattice parameter (alat)  =       7.4200  a.u.
     unit-cell volume          =     102.1296 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      30.0000  Ry
     charge density cut-off    =     250.0000  Ry
     convergence threshold     =      1.0E-16
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Non magnetic calculation with spin-orbit

     celldm(1)=    7.42000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Pt 195.0780   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   1.0000000   0.0000000   0.0000000 )

     17 Sym.Ops. (with q -> -q+G )


     G cutoff =  348.6487  (   1715 G-vectors)     FFT grid: ( 27, 27, 27)
     G cutoff =  167.3514  (    557 G-vectors)  smooth grid: ( 20, 20, 20)

     number of k points=     6  Marzari-Vanderbilt smearing, width (Ry)=  0.0200
                       cart. coord. in units 2pi/alat
        k(    1) = (  -0.2500000   0.2500000   0.2500000), wk =   0.2500000
        k(    2) = (   0.7500000   0.2500000   0.2500000), wk =   0.0000000
        k(    3) = (   0.2500000  -0.2500000   0.7500000), wk =   0.5000000
        k(    4) = (   1.2500000  -0.2500000   0.7500000), wk =   0.0000000
        k(    5) = (   0.7500000   0.2500000  -0.2500000), wk =   0.2500000
        k(    6) = (   1.7500000   0.2500000  -0.2500000), wk =   0.0000000

     PseudoPot. # 1 for Pt read from file:
     ./Pt.rel-pz-n-rrkjus.UPF
     MD5 check sum: 4baafe8ec1942611396c7a5466f52249
     Pseudo is Ultrasoft + core correction, Zval = 10.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of 1277 points,  6 beta functions with: 
                l(1) =   2
                l(2) =   2
                l(3) =   2
                l(4) =   2
                l(5) =   1
                l(6) =   1
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, D_4h(4/mmm) point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_2u X_4' M_4'  To be done

     Representation     2      2 modes -E_u  X_5' M_5'  To be done



     Alpha used in Ewald sum =   2.6000
     PHONON       :     1.96s CPU         2.34s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     2.5 secs   av.it.:   8.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.837E-04

      iter #   2 total cpu time :     2.7 secs   av.it.:  10.7
      thresh= 2.199E-03 alpha_mix =  0.700 |ddv_scf|^2 =  5.208E-04

      iter #   3 total cpu time :     2.9 secs   av.it.:   9.3
      thresh= 2.282E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.227E-08

      iter #   4 total cpu time :     3.1 secs   av.it.:  10.7
      thresh= 1.492E-05 alpha_mix =  0.700 |ddv_scf|^2 =  1.950E-10

      iter #   5 total cpu time :     3.4 secs   av.it.:   9.7
      thresh= 1.396E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.734E-12

      iter #   6 total cpu time :     3.6 secs   av.it.:  10.3
      thresh= 2.176E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.745E-15

      iter #   7 total cpu time :     3.8 secs   av.it.:  10.7
      thresh= 8.213E-09 alpha_mix =  0.700 |ddv_scf|^2 =  2.979E-16

      iter #   8 total cpu time :     4.0 secs   av.it.:  10.0
      thresh= 1.726E-09 alpha_mix =  0.700 |ddv_scf|^2 =  1.378E-18

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   2  3

     Self-consistent Calculation

      iter #   1 total cpu time :     4.4 secs   av.it.:   7.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  7.681E-06

      iter #   2 total cpu time :     4.7 secs   av.it.:  11.5
      thresh= 2.771E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.049E-06

      iter #   3 total cpu time :     5.1 secs   av.it.:  11.2
      thresh= 1.024E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.459E-09

      iter #   4 total cpu time :     5.5 secs   av.it.:  11.2
      thresh= 4.959E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.590E-12

      iter #   5 total cpu time :     5.8 secs   av.it.:  11.0
      thresh= 2.142E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.683E-14

      iter #   6 total cpu time :     6.2 secs   av.it.:  11.0
      thresh= 2.164E-08 alpha_mix =  0.700 |ddv_scf|^2 =  7.103E-16

      iter #   7 total cpu time :     6.5 secs   av.it.:  10.7
      thresh= 2.665E-09 alpha_mix =  0.700 |ddv_scf|^2 =  5.031E-19

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    3
     List of q in the star:
          1   1.000000000   0.000000000   0.000000000
          2   0.000000000   0.000000000   1.000000000
          3   0.000000000   1.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    1.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       3.691550 [THz] =     123.136841 [cm-1]
     freq (    2) =       3.691550 [THz] =     123.136841 [cm-1]
     freq (    3) =       5.815216 [THz] =     193.974741 [cm-1]
 **************************************************************************

     Mode symmetry, D_4h(4/mmm) point group:

     freq (  1 -  2) =        123.1  [cm-1]   --> E_u  X_5' M_5'     
     freq (  3 -  3) =        194.0  [cm-1]   --> A_2u X_4' M_4'     

     init_run     :      0.25s CPU      0.29s WALL (       1 calls)
     electrons    :      0.50s CPU      0.66s WALL (       1 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       1 calls)
     potinit      :      0.00s CPU      0.01s WALL (       1 calls)

     Called by electrons:
     c_bands      :      0.50s CPU      0.65s WALL (       1 calls)
     v_of_rho     :      0.00s CPU      0.00s WALL (       2 calls)
     newd         :      0.05s CPU      0.06s WALL (       2 calls)

     Called by c_bands:
     init_us_2    :      0.00s CPU      0.00s WALL (      69 calls)
     cegterg      :      0.45s CPU      0.59s WALL (       6 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :      1.16s CPU      1.47s WALL (    1007 calls)
     s_psi        :      0.27s CPU      0.34s WALL (    1977 calls)
     g_psi        :      0.00s CPU      0.00s WALL (      91 calls)
     cdiaghg      :      0.16s CPU      0.18s WALL (      97 calls)

     Called by h_psi:
     h_psi:pot    :      1.15s CPU      1.46s WALL (    1007 calls)
     h_psi:calbec :      0.15s CPU      0.17s WALL (    1007 calls)
     vloc_psi     :      0.87s CPU      1.13s WALL (    1007 calls)
     add_vuspsi   :      0.13s CPU      0.15s WALL (    1007 calls)

     General routines
     calbec       :      0.32s CPU      0.37s WALL (    2130 calls)
     fft          :      0.10s CPU      0.19s WALL (     221 calls)
     ffts         :      0.00s CPU      0.00s WALL (      73 calls)
     fftw         :      0.83s CPU      1.05s WALL (   41356 calls)
     interpolate  :      0.00s CPU      0.01s WALL (      55 calls)
     davcio       :      0.01s CPU      0.02s WALL (     455 calls)

     Parallel routines
     fft_scatter  :      0.48s CPU      0.67s WALL (   41650 calls)

     PHONON       :     5.56s CPU         6.52s WALL

     INITIALIZATION: 
     phq_setup    :      0.00s CPU      0.01s WALL (       1 calls)
     phq_init     :      0.65s CPU      0.72s WALL (       1 calls)

     phq_init     :      0.65s CPU      0.72s WALL (       1 calls)
     set_drhoc    :      0.16s CPU      0.16s WALL (       3 calls)
     init_vloc    :      0.01s CPU      0.01s WALL (       2 calls)
     init_us_1    :      0.45s CPU      0.58s WALL (       2 calls)
     newd         :      0.05s CPU      0.06s WALL (       2 calls)
     dvanqq       :      0.14s CPU      0.15s WALL (       1 calls)
     drho         :      0.26s CPU      0.31s WALL (       1 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.14s CPU      0.15s WALL (       1 calls)
     phqscf       :      3.59s CPU      4.17s WALL (       1 calls)
     dynmatrix    :      0.00s CPU      0.01s WALL (       1 calls)

     phqscf       :      3.59s CPU      4.17s WALL (       1 calls)
     solve_linter :      3.58s CPU      4.15s WALL (       2 calls)
     drhodv       :      0.02s CPU      0.02s WALL (       2 calls)

     dynmat0      :      0.14s CPU      0.15s WALL (       1 calls)
     dynmat_us    :      0.04s CPU      0.04s WALL (       1 calls)
     d2ionq       :      0.00s CPU      0.00s WALL (       1 calls)
     dynmatcc     :      0.10s CPU      0.10s WALL (       1 calls)

     dynmat_us    :      0.04s CPU      0.04s WALL (       1 calls)
     addusdynmat  :      0.01s CPU      0.01s WALL (       1 calls)

     phqscf       :      3.59s CPU      4.17s WALL (       1 calls)
     solve_linter :      3.58s CPU      4.15s WALL (       2 calls)

     solve_linter :      3.58s CPU      4.15s WALL (       2 calls)
     dvqpsi_us    :      0.07s CPU      0.09s WALL (       9 calls)
     ortho        :      0.04s CPU      0.07s WALL (      66 calls)
     cgsolve      :      1.58s CPU      1.91s WALL (      66 calls)
     incdrhoscf   :      0.11s CPU      0.10s WALL (      66 calls)
     addusddens   :      0.58s CPU      0.63s WALL (      17 calls)
     vpsifft      :      0.07s CPU      0.08s WALL (      57 calls)
     dv_of_drho   :      0.00s CPU      0.03s WALL (      22 calls)
     mix_pot      :      0.02s CPU      0.01s WALL (      15 calls)
     psymdvscf    :      0.39s CPU      0.42s WALL (      15 calls)
     newdq        :      0.68s CPU      0.76s WALL (      15 calls)
     adddvscf     :      0.03s CPU      0.03s WALL (      57 calls)
     drhodvus     :      0.00s CPU      0.00s WALL (       2 calls)

     dvqpsi_us    :      0.07s CPU      0.09s WALL (       9 calls)
     dvqpsi_us_on :      0.06s CPU      0.07s WALL (       9 calls)

     cgsolve      :      1.58s CPU      1.91s WALL (      66 calls)
     ch_psi       :      1.50s CPU      1.82s WALL (     904 calls)

     ch_psi       :      1.50s CPU      1.82s WALL (     904 calls)
     h_psi        :      1.16s CPU      1.47s WALL (    1007 calls)
     last         :      0.36s CPU      0.47s WALL (     904 calls)

     h_psi        :      1.16s CPU      1.47s WALL (    1007 calls)
     add_vuspsi   :      0.13s CPU      0.15s WALL (    1007 calls)

     incdrhoscf   :      0.11s CPU      0.10s WALL (      66 calls)

     drhodvus     :      0.00s CPU      0.00s WALL (       2 calls)

      General routines
     calbec       :      0.32s CPU      0.37s WALL (    2130 calls)
     fft          :      0.10s CPU      0.19s WALL (     221 calls)
     ffts         :      0.00s CPU      0.00s WALL (      73 calls)
     fftw         :      0.83s CPU      1.05s WALL (   41356 calls)
     davcio       :      0.01s CPU      0.02s WALL (     455 calls)
     write_rec    :      0.01s CPU      0.03s WALL (      17 calls)


     PHONON       :     5.56s CPU         6.52s WALL


   This run was terminated on:  13:12:10   7Dec2016            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_56 = """
     Program PHONON v.6.0 (svn rev. 13286) starts on  7Feb2017 at 14:28:57 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     path-images division:  nimage    =       2
     R & G space division:  proc/nbgrp/npool/nimage =       2

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/alas.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     229
     Max         121     121     43                 1224     1224     230
     Sum         241     241     85                 2445     2445     459



     Dynamical matrices for ( 4, 4, 4)  uniform grid of q-points
     (   8q-points):
       N         xq(1)         xq(2)         xq(3) 
       1   0.000000000   0.000000000   0.000000000
       2  -0.250000000   0.250000000  -0.250000000
       3   0.500000000  -0.500000000   0.500000000
       4   0.000000000   0.500000000   0.000000000
       5   0.750000000  -0.250000000   0.750000000
       6   0.500000000   0.000000000   0.500000000
       7   0.000000000  -1.000000000   0.000000000
       8  -0.500000000  -1.000000000   0.000000000

      Image parallelization. There are  2 images and    38 representations
      The estimated total work is   336 self-consistent (scf) runs
      I am image number     0 and my work is about  165 scf runs. I calculate: 
      q point number     1, representations:
       0 1 2
      q point number     2, representations:
       0 1 2 3 4
      q point number     3, representations:
       0 1 2 3 4
      q point number     4, representations:
       0 1 2 3 4 5 6
      q point number     5, representations:
       0 1 2 3 4

     Calculation of q =    0.0000000   0.0000000   0.0000000

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   0.0000000 )

     25 Sym.Ops. (with q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=     2

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, T_d (-43m)  point group:


     Electric field:
     Dielectric constant
     Born effective charges in two ways 


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      3 modes -T_2  G_15 P_4  To be done

     Representation     2      3 modes -T_2  G_15 P_4  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     0.25s CPU         0.26s WALL


     Electric Fields Calculation

      iter #   1 total cpu time :     0.3 secs   av.it.:   6.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.326E-06

      iter #   2 total cpu time :     0.4 secs   av.it.:   9.3
      thresh= 1.151E-04 alpha_mix =  0.700 |ddv_scf|^2 =  6.508E-08

      iter #   3 total cpu time :     0.4 secs   av.it.:   9.5
      thresh= 2.551E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.401E-10

      iter #   4 total cpu time :     0.5 secs   av.it.:   9.8
      thresh= 2.530E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.108E-12

      iter #   5 total cpu time :     0.5 secs   av.it.:   9.0
      thresh= 1.763E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.543E-14

     End of electric fields calculation

          Dielectric constant in cartesian axis 

          (      13.744216098      -0.000000000      -0.000000000 )
          (       0.000000000      13.744216098       0.000000000 )
          (      -0.000000000       0.000000000      13.744216098 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   Al 
      Ex  (        1.88265       -0.00000       -0.00000 )
      Ey  (       -0.00000        1.88265       -0.00000 )
      Ez  (        0.00000        0.00000        1.88265 )
           atom      2   As 
      Ex  (       -3.23374       -0.00000       -0.00000 )
      Ey  (        0.00000       -3.23374       -0.00000 )
      Ez  (       -0.00000       -0.00000       -3.23374 )


     Representation #  1 modes #   1  2  3

     Self-consistent Calculation

      iter #   1 total cpu time :     0.6 secs   av.it.:   5.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.662E-07

      iter #   2 total cpu time :     0.7 secs   av.it.:   9.7
      thresh= 6.828E-05 alpha_mix =  0.700 |ddv_scf|^2 =  2.273E-08

      iter #   3 total cpu time :     0.7 secs   av.it.:   9.7
      thresh= 1.508E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.797E-11

      iter #   4 total cpu time :     0.8 secs   av.it.:   9.5
      thresh= 6.162E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.182E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   4  5  6

     Self-consistent Calculation

      iter #   1 total cpu time :     0.8 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.910E-08

      iter #   2 total cpu time :     0.9 secs   av.it.:   9.8
      thresh= 1.706E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.259E-10

      iter #   3 total cpu time :     0.9 secs   av.it.:   9.5
      thresh= 1.805E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.012E-11

      iter #   4 total cpu time :     1.0 secs   av.it.:   9.5
      thresh= 5.488E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.306E-12

      iter #   5 total cpu time :     1.1 secs   av.it.:   9.5
      thresh= 1.143E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.628E-16

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

          Dielectric constant in cartesian axis 

          (      13.744216098      -0.000000000      -0.000000000 )
          (       0.000000000      13.744216098       0.000000000 )
          (      -0.000000000       0.000000000      13.744216098 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   Al 
      Ex  (        1.88265       -0.00000       -0.00000 )
      Ey  (       -0.00000        1.88265       -0.00000 )
      Ez  (        0.00000        0.00000        1.88265 )
           atom      2   As 
      Ex  (       -3.23374       -0.00000       -0.00000 )
      Ey  (        0.00000       -3.23374       -0.00000 )
      Ez  (       -0.00000       -0.00000       -3.23374 )

          Effective charges (d P / du) in cartesian axis 

           atom      1   Al 
      Px  (        1.88284       -0.00000       -0.00000 )
      Py  (       -0.00000        1.88284        0.00000 )
      Pz  (       -0.00000        0.00000        1.88284 )
           atom      2   As 
      Px  (       -3.23837       -0.00000        0.00000 )
      Py  (       -0.00000       -3.23837       -0.00000 )
      Pz  (        0.00000       -0.00000       -3.23837 )

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    2) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    3) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    4) =      11.258806 [THz] =     375.553348 [cm-1]
     freq (    5) =      11.258806 [THz] =     375.553348 [cm-1]
     freq (    6) =      11.258806 [THz] =     375.553348 [cm-1]
 **************************************************************************

     Mode symmetry, T_d (-43m)  point group:

     freq (  1 -  3) =          5.5  [cm-1]   --> T_2  G_15 P_4   I+R
     freq (  4 -  6) =        375.6  [cm-1]   --> T_2  G_15 P_4   I+R

     Calculation of q =   -0.2500000   0.2500000  -0.2500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     264
     Max         121     121     43                 1224     1224     267
     Sum         241     241     85                 2445     2445     531



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    20
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.1875000
        k(    2) = (   0.0000000   0.5000000   0.0000000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.3750000
        k(    4) = (   0.0000000   0.5000000   0.5000000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0625000
        k(    6) = (  -0.5000000   0.5000000  -0.5000000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.1875000
        k(    8) = (  -0.5000000   0.0000000  -0.5000000), wk =   0.0000000
        k(    9) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0625000
        k(   10) = (   0.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   11) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.1875000
        k(   12) = (  -1.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   13) = (  -0.7500000   0.2500000  -0.2500000), wk =   0.1875000
        k(   14) = (  -1.0000000   0.5000000  -0.5000000), wk =   0.0000000
        k(   15) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.3750000
        k(   16) = (  -0.5000000   0.0000000  -1.0000000), wk =   0.0000000
        k(   17) = (   0.2500000  -0.2500000   0.7500000), wk =   0.1875000
        k(   18) = (   0.0000000   0.0000000   0.5000000), wk =   0.0000000
        k(   19) = (  -0.2500000   0.2500000   0.7500000), wk =   0.1875000
        k(   20) = (  -0.5000000   0.5000000   0.5000000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.0

     total cpu time spent up to now is        0.2 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.0000 0.5000 0.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.0000 0.5000 0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.5000 0.5000-0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.5000 0.0000-0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k = 0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.0000 0.0000 0.0000 (   331 PWs)   bands (ev):

    -6.9797   5.1761   5.1761   5.1761

          k =-0.7500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-1.0000 0.0000 0.0000 (   302 PWs)   bands (ev):

    -4.8217  -0.4470   2.9274   2.9274

          k =-0.7500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-1.0000 0.5000-0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.5000 0.0000-1.0000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k = 0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.0000 0.0000 0.5000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.5000 0.5000 0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (  -0.2500000   0.2500000  -0.2500000 )

      6 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    20

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      1 modes -A_1  L_1  To be done

     Representation     3      2 modes -E    L_3  To be done

     Representation     4      2 modes -E    L_3  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     1.28s CPU         1.36s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     1.4 secs   av.it.:   6.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.084E-03

      iter #   2 total cpu time :     1.4 secs   av.it.:   7.6
      thresh= 5.554E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.965E-02

      iter #   3 total cpu time :     1.5 secs   av.it.:   6.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.334E-06

      iter #   4 total cpu time :     1.5 secs   av.it.:   7.2
      thresh= 2.517E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.380E-07

      iter #   5 total cpu time :     1.6 secs   av.it.:   7.6
      thresh= 3.714E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.454E-09

      iter #   6 total cpu time :     1.6 secs   av.it.:   7.0
      thresh= 6.674E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.717E-10

      iter #   7 total cpu time :     1.7 secs   av.it.:   7.2
      thresh= 2.172E-06 alpha_mix =  0.700 |ddv_scf|^2 =  3.731E-11

      iter #   8 total cpu time :     1.7 secs   av.it.:   7.2
      thresh= 6.108E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.132E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     1.8 secs   av.it.:   5.6
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  6.514E-04

      iter #   2 total cpu time :     1.8 secs   av.it.:   7.6
      thresh= 2.552E-03 alpha_mix =  0.700 |ddv_scf|^2 =  5.928E-03

      iter #   3 total cpu time :     1.9 secs   av.it.:   6.2
      thresh= 7.699E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.208E-07

      iter #   4 total cpu time :     1.9 secs   av.it.:   8.2
      thresh= 4.699E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.970E-09

      iter #   5 total cpu time :     2.0 secs   av.it.:   8.1
      thresh= 8.348E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.212E-10

      iter #   6 total cpu time :     2.0 secs   av.it.:   7.4
      thresh= 2.283E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.038E-09

      iter #   7 total cpu time :     2.1 secs   av.it.:   6.9
      thresh= 3.223E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.551E-11

      iter #   8 total cpu time :     2.1 secs   av.it.:   7.6
      thresh= 3.938E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.433E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 modes #   3  4

     Self-consistent Calculation

      iter #   1 total cpu time :     2.2 secs   av.it.:   5.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.311E-06

      iter #   2 total cpu time :     2.3 secs   av.it.:   9.2
      thresh= 1.145E-04 alpha_mix =  0.700 |ddv_scf|^2 =  9.088E-08

      iter #   3 total cpu time :     2.4 secs   av.it.:   9.2
      thresh= 3.015E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.906E-11

      iter #   4 total cpu time :     2.5 secs   av.it.:   9.2
      thresh= 9.437E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.651E-12

      iter #   5 total cpu time :     2.6 secs   av.it.:   9.0
      thresh= 1.285E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.874E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 modes #   5  6

     Self-consistent Calculation

      iter #   1 total cpu time :     2.7 secs   av.it.:   5.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.108E-07

      iter #   2 total cpu time :     2.8 secs   av.it.:   9.4
      thresh= 3.328E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.511E-09

      iter #   3 total cpu time :     2.9 secs   av.it.:   9.2
      thresh= 6.717E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.323E-10

      iter #   4 total cpu time :     3.0 secs   av.it.:   9.1
      thresh= 1.150E-06 alpha_mix =  0.700 |ddv_scf|^2 =  6.943E-12

      iter #   5 total cpu time :     3.1 secs   av.it.:   8.8
      thresh= 2.635E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.124E-15

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    4
     List of q in the star:
          1  -0.250000000   0.250000000  -0.250000000
          2   0.250000000  -0.250000000  -0.250000000
          3  -0.250000000  -0.250000000   0.250000000
          4   0.250000000   0.250000000   0.250000000
     In addition there is the -q list: 
          1   0.250000000  -0.250000000   0.250000000
          2  -0.250000000   0.250000000   0.250000000
          3   0.250000000   0.250000000  -0.250000000
          4  -0.250000000  -0.250000000  -0.250000000

     Diagonalizing the dynamical matrix

     q = (   -0.250000000   0.250000000  -0.250000000 ) 

 **************************************************************************
     freq (    1) =       1.761214 [THz] =      58.747782 [cm-1]
     freq (    2) =       1.761214 [THz] =      58.747782 [cm-1]
     freq (    3) =       4.534095 [THz] =     151.241127 [cm-1]
     freq (    4) =      11.004844 [THz] =     367.082097 [cm-1]
     freq (    5) =      11.004844 [THz] =     367.082097 [cm-1]
     freq (    6) =      12.136604 [THz] =     404.833529 [cm-1]
 **************************************************************************

     Mode symmetry, C_3v (3m)   point group:

     freq (  1 -  2) =         58.7  [cm-1]   --> E    L_3           
     freq (  3 -  3) =        151.2  [cm-1]   --> A_1  L_1           
     freq (  4 -  5) =        367.1  [cm-1]   --> E    L_3           
     freq (  6 -  6) =        404.8  [cm-1]   --> A_1  L_1           

     Calculation of q =    0.5000000  -0.5000000   0.5000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     267
     Max         121     121     43                 1224     1224     270
     Sum         241     241     85                 2445     2445     537



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    10
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.3750000
        k(    2) = (   0.7500000  -0.2500000   0.7500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.7500000
        k(    4) = (   0.7500000  -0.2500000   1.2500000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(    6) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0000000
        k(    7) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.3750000
        k(    8) = (  -0.2500000  -0.7500000   0.7500000), wk =   0.0000000
        k(    9) = (  -0.7500000   0.2500000  -0.2500000), wk =   0.3750000
        k(   10) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.6

     total cpu time spent up to now is        0.3 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.7500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.7500-0.2500 1.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.7500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.7500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.7500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.5000000  -0.5000000   0.5000000 )

      7 Sym.Ops. (with q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    10

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      1 modes -A_1  L_1  To be done

     Representation     2      1 modes -A_1  L_1  To be done

     Representation     3      2 modes -E    L_3  To be done

     Representation     4      2 modes -E    L_3  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     3.13s CPU         3.33s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     3.3 secs   av.it.:   6.2
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.569E-04

      iter #   2 total cpu time :     3.4 secs   av.it.:   8.2
      thresh= 1.889E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.022E-03

      iter #   3 total cpu time :     3.4 secs   av.it.:   7.4
      thresh= 3.197E-03 alpha_mix =  0.700 |ddv_scf|^2 =  5.259E-08

      iter #   4 total cpu time :     3.4 secs   av.it.:   8.0
      thresh= 2.293E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.120E-09

      iter #   5 total cpu time :     3.5 secs   av.it.:   7.4
      thresh= 9.011E-06 alpha_mix =  0.700 |ddv_scf|^2 =  4.293E-11

      iter #   6 total cpu time :     3.5 secs   av.it.:   8.4
      thresh= 6.552E-07 alpha_mix =  0.700 |ddv_scf|^2 =  4.553E-12

      iter #   7 total cpu time :     3.5 secs   av.it.:   8.0
      thresh= 2.134E-07 alpha_mix =  0.700 |ddv_scf|^2 =  8.070E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     3.5 secs   av.it.:   5.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  5.799E-05

      iter #   2 total cpu time :     3.6 secs   av.it.:   8.2
      thresh= 7.615E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.451E-04

      iter #   3 total cpu time :     3.6 secs   av.it.:   7.4
      thresh= 1.205E-03 alpha_mix =  0.700 |ddv_scf|^2 =  6.734E-07

      iter #   4 total cpu time :     3.6 secs   av.it.:   7.6
      thresh= 8.206E-05 alpha_mix =  0.700 |ddv_scf|^2 =  4.035E-09

      iter #   5 total cpu time :     3.7 secs   av.it.:   8.0
      thresh= 6.352E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.159E-11

      iter #   6 total cpu time :     3.7 secs   av.it.:   8.4
      thresh= 8.461E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.226E-12

      iter #   7 total cpu time :     3.7 secs   av.it.:   8.2
      thresh= 1.107E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.358E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 modes #   3  4

     Self-consistent Calculation

      iter #   1 total cpu time :     3.8 secs   av.it.:   6.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.539E-06

      iter #   2 total cpu time :     3.8 secs   av.it.:   9.2
      thresh= 1.241E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.304E-07

      iter #   3 total cpu time :     3.9 secs   av.it.:   9.0
      thresh= 3.611E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.951E-11

      iter #   4 total cpu time :     3.9 secs   av.it.:   9.2
      thresh= 9.461E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.026E-13

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 modes #   5  6

     Self-consistent Calculation

      iter #   1 total cpu time :     4.0 secs   av.it.:   4.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.480E-07

      iter #   2 total cpu time :     4.0 secs   av.it.:   9.0
      thresh= 3.847E-05 alpha_mix =  0.700 |ddv_scf|^2 =  8.827E-09

      iter #   3 total cpu time :     4.1 secs   av.it.:   9.0
      thresh= 9.395E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.470E-10

      iter #   4 total cpu time :     4.2 secs   av.it.:   9.1
      thresh= 1.212E-06 alpha_mix =  0.700 |ddv_scf|^2 =  7.522E-12

      iter #   5 total cpu time :     4.2 secs   av.it.:   8.3
      thresh= 2.743E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.741E-15

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    4
     List of q in the star:
          1   0.500000000  -0.500000000   0.500000000
          2  -0.500000000   0.500000000   0.500000000
          3   0.500000000   0.500000000  -0.500000000
          4  -0.500000000  -0.500000000  -0.500000000

     Diagonalizing the dynamical matrix

     q = (    0.500000000  -0.500000000   0.500000000 ) 

 **************************************************************************
     freq (    1) =       2.016390 [THz] =      67.259545 [cm-1]
     freq (    2) =       2.016390 [THz] =      67.259545 [cm-1]
     freq (    3) =       6.494357 [THz] =     216.628437 [cm-1]
     freq (    4) =      10.940872 [THz] =     364.948217 [cm-1]
     freq (    5) =      10.940872 [THz] =     364.948217 [cm-1]
     freq (    6) =      11.551694 [THz] =     385.323024 [cm-1]
 **************************************************************************

     Mode symmetry, C_3v (3m)   point group:

     freq (  1 -  2) =         67.3  [cm-1]   --> E    L_3           
     freq (  3 -  3) =        216.6  [cm-1]   --> A_1  L_1           
     freq (  4 -  5) =        364.9  [cm-1]   --> E    L_3           
     freq (  6 -  6) =        385.3  [cm-1]   --> A_1  L_1           

     Calculation of q =    0.0000000   0.5000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     264
     Max         121     121     43                 1224     1224     267
     Sum         241     241     85                 2445     2445     531



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    24
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.1250000
        k(    2) = (   0.2500000   0.7500000   0.2500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.2500000
        k(    4) = (   0.2500000   0.7500000   0.7500000), wk =   0.0000000
        k(    5) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.1250000
        k(    6) = (  -0.2500000   0.2500000   0.2500000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.1250000
        k(    8) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0000000
        k(    9) = (   0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(   10) = (   0.2500000   0.7500000  -0.2500000), wk =   0.0000000
        k(   11) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.2500000
        k(   12) = (  -0.2500000   0.2500000   0.7500000), wk =   0.0000000
        k(   13) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.1250000
        k(   14) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.0000000
        k(   15) = (  -0.2500000   0.7500000  -0.2500000), wk =   0.1250000
        k(   16) = (  -0.2500000   1.2500000  -0.2500000), wk =   0.0000000
        k(   17) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.2500000
        k(   18) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.0000000
        k(   19) = (   0.2500000   0.2500000  -0.7500000), wk =   0.2500000
        k(   20) = (   0.2500000   0.7500000  -0.7500000), wk =   0.0000000
        k(   21) = (  -0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   22) = (  -0.2500000   1.2500000   0.2500000), wk =   0.0000000
        k(   23) = (  -0.2500000  -0.7500000  -0.2500000), wk =   0.1250000
        k(   24) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.7

     total cpu time spent up to now is        0.6 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.7500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.2500 0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k =-0.2500 0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 1.2500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.2500 0.7500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500 1.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.5000000   0.0000000 )

      4 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    24

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  To be done

     Representation     2      1 modes -A_1  D_1  S_1  To be done

     Representation     3      1 modes -B_1  D_3  S_3  To be done

     Representation     4      1 modes -B_1  D_3  S_3  To be done

     Representation     5      1 modes -B_2  D_4  S_4  To be done

     Representation     6      1 modes -B_2  D_4  S_4  To be done



     Alpha used in Ewald sum =   0.7000
     PHONON       :     4.23s CPU         4.55s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     4.6 secs   av.it.:   6.5
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.919E-03

      iter #   2 total cpu time :     4.7 secs   av.it.:   8.0
      thresh= 4.381E-03 alpha_mix =  0.700 |ddv_scf|^2 =  1.597E-02

      iter #   3 total cpu time :     4.7 secs   av.it.:   7.1
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  2.108E-06

      iter #   4 total cpu time :     4.8 secs   av.it.:   8.3
      thresh= 1.452E-04 alpha_mix =  0.700 |ddv_scf|^2 =  2.578E-08

      iter #   5 total cpu time :     4.8 secs   av.it.:   8.7
      thresh= 1.606E-05 alpha_mix =  0.700 |ddv_scf|^2 =  7.158E-11

      iter #   6 total cpu time :     4.9 secs   av.it.:   8.2
      thresh= 8.460E-07 alpha_mix =  0.700 |ddv_scf|^2 =  9.928E-11

      iter #   7 total cpu time :     4.9 secs   av.it.:   7.1
      thresh= 9.964E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.226E-11

      iter #   8 total cpu time :     5.0 secs   av.it.:   7.2
      thresh= 5.679E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.114E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 mode #   2

     Self-consistent Calculation

      iter #   1 total cpu time :     5.0 secs   av.it.:   5.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  3.751E-04

      iter #   2 total cpu time :     5.1 secs   av.it.:   8.0
      thresh= 1.937E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.908E-03

      iter #   3 total cpu time :     5.2 secs   av.it.:   6.7
      thresh= 5.393E-03 alpha_mix =  0.700 |ddv_scf|^2 =  5.668E-07

      iter #   4 total cpu time :     5.2 secs   av.it.:   7.8
      thresh= 7.528E-05 alpha_mix =  0.700 |ddv_scf|^2 =  5.601E-09

      iter #   5 total cpu time :     5.3 secs   av.it.:   8.7
      thresh= 7.484E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.285E-11

      iter #   6 total cpu time :     5.3 secs   av.it.:   8.3
      thresh= 7.270E-07 alpha_mix =  0.700 |ddv_scf|^2 =  7.730E-12

      iter #   7 total cpu time :     5.4 secs   av.it.:   7.9
      thresh= 2.780E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.587E-11

      iter #   8 total cpu time :     5.4 secs   av.it.:   6.9
      thresh= 3.983E-07 alpha_mix =  0.700 |ddv_scf|^2 =  2.867E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  3 mode #   3

     Self-consistent Calculation

      iter #   1 total cpu time :     5.5 secs   av.it.:   5.7
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.766E-06

      iter #   2 total cpu time :     5.6 secs   av.it.:   8.4
      thresh= 2.961E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.179E-06

      iter #   3 total cpu time :     5.6 secs   av.it.:   8.2
      thresh= 1.086E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.850E-10

      iter #   4 total cpu time :     5.7 secs   av.it.:   8.0
      thresh= 1.962E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.528E-11

      iter #   5 total cpu time :     5.7 secs   av.it.:   8.2
      thresh= 3.908E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.631E-14

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  4 mode #   4

     Self-consistent Calculation

      iter #   1 total cpu time :     5.8 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.095E-06

      iter #   2 total cpu time :     5.8 secs   av.it.:   8.4
      thresh= 1.046E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.263E-07

      iter #   3 total cpu time :     5.9 secs   av.it.:   8.3
      thresh= 3.554E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.599E-10

      iter #   4 total cpu time :     6.0 secs   av.it.:   7.9
      thresh= 2.569E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.110E-11

      iter #   5 total cpu time :     6.0 secs   av.it.:   7.9
      thresh= 4.594E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.831E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  5 mode #   5

     Self-consistent Calculation

      iter #   1 total cpu time :     6.1 secs   av.it.:   4.9
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.095E-06

      iter #   2 total cpu time :     6.1 secs   av.it.:   8.4
      thresh= 1.046E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.261E-07

      iter #   3 total cpu time :     6.2 secs   av.it.:   8.2
      thresh= 3.552E-05 alpha_mix =  0.700 |ddv_scf|^2 =  6.596E-10

      iter #   4 total cpu time :     6.2 secs   av.it.:   7.9
      thresh= 2.568E-06 alpha_mix =  0.700 |ddv_scf|^2 =  2.115E-11

      iter #   5 total cpu time :     6.3 secs   av.it.:   7.8
      thresh= 4.599E-07 alpha_mix =  0.700 |ddv_scf|^2 =  3.785E-15

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  6 mode #   6

     Self-consistent Calculation

      iter #   1 total cpu time :     6.3 secs   av.it.:   5.8
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  8.766E-06

      iter #   2 total cpu time :     6.4 secs   av.it.:   8.4
      thresh= 2.961E-04 alpha_mix =  0.700 |ddv_scf|^2 =  1.179E-06

      iter #   3 total cpu time :     6.5 secs   av.it.:   8.1
      thresh= 1.086E-04 alpha_mix =  0.700 |ddv_scf|^2 =  3.848E-10

      iter #   4 total cpu time :     6.5 secs   av.it.:   8.0
      thresh= 1.962E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.525E-11

      iter #   5 total cpu time :     6.6 secs   av.it.:   8.2
      thresh= 3.905E-07 alpha_mix =  0.700 |ddv_scf|^2 =  6.624E-14

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    6
     List of q in the star:
          1   0.000000000   0.500000000   0.000000000
          2  -0.500000000   0.000000000   0.000000000
          3   0.000000000  -0.500000000   0.000000000
          4   0.000000000   0.000000000   0.500000000
          5   0.000000000   0.000000000  -0.500000000
          6   0.500000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.500000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       2.421101 [THz] =      80.759230 [cm-1]
     freq (    2) =       2.421101 [THz] =      80.759230 [cm-1]
     freq (    3) =       4.606324 [THz] =     153.650423 [cm-1]
     freq (    4) =      10.666710 [THz] =     355.803149 [cm-1]
     freq (    5) =      10.666710 [THz] =     355.803149 [cm-1]
     freq (    6) =      12.371391 [THz] =     412.665187 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =         80.8  [cm-1]   --> B_1  D_3  S_3      
     freq (  2 -  2) =         80.8  [cm-1]   --> B_2  D_4  S_4      
     freq (  3 -  3) =        153.7  [cm-1]   --> A_1  D_1  S_1      
     freq (  4 -  4) =        355.8  [cm-1]   --> B_1  D_3  S_3      
     freq (  5 -  5) =        355.8  [cm-1]   --> B_2  D_4  S_4      
     freq (  6 -  6) =        412.7  [cm-1]   --> A_1  D_1  S_1      

     Calculation of q =    0.7500000  -0.2500000   0.7500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     48                 1221     1221     322
     Max         121     121     49                 1224     1224     323
     Sum         241     241     97                 2445     2445     645



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    40
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.0625000
        k(    2) = (   1.0000000   0.0000000   1.0000000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(    4) = (   1.0000000   0.0000000   1.5000000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0625000
        k(    6) = (   0.5000000   0.0000000   0.5000000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.1250000
        k(    8) = (   0.5000000  -0.5000000   1.0000000), wk =   0.0000000
        k(    9) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.0625000
        k(   10) = (   0.5000000  -0.5000000   0.5000000), wk =   0.0000000
        k(   11) = (   0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(   12) = (   1.0000000   0.0000000   0.5000000), wk =   0.0000000
        k(   13) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0625000
        k(   14) = (   1.0000000  -0.5000000   1.0000000), wk =   0.0000000
        k(   15) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   16) = (   0.5000000   0.0000000   0.0000000), wk =   0.0000000
        k(   17) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   18) = (   0.5000000  -0.5000000   1.5000000), wk =   0.0000000
        k(   19) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.1250000
        k(   20) = (   0.5000000  -1.0000000   1.0000000), wk =   0.0000000
        k(   21) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.1250000
        k(   22) = (   0.0000000  -0.5000000   1.0000000), wk =   0.0000000
        k(   23) = (  -0.2500000   0.7500000  -0.2500000), wk =   0.0625000
        k(   24) = (   0.5000000   0.5000000   0.5000000), wk =   0.0000000
        k(   25) = (   0.2500000   0.7500000   0.2500000), wk =   0.0625000
        k(   26) = (   1.0000000   0.5000000   1.0000000), wk =   0.0000000
        k(   27) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.1250000
        k(   28) = (   0.5000000  -0.5000000   0.0000000), wk =   0.0000000
        k(   29) = (   0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   30) = (   1.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   31) = (   0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   32) = (   1.0000000  -0.5000000   1.5000000), wk =   0.0000000
        k(   33) = (  -0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(   34) = (   0.5000000   0.0000000   1.5000000), wk =   0.0000000
        k(   35) = (  -0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   36) = (   0.5000000   0.5000000   1.0000000), wk =   0.0000000
        k(   37) = (  -0.2500000  -0.7500000  -0.2500000), wk =   0.0625000
        k(   38) = (   0.5000000  -1.0000000   0.5000000), wk =   0.0000000
        k(   39) = (   0.2500000  -0.7500000   0.2500000), wk =   0.0625000
        k(   40) = (   1.0000000  -1.0000000   1.0000000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.25E-10,  avg # of iterations = 11.0

     total cpu time spent up to now is        1.0 secs

     End of band structure calculation

          k = 0.2500 0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 1.0000 0.0000 1.0000 (   302 PWs)   bands (ev):

    -4.8217  -0.4470   2.9274   2.9274

          k = 0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000 0.0000 1.5000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k =-0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.5000 0.0000 0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.5000-0.5000 1.0000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 0.5000-0.5000 0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k = 0.2500 0.2500-0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 1.0000 0.0000 0.5000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k = 0.2500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -6.3575   1.7035   4.6970   4.6970

          k = 1.0000-0.5000 1.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.0000 0.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-0.5000 1.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k =-0.2500-0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-1.0000 1.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.7500-0.2500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.0000-0.5000 1.0000 (   308 PWs)   bands (ev):

    -4.7852  -0.0517   1.7949   2.1910

          k =-0.2500 0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.5000 0.5000 (   302 PWs)   bands (ev):

    -5.4218  -0.6403   4.3483   4.3483

          k = 0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000 0.5000 1.0000 (   311 PWs)   bands (ev):

    -6.1430   1.9396   3.7847   3.7847

          k =-0.2500-0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-0.5000 0.0000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k = 0.2500 0.2500-0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000 0.0000 0.0000 (   302 PWs)   bands (ev):

    -4.8217  -0.4470   2.9274   2.9274

          k = 0.2500-0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000-0.5000 1.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500 0.2500 0.7500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.0000 1.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500 0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000 0.5000 1.0000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k =-0.2500-0.7500-0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 0.5000-1.0000 0.5000 (   315 PWs)   bands (ev):

    -5.5287   0.5005   2.1485   4.2663

          k = 0.2500-0.7500 0.2500 (   311 PWs)   bands (ev):

    -5.1819  -0.0415   2.3125   3.5086

          k = 1.0000-1.0000 1.0000 (   331 PWs)   bands (ev):

    -6.9797   5.1761   5.1761   5.1761

     highest occupied level (ev):     4.6970

     Writing output data file alas.save

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.7500000  -0.2500000   0.7500000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    40

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  To be done

     Representation     2      1 modes -A'  To be done

     Representation     3      1 modes -A'  To be done

     Representation     4      1 modes -A'  To be done

     Representation     5      1 modes -A''  Not done in this run

     Representation     6      1 modes -A''  Not done in this run

     Compute atoms:     1,    2,



     Alpha used in Ewald sum =   0.7000
     PHONON       :     6.58s CPU         7.05s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     7.1 secs   av.it.:   6.3
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.089E-04

      iter #   2 total cpu time :     7.2 secs   av.it.:   8.7
      thresh= 1.044E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.311E-04

     Maximum CPU time exceeded

     max_seconds     =       7.00
     elapsed seconds =       7.09

     PHONON       :     6.74s CPU         7.24s WALL

     INITIALIZATION: 
     phq_setup    :      0.02s CPU      0.02s WALL (       5 calls)
     phq_init     :      0.17s CPU      0.17s WALL (       5 calls)

     phq_init     :      0.17s CPU      0.17s WALL (       5 calls)
     init_vloc    :      0.01s CPU      0.01s WALL (       5 calls)
     init_us_1    :      0.04s CPU      0.05s WALL (       5 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.06s CPU      0.06s WALL (       5 calls)
     phqscf       :      5.01s CPU      5.42s WALL (       5 calls)
     dynmatrix    :      0.01s CPU      0.01s WALL (       4 calls)

     phqscf       :      5.01s CPU      5.42s WALL (       6 calls)
     solve_linter :      4.94s CPU      5.34s WALL (      17 calls)
     drhodv       :      0.02s CPU      0.04s WALL (      16 calls)

     dynmat0      :      0.06s CPU      0.06s WALL (       5 calls)
     dynmat_us    :      0.04s CPU      0.04s WALL (       5 calls)
     d2ionq       :      0.01s CPU      0.01s WALL (       5 calls)

     dynmat_us    :      0.04s CPU      0.04s WALL (       5 calls)

     phqscf       :      5.01s CPU      5.42s WALL (       7 calls)
     solve_linter :      4.94s CPU      5.34s WALL (      18 calls)

     solve_linter :      4.94s CPU      5.34s WALL (      19 calls)
     dvqpsi_us    :      0.10s CPU      0.10s WALL (     206 calls)
     ortho        :      0.02s CPU      0.03s WALL (    1082 calls)
     cgsolve      :      3.60s CPU      3.95s WALL (    1082 calls)
     incdrhoscf   :      0.39s CPU      0.35s WALL (    1076 calls)
     vpsifft      :      0.27s CPU      0.29s WALL (     852 calls)
     dv_of_drho   :      0.02s CPU      0.03s WALL (     148 calls)
     mix_pot      :      0.04s CPU      0.05s WALL (     101 calls)
     psymdvscf    :      0.48s CPU      0.48s WALL (      96 calls)

     dvqpsi_us    :      0.10s CPU      0.10s WALL (     206 calls)
     dvqpsi_us_on :      0.02s CPU      0.01s WALL (     206 calls)

     cgsolve      :      3.60s CPU      3.95s WALL (    1082 calls)
     ch_psi       :      3.40s CPU      3.74s WALL (   10065 calls)

     ch_psi       :      3.40s CPU      3.74s WALL (   10065 calls)
     h_psi        :      3.32s CPU      3.66s WALL (   11309 calls)
     last         :      0.42s CPU      0.44s WALL (   10065 calls)

     h_psi        :      3.32s CPU      3.66s WALL (   11309 calls)
     add_vuspsi   :      0.20s CPU      0.18s WALL (   11309 calls)

     incdrhoscf   :      0.39s CPU      0.35s WALL (    1076 calls)


      General routines
     calbec       :      0.29s CPU      0.43s WALL (   22566 calls)
     fft          :      0.02s CPU      0.03s WALL (     464 calls)
     ffts         :      0.01s CPU      0.01s WALL (     302 calls)
     fftw         :      3.06s CPU      3.18s WALL (   98100 calls)
     davcio       :      0.04s CPU      0.06s WALL (    5899 calls)
     write_rec    :      0.17s CPU      0.17s WALL (     117 calls)


     PHONON       :     6.74s CPU         7.24s WALL


   This run was terminated on:  14:29: 4   7Feb2017            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
-------------------------------------------------------
Primary job  terminated normally, but 1 process returned
a non-zero exit code.. Per user-direction, the job has been aborted.
-------------------------------------------------------
--------------------------------------------------------------------------
mpirun detected that one or more processes exited with non-zero status, thus causing
the job to be terminated. The first process to do so was:

  Process name: [[33520,1],0]
  Exit code:    1
--------------------------------------------------------------------------
"""


test_61 = """
     Program PHONON v.6.0 (svn rev. 13286) starts on  7Feb2017 at 14: 1:20 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     2 processors
     R & G space division:  proc/nbgrp/npool/nimage =       2

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/alas.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     229
     Max         121     121     43                 1224     1224     230
     Sum         241     241     85                 2445     2445     459

        8 /   8 q-points for this run, from  1 to  8:
       N       xq(1)         xq(2)         xq(3) 
       1   0.000000000   0.000000000   0.000000000
       2  -0.250000000   0.250000000  -0.250000000
       3   0.500000000  -0.500000000   0.500000000
       4   0.000000000   0.500000000   0.000000000
       5   0.750000000  -0.250000000   0.750000000
       6   0.500000000   0.000000000   0.500000000
       7   0.000000000  -1.000000000   0.000000000
       8  -0.500000000  -1.000000000   0.000000000


     Calculation of q =    0.0000000   0.0000000   0.0000000

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   0.0000000 )

     25 Sym.Ops. (with q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=     2

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, T_d (-43m)  point group:


     Electric field:
     Dielectric constant
     Born effective charges in two ways 


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      3 modes -T_2  G_15 P_4  Done

     Representation     2      3 modes -T_2  G_15 P_4  Done


     PHONON       :     0.12s CPU         0.13s WALL


          Dielectric constant in cartesian axis 

          (      13.744216098       0.000000000       0.000000000 )
          (      -0.000000000      13.744216098      -0.000000000 )
          (       0.000000000       0.000000000      13.744216098 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   Al 
      Ex  (        1.88265       -0.00000       -0.00000 )
      Ey  (       -0.00000        1.88265        0.00000 )
      Ez  (       -0.00000        0.00000        1.88265 )
           atom      2   As 
      Ex  (       -3.23374       -0.00000       -0.00000 )
      Ey  (       -0.00000       -3.23374       -0.00000 )
      Ez  (       -0.00000       -0.00000       -3.23374 )

     Number of q in the star =    1
     List of q in the star:
          1   0.000000000   0.000000000   0.000000000

          Dielectric constant in cartesian axis 

          (      13.744216098       0.000000000       0.000000000 )
          (      -0.000000000      13.744216098      -0.000000000 )
          (       0.000000000       0.000000000      13.744216098 )

          Effective charges (d Force / dE) in cartesian axis

           atom      1   Al 
      Ex  (        1.88265       -0.00000       -0.00000 )
      Ey  (       -0.00000        1.88265        0.00000 )
      Ez  (       -0.00000        0.00000        1.88265 )
           atom      2   As 
      Ex  (       -3.23374       -0.00000       -0.00000 )
      Ey  (       -0.00000       -3.23374       -0.00000 )
      Ez  (       -0.00000       -0.00000       -3.23374 )

          Effective charges (d P / du) in cartesian axis 

           atom      1   Al 
      Px  (        1.88284       -0.00000        0.00000 )
      Py  (       -0.00000        1.88284        0.00000 )
      Pz  (        0.00000        0.00000        1.88284 )
           atom      2   As 
      Px  (       -3.23837       -0.00000       -0.00000 )
      Py  (       -0.00000       -3.23837        0.00000 )
      Pz  (       -0.00000        0.00000       -3.23837 )

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    2) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    3) =       0.164575 [THz] =       5.489636 [cm-1]
     freq (    4) =      11.258806 [THz] =     375.553348 [cm-1]
     freq (    5) =      11.258806 [THz] =     375.553348 [cm-1]
     freq (    6) =      11.258806 [THz] =     375.553348 [cm-1]
 **************************************************************************

     Mode symmetry, T_d (-43m)  point group:

     freq (  1 -  3) =          5.5  [cm-1]   --> T_2  G_15 P_4   I+R
     freq (  4 -  6) =        375.6  [cm-1]   --> T_2  G_15 P_4   I+R

     Calculation of q =   -0.2500000   0.2500000  -0.2500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     264
     Max         121     121     43                 1224     1224     267
     Sum         241     241     85                 2445     2445     531



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    20
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.1875000
        k(    2) = (   0.0000000   0.5000000   0.0000000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.3750000
        k(    4) = (   0.0000000   0.5000000   0.5000000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0625000
        k(    6) = (  -0.5000000   0.5000000  -0.5000000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.1875000
        k(    8) = (  -0.5000000   0.0000000  -0.5000000), wk =   0.0000000
        k(    9) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0625000
        k(   10) = (   0.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   11) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.1875000
        k(   12) = (  -0.5000000   0.5000000  -1.0000000), wk =   0.0000000
        k(   13) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.1875000
        k(   14) = (  -1.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   15) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.3750000
        k(   16) = (  -0.5000000   0.0000000  -1.0000000), wk =   0.0000000
        k(   17) = (   0.2500000  -0.2500000   0.7500000), wk =   0.1875000
        k(   18) = (   0.0000000   0.0000000   0.5000000), wk =   0.0000000
        k(   19) = (  -0.2500000   0.2500000   0.7500000), wk =   0.1875000
        k(   20) = (  -0.5000000   0.5000000   0.5000000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (  -0.2500000   0.2500000  -0.2500000 )

      6 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    20

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      1 modes -A_1  L_1  Done

     Representation     2      1 modes -A_1  L_1  Done

     Representation     3      2 modes -E    L_3  Done

     Representation     4      2 modes -E    L_3  Done


     PHONON       :     0.15s CPU         0.15s WALL


     Number of q in the star =    4
     List of q in the star:
          1  -0.250000000   0.250000000  -0.250000000
          2   0.250000000   0.250000000   0.250000000
          3   0.250000000  -0.250000000  -0.250000000
          4  -0.250000000  -0.250000000   0.250000000
     In addition there is the -q list: 
          1   0.250000000  -0.250000000   0.250000000
          2  -0.250000000  -0.250000000  -0.250000000
          3  -0.250000000   0.250000000   0.250000000
          4   0.250000000   0.250000000  -0.250000000

     Diagonalizing the dynamical matrix

     q = (   -0.250000000   0.250000000  -0.250000000 ) 

 **************************************************************************
     freq (    1) =       1.761219 [THz] =      58.747948 [cm-1]
     freq (    2) =       1.761219 [THz] =      58.747948 [cm-1]
     freq (    3) =       4.534095 [THz] =     151.241131 [cm-1]
     freq (    4) =      11.004845 [THz] =     367.082131 [cm-1]
     freq (    5) =      11.004845 [THz] =     367.082131 [cm-1]
     freq (    6) =      12.136604 [THz] =     404.833528 [cm-1]
 **************************************************************************

     Mode symmetry, C_3v (3m)   point group:

     freq (  1 -  2) =         58.7  [cm-1]   --> E    L_3           
     freq (  3 -  3) =        151.2  [cm-1]   --> A_1  L_1           
     freq (  4 -  5) =        367.1  [cm-1]   --> E    L_3           
     freq (  6 -  6) =        404.8  [cm-1]   --> A_1  L_1           

     Calculation of q =    0.5000000  -0.5000000   0.5000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     267
     Max         121     121     43                 1224     1224     270
     Sum         241     241     85                 2445     2445     537



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    10
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.3750000
        k(    2) = (   0.7500000  -0.2500000   0.7500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.7500000
        k(    4) = (   0.7500000  -0.2500000   1.2500000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(    6) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0000000
        k(    7) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.3750000
        k(    8) = (   0.2500000  -0.2500000  -0.2500000), wk =   0.0000000
        k(    9) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.3750000
        k(   10) = (  -0.2500000  -0.7500000   0.7500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.5000000  -0.5000000   0.5000000 )

      7 Sym.Ops. (with q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    10

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_3v (3m)   point group:


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      1 modes -A_1  L_1  Done

     Representation     2      1 modes -A_1  L_1  Done

     Representation     3      2 modes -E    L_3  Done

     Representation     4      2 modes -E    L_3  Done


     PHONON       :     0.16s CPU         0.18s WALL


     Number of q in the star =    4
     List of q in the star:
          1   0.500000000  -0.500000000   0.500000000
          2  -0.500000000  -0.500000000  -0.500000000
          3  -0.500000000   0.500000000   0.500000000
          4   0.500000000   0.500000000  -0.500000000

     Diagonalizing the dynamical matrix

     q = (    0.500000000  -0.500000000   0.500000000 ) 

 **************************************************************************
     freq (    1) =       2.016459 [THz] =      67.261839 [cm-1]
     freq (    2) =       2.016459 [THz] =      67.261839 [cm-1]
     freq (    3) =       6.494357 [THz] =     216.628431 [cm-1]
     freq (    4) =      10.940855 [THz] =     364.947641 [cm-1]
     freq (    5) =      10.940855 [THz] =     364.947641 [cm-1]
     freq (    6) =      11.551694 [THz] =     385.323027 [cm-1]
 **************************************************************************

     Mode symmetry, C_3v (3m)   point group:

     freq (  1 -  2) =         67.3  [cm-1]   --> E    L_3           
     freq (  3 -  3) =        216.6  [cm-1]   --> A_1  L_1           
     freq (  4 -  5) =        364.9  [cm-1]   --> E    L_3           
     freq (  6 -  6) =        385.3  [cm-1]   --> A_1  L_1           

     Calculation of q =    0.0000000   0.5000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     264
     Max         121     121     43                 1224     1224     267
     Sum         241     241     85                 2445     2445     531



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    24
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.1250000
        k(    2) = (   0.2500000   0.7500000   0.2500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.2500000
        k(    4) = (   0.2500000   0.7500000   0.7500000), wk =   0.0000000
        k(    5) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.1250000
        k(    6) = (  -0.2500000   0.2500000   0.2500000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.1250000
        k(    8) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0000000
        k(    9) = (   0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(   10) = (   0.2500000   0.7500000  -0.2500000), wk =   0.0000000
        k(   11) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.1250000
        k(   12) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.0000000
        k(   13) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.2500000
        k(   14) = (  -0.2500000   0.2500000   0.7500000), wk =   0.0000000
        k(   15) = (   0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   16) = (   0.2500000   1.2500000   0.2500000), wk =   0.0000000
        k(   17) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.2500000
        k(   18) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.0000000
        k(   19) = (   0.2500000   0.2500000  -0.7500000), wk =   0.2500000
        k(   20) = (   0.2500000   0.7500000  -0.7500000), wk =   0.0000000
        k(   21) = (  -0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   22) = (  -0.2500000   1.2500000   0.2500000), wk =   0.0000000
        k(   23) = (  -0.2500000  -0.7500000  -0.2500000), wk =   0.1250000
        k(   24) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.5000000   0.0000000 )

      4 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    24

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_2v (mm2)  point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A_1  D_1  S_1  Done

     Representation     2      1 modes -A_1  D_1  S_1  Done

     Representation     3      1 modes -B_1  D_3  S_3  Done

     Representation     4      1 modes -B_1  D_3  S_3  Done

     Representation     5      1 modes -B_2  D_4  S_4  Done

     Representation     6      1 modes -B_2  D_4  S_4  Done


     PHONON       :     0.19s CPU         0.20s WALL


     Number of q in the star =    6
     List of q in the star:
          1   0.000000000   0.500000000   0.000000000
          2  -0.500000000   0.000000000   0.000000000
          3   0.000000000   0.000000000  -0.500000000
          4   0.500000000   0.000000000   0.000000000
          5   0.000000000  -0.500000000   0.000000000
          6   0.000000000   0.000000000   0.500000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.500000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       2.421101 [THz] =      80.759234 [cm-1]
     freq (    2) =       2.421101 [THz] =      80.759234 [cm-1]
     freq (    3) =       4.606433 [THz] =     153.654058 [cm-1]
     freq (    4) =      10.666709 [THz] =     355.803105 [cm-1]
     freq (    5) =      10.666709 [THz] =     355.803105 [cm-1]
     freq (    6) =      12.371392 [THz] =     412.665230 [cm-1]
 **************************************************************************

     Mode symmetry, C_2v (mm2)  point group:

     freq (  1 -  1) =         80.8  [cm-1]   --> B_1  D_3  S_3      
     freq (  2 -  2) =         80.8  [cm-1]   --> B_2  D_4  S_4      
     freq (  3 -  3) =        153.7  [cm-1]   --> A_1  D_1  S_1      
     freq (  4 -  4) =        355.8  [cm-1]   --> B_1  D_3  S_3      
     freq (  5 -  5) =        355.8  [cm-1]   --> B_2  D_4  S_4      
     freq (  6 -  6) =        412.7  [cm-1]   --> A_1  D_1  S_1      

     Calculation of q =    0.7500000  -0.2500000   0.7500000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     48                 1221     1221     322
     Max         121     121     49                 1224     1224     323
     Sum         241     241     97                 2445     2445     645



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    40
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.0625000
        k(    2) = (   1.0000000   0.0000000   1.0000000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(    4) = (   1.0000000   0.0000000   1.5000000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0625000
        k(    6) = (   0.5000000   0.0000000   0.5000000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.1250000
        k(    8) = (   0.5000000  -0.5000000   1.0000000), wk =   0.0000000
        k(    9) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.0625000
        k(   10) = (   0.5000000  -0.5000000   0.5000000), wk =   0.0000000
        k(   11) = (   0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(   12) = (   1.0000000   0.0000000   0.5000000), wk =   0.0000000
        k(   13) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0625000
        k(   14) = (   1.0000000  -0.5000000   1.0000000), wk =   0.0000000
        k(   15) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   16) = (   0.5000000   0.0000000   0.0000000), wk =   0.0000000
        k(   17) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.1250000
        k(   18) = (   0.5000000  -1.0000000   1.0000000), wk =   0.0000000
        k(   19) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   20) = (   0.5000000  -0.5000000   1.5000000), wk =   0.0000000
        k(   21) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.1250000
        k(   22) = (   0.0000000  -0.5000000   1.0000000), wk =   0.0000000
        k(   23) = (   0.2500000   0.7500000   0.2500000), wk =   0.0625000
        k(   24) = (   1.0000000   0.5000000   1.0000000), wk =   0.0000000
        k(   25) = (  -0.2500000   0.7500000  -0.2500000), wk =   0.0625000
        k(   26) = (   0.5000000   0.5000000   0.5000000), wk =   0.0000000
        k(   27) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.1250000
        k(   28) = (   0.5000000  -0.5000000   0.0000000), wk =   0.0000000
        k(   29) = (   0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   30) = (   1.0000000   0.0000000   0.0000000), wk =   0.0000000
        k(   31) = (   0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   32) = (   1.0000000  -0.5000000   1.5000000), wk =   0.0000000
        k(   33) = (  -0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(   34) = (   0.5000000   0.0000000   1.5000000), wk =   0.0000000
        k(   35) = (  -0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   36) = (   0.5000000   0.5000000   1.0000000), wk =   0.0000000
        k(   37) = (  -0.2500000  -0.7500000  -0.2500000), wk =   0.0625000
        k(   38) = (   0.5000000  -1.0000000   0.5000000), wk =   0.0000000
        k(   39) = (   0.2500000  -0.7500000   0.2500000), wk =   0.0625000
        k(   40) = (   1.0000000  -1.0000000   1.0000000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.7500000  -0.2500000   0.7500000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    40

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  Done

     Representation     2      1 modes -A'  Done

     Representation     3      1 modes -A'  Done

     Representation     4      1 modes -A'  Done

     Representation     5      1 modes -A''  Done

     Representation     6      1 modes -A''  Done


     PHONON       :     0.22s CPU         0.23s WALL


     Number of q in the star =   12
     List of q in the star:
          1   0.750000000  -0.250000000   0.750000000
          2  -0.750000000  -0.250000000  -0.750000000
          3   0.250000000  -0.750000000   0.750000000
          4   0.750000000  -0.750000000   0.250000000
          5  -0.250000000  -0.750000000  -0.750000000
          6  -0.750000000   0.250000000   0.750000000
          7   0.750000000   0.750000000  -0.250000000
          8  -0.750000000  -0.750000000  -0.250000000
          9  -0.750000000   0.750000000   0.250000000
         10   0.750000000   0.250000000  -0.750000000
         11  -0.250000000   0.750000000   0.750000000
         12   0.250000000   0.750000000  -0.750000000
     In addition there is the -q list: 
          1  -0.750000000   0.250000000  -0.750000000
          2   0.750000000   0.250000000   0.750000000
          3  -0.250000000   0.750000000  -0.750000000
          4  -0.750000000   0.750000000  -0.250000000
          5   0.250000000   0.750000000   0.750000000
          6   0.750000000  -0.250000000  -0.750000000
          7  -0.750000000  -0.750000000   0.250000000
          8   0.750000000   0.750000000   0.250000000
          9   0.750000000  -0.750000000  -0.250000000
         10  -0.750000000  -0.250000000   0.750000000
         11   0.250000000  -0.750000000  -0.750000000
         12  -0.250000000  -0.750000000   0.750000000

     Diagonalizing the dynamical matrix

     q = (    0.750000000  -0.250000000   0.750000000 ) 

 **************************************************************************
     freq (    1) =       2.620989 [THz] =      87.426770 [cm-1]
     freq (    2) =       3.804579 [THz] =     126.907082 [cm-1]
     freq (    3) =       5.902630 [THz] =     196.890533 [cm-1]
     freq (    4) =      10.568988 [THz] =     352.543493 [cm-1]
     freq (    5) =      10.588928 [THz] =     353.208622 [cm-1]
     freq (    6) =      11.478004 [THz] =     382.865013 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =         87.4  [cm-1]   --> A''                
     freq (  2 -  2) =        126.9  [cm-1]   --> A'                 
     freq (  3 -  3) =        196.9  [cm-1]   --> A'                 
     freq (  4 -  4) =        352.5  [cm-1]   --> A''                
     freq (  5 -  5) =        353.2  [cm-1]   --> A'                 
     freq (  6 -  6) =        382.9  [cm-1]   --> A'                 

     Calculation of q =    0.5000000   0.0000000   0.5000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     42                 1221     1221     267
     Max         121     121     43                 1224     1224     270
     Sum         241     241     85                 2445     2445     537



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    40
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.0625000
        k(    2) = (   0.7500000   0.2500000   0.7500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(    4) = (   0.7500000   0.2500000   1.2500000), wk =   0.0000000
        k(    5) = (  -0.2500000   0.2500000  -0.2500000), wk =   0.0625000
        k(    6) = (   0.2500000   0.2500000   0.2500000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000   0.2500000), wk =   0.1250000
        k(    8) = (   0.2500000  -0.2500000   0.7500000), wk =   0.0000000
        k(    9) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.0625000
        k(   10) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0000000
        k(   11) = (   0.2500000   0.2500000  -0.2500000), wk =   0.1250000
        k(   12) = (   0.7500000   0.2500000   0.2500000), wk =   0.0000000
        k(   13) = (   0.2500000  -0.2500000   0.2500000), wk =   0.0625000
        k(   14) = (   0.7500000  -0.2500000   0.7500000), wk =   0.0000000
        k(   15) = (  -0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   16) = (   0.2500000   0.2500000  -0.2500000), wk =   0.0000000
        k(   17) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.1250000
        k(   18) = (   0.2500000  -0.7500000   0.7500000), wk =   0.0000000
        k(   19) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   20) = (   0.2500000  -0.2500000   1.2500000), wk =   0.0000000
        k(   21) = (  -0.7500000  -0.2500000   0.2500000), wk =   0.1250000
        k(   22) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.0000000
        k(   23) = (   0.2500000   0.7500000   0.2500000), wk =   0.0625000
        k(   24) = (   0.7500000   0.7500000   0.7500000), wk =   0.0000000
        k(   25) = (  -0.2500000   0.7500000  -0.2500000), wk =   0.0625000
        k(   26) = (   0.2500000   0.7500000   0.2500000), wk =   0.0000000
        k(   27) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.1250000
        k(   28) = (   0.2500000  -0.2500000  -0.2500000), wk =   0.0000000
        k(   29) = (   0.2500000   0.2500000  -0.7500000), wk =   0.1250000
        k(   30) = (   0.7500000   0.2500000  -0.2500000), wk =   0.0000000
        k(   31) = (   0.2500000  -0.2500000   0.7500000), wk =   0.1250000
        k(   32) = (   0.7500000  -0.2500000   1.2500000), wk =   0.0000000
        k(   33) = (  -0.2500000   0.2500000   0.7500000), wk =   0.1250000
        k(   34) = (   0.2500000   0.2500000   1.2500000), wk =   0.0000000
        k(   35) = (  -0.2500000   0.7500000   0.2500000), wk =   0.1250000
        k(   36) = (   0.2500000   0.7500000   0.7500000), wk =   0.0000000
        k(   37) = (  -0.2500000  -0.7500000  -0.2500000), wk =   0.0625000
        k(   38) = (   0.2500000  -0.7500000   0.2500000), wk =   0.0000000
        k(   39) = (   0.2500000  -0.7500000   0.2500000), wk =   0.0625000
        k(   40) = (   0.7500000  -0.7500000   0.7500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.5000000   0.0000000   0.5000000 )

      2 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    40

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, C_s (m)     point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A'  Done

     Representation     2      1 modes -A'  Done

     Representation     3      1 modes -A'  Done

     Representation     4      1 modes -A'  Done

     Representation     5      1 modes -A''  Done

     Representation     6      1 modes -A''  Done


     PHONON       :     0.24s CPU         0.26s WALL


     Number of q in the star =   12
     List of q in the star:
          1   0.500000000   0.000000000   0.500000000
          2  -0.500000000   0.000000000  -0.500000000
          3   0.000000000  -0.500000000   0.500000000
          4   0.500000000  -0.500000000   0.000000000
          5   0.000000000  -0.500000000  -0.500000000
          6  -0.500000000   0.000000000   0.500000000
          7   0.500000000   0.500000000   0.000000000
          8  -0.500000000  -0.500000000   0.000000000
          9  -0.500000000   0.500000000   0.000000000
         10   0.500000000   0.000000000  -0.500000000
         11   0.000000000   0.500000000   0.500000000
         12   0.000000000   0.500000000  -0.500000000

     Diagonalizing the dynamical matrix

     q = (    0.500000000   0.000000000   0.500000000 ) 

 **************************************************************************
     freq (    1) =       2.514962 [THz] =      83.890086 [cm-1]
     freq (    2) =       3.826519 [THz] =     127.638919 [cm-1]
     freq (    3) =       5.424240 [THz] =     180.933178 [cm-1]
     freq (    4) =      10.719166 [THz] =     357.552875 [cm-1]
     freq (    5) =      10.737649 [THz] =     358.169410 [cm-1]
     freq (    6) =      11.302704 [THz] =     377.017617 [cm-1]
 **************************************************************************

     Mode symmetry, C_s (m)     point group:

     freq (  1 -  1) =         83.9  [cm-1]   --> A''                
     freq (  2 -  2) =        127.6  [cm-1]   --> A'                 
     freq (  3 -  3) =        180.9  [cm-1]   --> A'                 
     freq (  4 -  4) =        357.6  [cm-1]   --> A'                 
     freq (  5 -  5) =        358.2  [cm-1]   --> A''                
     freq (  6 -  6) =        377.0  [cm-1]   --> A'                 

     Calculation of q =    0.0000000  -1.0000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     45                 1221     1221     304
     Max         121     121     46                 1224     1224     305
     Sum         241     241     91                 2445     2445     609



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=     6
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.5000000
        k(    2) = (   0.2500000  -0.7500000   0.2500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   1.0000000
        k(    4) = (   0.2500000  -0.7500000   0.7500000), wk =   0.0000000
        k(    5) = (   0.2500000  -0.7500000  -0.2500000), wk =   0.5000000
        k(    6) = (   0.2500000  -1.7500000  -0.2500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (   0.0000000  -1.0000000   0.0000000 )

      9 Sym.Ops. (with q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=     6

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, D_2d (-42m) point group:


     Atomic displacements:
     There are   4 irreducible representations

     Representation     1      1 modes -A_1  X_1  W_1  Done

     Representation     2      1 modes -B_2  X_3  W_2  Done

     Representation     3      2 modes -E    X_5  W_3  Done

     Representation     4      2 modes -E    X_5  W_3  Done


     PHONON       :     0.26s CPU         0.28s WALL


     Number of q in the star =    3
     List of q in the star:
          1   0.000000000  -1.000000000   0.000000000
          2   0.000000000   0.000000000  -1.000000000
          3  -1.000000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000  -1.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       2.844774 [THz] =      94.891459 [cm-1]
     freq (    2) =       2.844774 [THz] =      94.891459 [cm-1]
     freq (    3) =       6.564967 [THz] =     218.983735 [cm-1]
     freq (    4) =      10.442955 [THz] =     348.339475 [cm-1]
     freq (    5) =      10.442955 [THz] =     348.339475 [cm-1]
     freq (    6) =      12.210582 [THz] =     407.301190 [cm-1]
 **************************************************************************

     Mode symmetry, D_2d (-42m) point group:

     freq (  1 -  2) =         94.9  [cm-1]   --> E    X_5  W_3      
     freq (  3 -  3) =        219.0  [cm-1]   --> A_1  X_1  W_1      
     freq (  4 -  5) =        348.3  [cm-1]   --> E    X_5  W_3      
     freq (  6 -  6) =        407.3  [cm-1]   --> B_2  X_3  W_2      

     Calculation of q =   -0.5000000  -1.0000000   0.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         120     120     45                 1221     1221     304
     Max         121     121     46                 1224     1224     305
     Sum         241     241     91                 2445     2445     609



     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      16.0000  Ry
     charge density cutoff     =      64.0000  Ry
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)

     celldm(1)=  10.500000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     atomic species   valence    mass     pseudopotential
        Al             3.00    26.98000     Al( 1.00)
        As             5.00    74.92000     As( 1.00)

     24 Sym. Ops. (no inversion) found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Al  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           As  tau(   2) = (   0.2500000   0.2500000   0.2500000  )

     number of k points=    16
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.2500000   0.2500000   0.2500000), wk =   0.2500000
        k(    2) = (  -0.2500000  -0.7500000   0.2500000), wk =   0.0000000
        k(    3) = (   0.2500000   0.2500000   0.7500000), wk =   0.2500000
        k(    4) = (  -0.2500000  -0.7500000   0.7500000), wk =   0.0000000
        k(    5) = (  -0.2500000  -0.2500000  -0.2500000), wk =   0.2500000
        k(    6) = (  -0.7500000  -1.2500000  -0.2500000), wk =   0.0000000
        k(    7) = (  -0.2500000  -0.2500000   0.7500000), wk =   0.2500000
        k(    8) = (  -0.7500000  -1.2500000   0.7500000), wk =   0.0000000
        k(    9) = (   0.7500000   0.2500000   0.2500000), wk =   0.2500000
        k(   10) = (   0.2500000  -0.7500000   0.2500000), wk =   0.0000000
        k(   11) = (  -0.2500000  -0.2500000  -0.7500000), wk =   0.2500000
        k(   12) = (  -0.7500000  -1.2500000  -0.7500000), wk =   0.0000000
        k(   13) = (   0.2500000   0.2500000  -0.7500000), wk =   0.2500000
        k(   14) = (  -0.2500000  -0.7500000  -0.7500000), wk =   0.0000000
        k(   15) = (   0.7500000  -0.2500000   0.2500000), wk =   0.2500000
        k(   16) = (   0.2500000  -1.2500000   0.2500000), wk =   0.0000000

     Dense  grid:     2445 G-vectors     FFT dimensions: (  20,  20,  20)

     Estimated max dynamical RAM per process >       0.51MB

     Estimated total allocated dynamical RAM >       1.02MB

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/alas.save/charge-density.dat

     Starting wfc are    8 atomic wfcs

                                                                                

     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.5000  a.u.
     unit-cell volume          =     289.4062 (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            2
     kinetic-energy cut-off    =      16.0000  Ry
     charge density cut-off    =      64.0000  Ry
     convergence threshold     =      1.0E-12
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PZ   NOGX NOGC ( 1  1  0  0 0 0)


     celldm(1)=   10.50000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Al  26.9800   tau(    1) = (    0.00000    0.00000    0.00000  )
        2     As  74.9200   tau(    2) = (    0.25000    0.25000    0.25000  )

     Computing dynamical matrix for 
                    q = (  -0.5000000  -1.0000000   0.0000000 )

      4 Sym.Ops. (no q -> -q+G )


     G cutoff =  178.7306  (   1224 G-vectors)     FFT grid: ( 20, 20, 20)
     number of k points=    16

     PseudoPot. # 1 for Al read from file:
     /home/pietro/espresso-svn/pseudo/Al.pz-vbc.UPF
     MD5 check sum: 614279c88ff8d45c90147292d03ed420
     Pseudo is Norm-conserving, Zval =  3.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  171 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     PseudoPot. # 2 for As read from file:
     /home/pietro/espresso-svn/pseudo/As.pz-bhs.UPF
     MD5 check sum: 451cd3365afcfc94d28b1934951c34a8
     Pseudo is Norm-conserving, Zval =  5.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of  525 points,  2 beta functions with: 
                l(1) =   0
                l(2) =   1

     Mode symmetry, S_4 (-4)    point group:


     Atomic displacements:
     There are   6 irreducible representations

     Representation     1      1 modes -A    W_1  Done

     Representation     2      1 modes -B    W_3  Done

     Representation     3      1 modes -B    W_3  Done

     Representation     4      1 modes -E    W_4  Done

     Representation     5      1 modes -E    W_4  Done

     Representation     6      1 modes -E*   W_2  Done


     PHONON       :     0.28s CPU         0.31s WALL


     Number of q in the star =    6
     List of q in the star:
          1  -0.500000000  -1.000000000   0.000000000
          2   0.500000000   1.000000000   0.000000000
          3   0.000000000  -1.000000000  -0.500000000
          4   0.000000000   1.000000000   0.500000000
          5   0.000000000  -0.500000000  -1.000000000
          6   0.000000000   0.500000000   1.000000000

     Diagonalizing the dynamical matrix

     q = (   -0.500000000  -1.000000000   0.000000000 ) 

 **************************************************************************
     freq (    1) =       3.747049 [THz] =     124.988106 [cm-1]
     freq (    2) =       4.016743 [THz] =     133.984119 [cm-1]
     freq (    3) =       5.965593 [THz] =     198.990748 [cm-1]
     freq (    4) =      10.537211 [THz] =     351.483535 [cm-1]
     freq (    5) =      10.644715 [THz] =     355.069463 [cm-1]
     freq (    6) =      10.758903 [THz] =     358.878368 [cm-1]
 **************************************************************************

     Mode symmetry, S_4 (-4)    point group:

     freq (  1 -  1) =        125.0  [cm-1]   --> B    W_3           
     freq (  2 -  2) =        134.0  [cm-1]   --> E    W_4           
     freq (  3 -  3) =        199.0  [cm-1]   --> A    W_1           
     freq (  4 -  4) =        351.5  [cm-1]   --> B    W_3           
     freq (  5 -  5) =        355.1  [cm-1]   --> E*   W_2           
     freq (  6 -  6) =        358.9  [cm-1]   --> E    W_4           

     init_run     :      0.12s CPU      0.12s WALL (       7 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       7 calls)
     potinit      :      0.01s CPU      0.01s WALL (       7 calls)

     Called by electrons:
     v_of_rho     :      0.00s CPU      0.00s WALL (       8 calls)

     Called by c_bands:

     Called by sum_band:

     Called by *egterg:

     Called by h_psi:

     General routines
     fft          :      0.00s CPU      0.00s WALL (      24 calls)

     Parallel routines
     fft_scatter  :      0.00s CPU      0.00s WALL (      24 calls)

     PHONON       :     0.28s CPU         0.31s WALL

     INITIALIZATION: 
     phq_setup    :      0.02s CPU      0.02s WALL (       8 calls)

     init_vloc    :      0.00s CPU      0.01s WALL (       8 calls)
     init_us_1    :      0.07s CPU      0.06s WALL (       8 calls)

     DYNAMICAL MATRIX:
     phqscf       :      0.00s CPU      0.00s WALL (       8 calls)
     dynmatrix    :      0.01s CPU      0.01s WALL (       8 calls)

     phqscf       :      0.00s CPU      0.00s WALL (       8 calls)



     phqscf       :      0.00s CPU      0.00s WALL (       8 calls)








      General routines
     fft          :      0.00s CPU      0.00s WALL (      24 calls)


     PHONON       :     0.28s CPU         0.31s WALL


   This run was terminated on:  14: 1:20   7Feb2017            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


test_66 = """
     Program PHONON v.6.0 (svn rev. 13188M) starts on  7Dec2016 at  0: 1:34 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     4 processors
     R & G space division:  proc/nbgrp/npool/nimage =       4

     Reading data from directory:
     /home/pietro/espresso-svn/tempdir/nickel.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      =  SLA  PW   PBE  PBE ( 1  4  3  4 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want

               file Ni.pbe-nd-rrkjus.UPF: wavefunction(s)  4S renormalized

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         112      40     15                 1604      351      82
     Max         113      41     16                 1607      354      83
     Sum         451     163     61                 6423     1411     331

     Generating pointlists ...
     new r_m :   0.2917 (alat units)  1.9397 (a.u.) for type    1

     Check: negative/imaginary core charge=   -0.000021    0.000000

     Calculation of q =    0.0000000   0.0000000   1.0000000

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         112      40     21                 1604      351     132
     Max         113      41     22                 1607      354     135
     Sum         451     163     85                 6423     1411     531


     Title: 
     phonons of Ni at X                                                         


     bravais-lattice index     =            2
     lattice parameter (alat)  =       6.6500  a.u.
     unit-cell volume          =      73.5199 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     number of electrons       =        10.00
     number of Kohn-Sham states=            9
     kinetic-energy cutoff     =      27.0000  Ry
     charge density cutoff     =     300.0000  Ry
     Exchange-correlation      =  SLA  PW   PBE  PBE ( 1  4  3  4 0 0)

     celldm(1)=   6.650000  celldm(2)=   0.000000  celldm(3)=   0.000000
     celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000

     crystal axes: (cart. coord. in units of alat)
               a(1) = (  -0.500000   0.000000   0.500000 )  
               a(2) = (   0.000000   0.500000   0.500000 )  
               a(3) = (  -0.500000   0.500000   0.000000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.000000 -1.000000  1.000000 )  
               b(2) = (  1.000000  1.000000  1.000000 )  
               b(3) = ( -1.000000  1.000000 -1.000000 )  


     PseudoPot. # 1 for Ni read from file:
     ./Ni.pbe-nd-rrkjus.UPF
     MD5 check sum: 8081f0a005c9a5470caab1a58e82ecb2
     Pseudo is Ultrasoft + core correction, Zval = 10.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of 1203 points,  6 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
                l(5) =   2
                l(6) =   2
     Q(r) pseudized with 0 coefficients 


     atomic species   valence    mass     pseudopotential
        Ni            10.00    58.69340     Ni( 1.00)

     Starting magnetic structure 
     atomic species   magnetization
        Ni           0.500

     48 Sym. Ops., with inversion, found



   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Ni  tau(   1) = (   0.0000000   0.0000000   0.0000000  )

     number of k points=    40  Marzari-Vanderbilt smearing, width (Ry)=  0.0200
                       cart. coord. in units 2pi/alat
        k(    1) = (  -0.1250000   0.1250000   0.1250000), wk =   0.0312500
        k(    2) = (  -0.1250000   0.1250000   1.1250000), wk =   0.0000000
        k(    3) = (  -0.3750000   0.3750000  -0.1250000), wk =   0.0312500
        k(    4) = (  -0.3750000   0.3750000   0.8750000), wk =   0.0000000
        k(    5) = (   0.3750000  -0.3750000   0.6250000), wk =   0.0312500
        k(    6) = (   0.3750000  -0.3750000   1.6250000), wk =   0.0000000
        k(    7) = (   0.1250000  -0.1250000   0.3750000), wk =   0.0312500
        k(    8) = (   0.1250000  -0.1250000   1.3750000), wk =   0.0000000
        k(    9) = (  -0.1250000   0.6250000   0.1250000), wk =   0.0625000
        k(   10) = (  -0.1250000   0.6250000   1.1250000), wk =   0.0000000
        k(   11) = (   0.6250000  -0.1250000   0.8750000), wk =   0.0625000
        k(   12) = (   0.6250000  -0.1250000   1.8750000), wk =   0.0000000
        k(   13) = (   0.3750000   0.1250000   0.6250000), wk =   0.0625000
        k(   14) = (   0.3750000   0.1250000   1.6250000), wk =   0.0000000
        k(   15) = (  -0.1250000  -0.8750000   0.1250000), wk =   0.0625000
        k(   16) = (  -0.1250000  -0.8750000   1.1250000), wk =   0.0000000
        k(   17) = (  -0.3750000   0.3750000   0.3750000), wk =   0.0312500
        k(   18) = (  -0.3750000   0.3750000   1.3750000), wk =   0.0000000
        k(   19) = (   0.3750000  -0.3750000   1.1250000), wk =   0.0312500
        k(   20) = (   0.3750000  -0.3750000   2.1250000), wk =   0.0000000
        k(   21) = (   0.3750000  -0.1250000  -0.3750000), wk =   0.0625000
        k(   22) = (   0.3750000  -0.1250000   0.6250000), wk =   0.0000000
        k(   23) = (  -0.3750000   0.6250000   0.3750000), wk =   0.0625000
        k(   24) = (  -0.3750000   0.6250000   1.3750000), wk =   0.0000000
        k(   25) = (  -0.1250000   0.3750000   0.1250000), wk =   0.0625000
        k(   26) = (  -0.1250000   0.3750000   1.1250000), wk =   0.0000000
        k(   27) = (   0.1250000  -0.1250000   0.6250000), wk =   0.0312500
        k(   28) = (   0.1250000  -0.1250000   1.6250000), wk =   0.0000000
        k(   29) = (  -0.1250000   0.8750000   0.6250000), wk =   0.0625000
        k(   30) = (  -0.1250000   0.8750000   1.6250000), wk =   0.0000000
        k(   31) = (   0.8750000   0.6250000  -0.1250000), wk =   0.0625000
        k(   32) = (   0.8750000   0.6250000   0.8750000), wk =   0.0000000
        k(   33) = (   0.1250000   0.6250000   0.3750000), wk =   0.0625000
        k(   34) = (   0.1250000   0.6250000   1.3750000), wk =   0.0000000
        k(   35) = (   0.6250000   0.3750000   0.1250000), wk =   0.0625000
        k(   36) = (   0.6250000   0.3750000   1.1250000), wk =   0.0000000
        k(   37) = (   0.1250000  -0.1250000  -0.8750000), wk =   0.0312500
        k(   38) = (   0.1250000  -0.1250000   0.1250000), wk =   0.0000000
        k(   39) = (  -0.3750000   1.1250000   0.3750000), wk =   0.0625000
        k(   40) = (  -0.3750000   1.1250000   1.3750000), wk =   0.0000000

     Dense  grid:     6423 G-vectors     FFT dimensions: (  25,  25,  25)

     Smooth grid:     1411 G-vectors     FFT dimensions: (  15,  15,  15)

     Estimated max dynamical RAM per process >       0.80Mb

     Estimated total allocated dynamical RAM >       3.19Mb
     Generating pointlists ...

     Check: negative/imaginary core charge=   -0.000021    0.000000

     The potential is recalculated from file :
     /home/pietro/espresso-svn/tempdir/_ph0/nickel.save/charge-density.dat

     Starting wfc are    6 atomic +    3 random wfc

     Band Structure Calculation
     Davidson diagonalization with overlap

     ethr =  1.00E-10,  avg # of iterations = 14.1

     total cpu time spent up to now is        1.1 secs

     End of band structure calculation

 ------ SPIN UP ------------


          k =-0.1250 0.1250 0.1250 (   172 PWs)   bands (ev):

     5.8697  11.5739  11.8319  11.8319  12.8613  12.8613  35.2151  39.1169
    41.0563

          k =-0.1250 0.1250 1.1250 (   176 PWs)   bands (ev):

     9.7774  10.1644  12.8693  13.3037  13.6224  16.7880  24.9790  26.3754
    30.0885

          k =-0.3750 0.3750-0.1250 (   171 PWs)   bands (ev):

     8.5750  11.2500  11.8343  12.1287  12.7521  13.6728  27.1042  32.6458
    39.6758

          k =-0.3750 0.3750 0.8750 (   176 PWs)   bands (ev):

    10.3654  11.0168  11.5578  12.5018  13.2687  17.7551  21.2364  27.2375
    34.3327

          k = 0.3750-0.3750 0.6250 (   172 PWs)   bands (ev):

     9.6624  11.5167  11.9827  12.1971  13.5534  15.4848  20.4983  33.7468
    36.0279

          k = 0.3750-0.3750 1.6250 (   174 PWs)   bands (ev):

     9.0450  11.8255  11.8255  12.3363  13.3397  13.3397  23.0016  37.0669
    39.2790

          k = 0.1250-0.1250 0.3750 (   169 PWs)   bands (ev):

     7.3631  11.1757  12.0273  12.1377  12.6935  13.1369  31.2698  36.2535
    36.8262

          k = 0.1250-0.1250 1.3750 (   178 PWs)   bands (ev):

     9.3859  10.5803  12.0474  12.7104  13.4797  13.7867  28.1568  31.5072
    32.3294

          k =-0.1250 0.6250 0.1250 (   178 PWs)   bands (ev):

     9.3859  10.5803  12.0474  12.7104  13.4797  13.7867  28.1568  31.5072
    32.3294

          k =-0.1250 0.6250 1.1250 (   179 PWs)   bands (ev):

    10.3860  10.6412  11.6238  12.9147  13.5152  19.0390  22.3265  26.0101
    28.3110

          k = 0.6250-0.1250 0.8750 (   179 PWs)   bands (ev):

    10.3860  10.6412  11.6238  12.9147  13.5152  19.0390  22.3265  26.0101
    28.3110

          k = 0.6250-0.1250 1.8750 (   178 PWs)   bands (ev):

     9.3859  10.5803  12.0474  12.7104  13.4797  13.7867  28.1568  31.5072
    32.3294

          k = 0.3750 0.1250 0.6250 (   174 PWs)   bands (ev):

    10.0139  11.0553  11.4269  12.4907  13.2324  15.3091  24.0931  29.7562
    32.8980

          k = 0.3750 0.1250 1.6250 (   171 PWs)   bands (ev):

     8.5750  11.2500  11.8343  12.1287  12.7521  13.6728  27.1042  32.6458
    39.6758

          k =-0.1250-0.8750 0.1250 (   176 PWs)   bands (ev):

     9.7774  10.1644  12.8693  13.3037  13.6224  16.7880  24.9790  26.3754
    30.0885

          k =-0.1250-0.8750 1.1250 (   176 PWs)   bands (ev):

     9.7774  10.1644  12.8693  13.3037  13.6224  16.7880  24.9790  26.3754
    30.0885

          k =-0.3750 0.3750 0.3750 (   174 PWs)   bands (ev):

     9.0450  11.8255  11.8255  12.3363  13.3397  13.3397  23.0016  37.0669
    39.2790

          k =-0.3750 0.3750 1.3750 (   172 PWs)   bands (ev):

     9.6624  11.5167  11.9827  12.1971  13.5534  15.4848  20.4983  33.7468
    36.0279

          k = 0.3750-0.3750 1.1250 (   176 PWs)   bands (ev):

    10.3654  11.0168  11.5578  12.5018  13.2687  17.7551  21.2364  27.2375
    34.3327

          k = 0.3750-0.3750 2.1250 (   171 PWs)   bands (ev):

     8.5750  11.2500  11.8343  12.1287  12.7521  13.6728  27.1042  32.6458
    39.6758

          k = 0.3750-0.1250-0.3750 (   171 PWs)   bands (ev):

     8.5750  11.2500  11.8343  12.1287  12.7521  13.6728  27.1042  32.6458
    39.6758

          k = 0.3750-0.1250 0.6250 (   174 PWs)   bands (ev):

    10.0139  11.0553  11.4269  12.4907  13.2324  15.3091  24.0931  29.7562
    32.8980

          k =-0.3750 0.6250 0.3750 (   172 PWs)   bands (ev):

     9.6624  11.5167  11.9827  12.1971  13.5534  15.4848  20.4983  33.7468
    36.0279

          k =-0.3750 0.6250 1.3750 (   172 PWs)   bands (ev):

     9.6624  11.5167  11.9827  12.1971  13.5534  15.4848  20.4983  33.7468
    36.0279

          k =-0.1250 0.3750 0.1250 (   169 PWs)   bands (ev):

     7.3631  11.1757  12.0273  12.1377  12.6935  13.1369  31.2698  36.2535
    36.8262

          k =-0.1250 0.3750 1.1250 (   179 PWs)   bands (ev):

    10.3860  10.6412  11.6238  12.9147  13.5152  19.0390  22.3265  26.0101
    28.3110

          k = 0.1250-0.1250 0.6250 (   178 PWs)   bands (ev):

     9.3859  10.5803  12.0474  12.7104  13.4797  13.7867  28.1568  31.5072
    32.3294

          k = 0.1250-0.1250 1.6250 (   169 PWs)   bands (ev):

     7.3631  11.1757  12.0273  12.1377  12.6935  13.1369  31.2698  36.2535
    36.8262

          k =-0.1250 0.8750 0.6250 (   179 PWs)   bands (ev):

    10.3860  10.6412  11.6238  12.9147  13.5152  19.0390  22.3265  26.0101
    28.3110

          k =-0.1250 0.8750 1.6250 (   179 PWs)   bands (ev):

    10.3860  10.6412  11.6238  12.9147  13.5152  19.0390  22.3265  26.0101
    28.3110

          k = 0.8750 0.6250-0.1250 (   179 PWs)   bands (ev):

    10.3860  10.6412  11.6238  12.9147  13.5152  19.0390  22.3265  26.0101
    28.3110

          k = 0.8750 0.6250 0.8750 (   169 PWs)   bands (ev):

     7.3631  11.1757  12.0273  12.1377  12.6935  13.1369  31.2698  36.2535
    36.8262

          k = 0.1250 0.6250 0.3750 (   174 PWs)   bands (ev):

    10.0139  11.0553  11.4269  12.4907  13.2324  15.3091  24.0931  29.7562
    32.8980

          k = 0.1250 0.6250 1.3750 (   176 PWs)   bands (ev):

    10.3654  11.0168  11.5578  12.5018  13.2687  17.7551  21.2364  27.2375
    34.3327

          k = 0.6250 0.3750 0.1250 (   174 PWs)   bands (ev):

    10.0139  11.0553  11.4269  12.4907  13.2324  15.3091  24.0931  29.7562
    32.8980

          k = 0.6250 0.3750 1.1250 (   174 PWs)   bands (ev):

    10.0139  11.0553  11.4269  12.4907  13.2324  15.3091  24.0931  29.7562
    32.8980

          k = 0.1250-0.1250-0.8750 (   176 PWs)   bands (ev):

     9.7774  10.1644  12.8693  13.3037  13.6224  16.7880  24.9790  26.3754
    30.0885

          k = 0.1250-0.1250 0.1250 (   172 PWs)   bands (ev):

     5.8697  11.5739  11.8319  11.8319  12.8613  12.8613  35.2151  39.1169
    41.0563

          k =-0.3750 1.1250 0.3750 (   176 PWs)   bands (ev):

    10.3654  11.0168  11.5578  12.5018  13.2687  17.7551  21.2364  27.2375
    34.3327

          k =-0.3750 1.1250 1.3750 (   174 PWs)   bands (ev):

    10.0139  11.0553  11.4269  12.4907  13.2324  15.3091  24.0931  29.7562
    32.8980

 ------ SPIN DOWN ----------


          k =-0.1250 0.1250 0.1250 (   172 PWs)   bands (ev):

     5.8234  12.4452  12.7307  12.7307  13.5994  13.5994  35.2386  38.9838
    41.0911

          k =-0.1250 0.1250 1.1250 (   176 PWs)   bands (ev):

    10.2090  10.8957  13.6527  14.1097  14.5846  17.0385  25.1835  26.4722
    30.1022

          k =-0.3750 0.3750-0.1250 (   171 PWs)   bands (ev):

     8.6207  11.9920  12.5953  12.9299  13.5959  14.4987  27.2785  32.7142
    39.6076

          k =-0.3750 0.3750 0.8750 (   176 PWs)   bands (ev):

    10.9698  11.5109  12.2800  13.2468  14.2186  18.1064  21.5401  27.3703
    34.3959

          k = 0.3750-0.3750 0.6250 (   172 PWs)   bands (ev):

    10.1825  12.1398  12.7501  12.7926  14.4701  15.8906  20.9029  33.7520
    36.0974

          k = 0.3750-0.3750 1.6250 (   174 PWs)   bands (ev):

     9.3305  12.6014  12.6014  12.6765  14.2265  14.2265  23.2891  36.8995
    39.3685

          k = 0.1250-0.1250 0.3750 (   169 PWs)   bands (ev):

     7.3327  11.9983  12.8358  13.0200  13.4875  13.9186  31.3755  36.3333
    36.7643

          k = 0.1250-0.1250 1.3750 (   178 PWs)   bands (ev):

     9.5396  11.3428  12.7065  13.5760  14.3301  14.5163  28.2785  31.5780
    32.3841

          k =-0.1250 0.6250 0.1250 (   178 PWs)   bands (ev):

     9.5396  11.3428  12.7065  13.5760  14.3301  14.5163  28.2785  31.5780
    32.3841

          k =-0.1250 0.6250 1.1250 (   179 PWs)   bands (ev):

    10.8818  11.3221  12.3443  13.6454  14.5133  19.3212  22.5348  26.1704
    28.4085

          k = 0.6250-0.1250 0.8750 (   179 PWs)   bands (ev):

    10.8818  11.3221  12.3443  13.6454  14.5133  19.3212  22.5348  26.1704
    28.4085

          k = 0.6250-0.1250 1.8750 (   178 PWs)   bands (ev):

     9.5396  11.3428  12.7065  13.5760  14.3301  14.5163  28.2785  31.5780
    32.3841

          k = 0.3750 0.1250 0.6250 (   174 PWs)   bands (ev):

    10.3493  11.6766  12.1579  13.2575  14.1340  15.9186  24.3093  29.8490
    32.9693

          k = 0.3750 0.1250 1.6250 (   171 PWs)   bands (ev):

     8.6207  11.9920  12.5953  12.9299  13.5959  14.4987  27.2785  32.7142
    39.6076

          k =-0.1250-0.8750 0.1250 (   176 PWs)   bands (ev):

    10.2090  10.8957  13.6527  14.1097  14.5846  17.0385  25.1835  26.4722
    30.1022

          k =-0.1250-0.8750 1.1250 (   176 PWs)   bands (ev):

    10.2090  10.8957  13.6527  14.1097  14.5846  17.0385  25.1835  26.4722
    30.1022

          k =-0.3750 0.3750 0.3750 (   174 PWs)   bands (ev):

     9.3305  12.6014  12.6014  12.6765  14.2265  14.2265  23.2891  36.8995
    39.3685

          k =-0.3750 0.3750 1.3750 (   172 PWs)   bands (ev):

    10.1825  12.1398  12.7501  12.7926  14.4701  15.8906  20.9029  33.7520
    36.0974

          k = 0.3750-0.3750 1.1250 (   176 PWs)   bands (ev):

    10.9698  11.5109  12.2800  13.2468  14.2186  18.1064  21.5401  27.3703
    34.3959

          k = 0.3750-0.3750 2.1250 (   171 PWs)   bands (ev):

     8.6207  11.9920  12.5953  12.9299  13.5959  14.4987  27.2785  32.7142
    39.6076

          k = 0.3750-0.1250-0.3750 (   171 PWs)   bands (ev):

     8.6207  11.9920  12.5953  12.9299  13.5959  14.4987  27.2785  32.7142
    39.6076

          k = 0.3750-0.1250 0.6250 (   174 PWs)   bands (ev):

    10.3493  11.6766  12.1579  13.2575  14.1340  15.9186  24.3093  29.8490
    32.9693

          k =-0.3750 0.6250 0.3750 (   172 PWs)   bands (ev):

    10.1825  12.1398  12.7501  12.7926  14.4701  15.8906  20.9029  33.7520
    36.0974

          k =-0.3750 0.6250 1.3750 (   172 PWs)   bands (ev):

    10.1825  12.1398  12.7501  12.7926  14.4701  15.8906  20.9029  33.7520
    36.0974

          k =-0.1250 0.3750 0.1250 (   169 PWs)   bands (ev):

     7.3327  11.9983  12.8358  13.0200  13.4875  13.9186  31.3755  36.3333
    36.7643

          k =-0.1250 0.3750 1.1250 (   179 PWs)   bands (ev):

    10.8818  11.3221  12.3443  13.6454  14.5133  19.3212  22.5348  26.1704
    28.4085

          k = 0.1250-0.1250 0.6250 (   178 PWs)   bands (ev):

     9.5396  11.3428  12.7065  13.5760  14.3301  14.5163  28.2785  31.5780
    32.3841

          k = 0.1250-0.1250 1.6250 (   169 PWs)   bands (ev):

     7.3327  11.9983  12.8358  13.0200  13.4875  13.9186  31.3755  36.3333
    36.7643

          k =-0.1250 0.8750 0.6250 (   179 PWs)   bands (ev):

    10.8818  11.3221  12.3443  13.6454  14.5133  19.3212  22.5348  26.1704
    28.4085

          k =-0.1250 0.8750 1.6250 (   179 PWs)   bands (ev):

    10.8818  11.3221  12.3443  13.6454  14.5133  19.3212  22.5348  26.1704
    28.4085

          k = 0.8750 0.6250-0.1250 (   179 PWs)   bands (ev):

    10.8818  11.3221  12.3443  13.6454  14.5133  19.3212  22.5348  26.1704
    28.4085

          k = 0.8750 0.6250 0.8750 (   169 PWs)   bands (ev):

     7.3327  11.9983  12.8358  13.0200  13.4875  13.9186  31.3755  36.3333
    36.7643

          k = 0.1250 0.6250 0.3750 (   174 PWs)   bands (ev):

    10.3493  11.6766  12.1579  13.2575  14.1340  15.9186  24.3093  29.8490
    32.9693

          k = 0.1250 0.6250 1.3750 (   176 PWs)   bands (ev):

    10.9698  11.5109  12.2800  13.2468  14.2186  18.1064  21.5401  27.3703
    34.3959

          k = 0.6250 0.3750 0.1250 (   174 PWs)   bands (ev):

    10.3493  11.6766  12.1579  13.2575  14.1340  15.9186  24.3093  29.8490
    32.9693

          k = 0.6250 0.3750 1.1250 (   174 PWs)   bands (ev):

    10.3493  11.6766  12.1579  13.2575  14.1340  15.9186  24.3093  29.8490
    32.9693

          k = 0.1250-0.1250-0.8750 (   176 PWs)   bands (ev):

    10.2090  10.8957  13.6527  14.1097  14.5846  17.0385  25.1835  26.4722
    30.1022

          k = 0.1250-0.1250 0.1250 (   172 PWs)   bands (ev):

     5.8234  12.4452  12.7307  12.7307  13.5994  13.5994  35.2386  38.9838
    41.0911

          k =-0.3750 1.1250 0.3750 (   176 PWs)   bands (ev):

    10.9698  11.5109  12.2800  13.2468  14.2186  18.1064  21.5401  27.3703
    34.3959

          k =-0.3750 1.1250 1.3750 (   174 PWs)   bands (ev):

    10.3493  11.6766  12.1579  13.2575  14.1340  15.9186  24.3093  29.8490
    32.9693

     the Fermi energy is    14.2874 ev

     Writing output data file nickel.save

     phonons of Ni at X                                                         

     bravais-lattice index     =            2
     lattice parameter (alat)  =       6.6500  a.u.
     unit-cell volume          =      73.5199 (a.u.)^3
     number of atoms/cell      =            1
     number of atomic types    =            1
     kinetic-energy cut-off    =      27.0000  Ry
     charge density cut-off    =     300.0000  Ry
     convergence threshold     =      1.0E-14
     beta                      =       0.7000
     number of iterations used =            4
     Exchange-correlation      =  SLA  PW   PBE  PBE ( 1  4  3  4 0 0)


     celldm(1)=    6.65000  celldm(2)=    0.00000  celldm(3)=    0.00000
     celldm(4)=    0.00000  celldm(5)=    0.00000  celldm(6)=    0.00000

     crystal axes: (cart. coord. in units of alat)
               a(1) = ( -0.5000  0.0000  0.5000 )  
               a(2) = (  0.0000  0.5000  0.5000 )  
               a(3) = ( -0.5000  0.5000  0.0000 )  

     reciprocal axes: (cart. coord. in units 2 pi/alat)
               b(1) = ( -1.0000 -1.0000  1.0000 )  
               b(2) = (  1.0000  1.0000  1.0000 )  
               b(3) = ( -1.0000  1.0000 -1.0000 )  


     Atoms inside the unit cell: 

     Cartesian axes

     site n.  atom      mass           positions (alat units)
        1     Ni  58.6934   tau(    1) = (    0.00000    0.00000    0.00000  )

     Computing dynamical matrix for 
                    q = (   0.0000000   0.0000000   1.0000000 )

     17 Sym.Ops. (with q -> -q+G )


     G cutoff =  336.0507  (   1607 G-vectors)     FFT grid: ( 25, 25, 25)
     G cutoff =  120.9783  (    354 G-vectors)  smooth grid: ( 15, 15, 15)

     number of k points=    80  Marzari-Vanderbilt smearing, width (Ry)=  0.0200
                       cart. coord. in units 2pi/alat
        k(    1) = (  -0.1250000   0.1250000   0.1250000), wk =   0.0312500
        k(    2) = (  -0.1250000   0.1250000   1.1250000), wk =   0.0000000
        k(    3) = (  -0.3750000   0.3750000  -0.1250000), wk =   0.0312500
        k(    4) = (  -0.3750000   0.3750000   0.8750000), wk =   0.0000000
        k(    5) = (   0.3750000  -0.3750000   0.6250000), wk =   0.0312500
        k(    6) = (   0.3750000  -0.3750000   1.6250000), wk =   0.0000000
        k(    7) = (   0.1250000  -0.1250000   0.3750000), wk =   0.0312500
        k(    8) = (   0.1250000  -0.1250000   1.3750000), wk =   0.0000000
        k(    9) = (  -0.1250000   0.6250000   0.1250000), wk =   0.0625000
        k(   10) = (  -0.1250000   0.6250000   1.1250000), wk =   0.0000000
        k(   11) = (   0.6250000  -0.1250000   0.8750000), wk =   0.0625000
        k(   12) = (   0.6250000  -0.1250000   1.8750000), wk =   0.0000000
        k(   13) = (   0.3750000   0.1250000   0.6250000), wk =   0.0625000
        k(   14) = (   0.3750000   0.1250000   1.6250000), wk =   0.0000000
        k(   15) = (  -0.1250000  -0.8750000   0.1250000), wk =   0.0625000
        k(   16) = (  -0.1250000  -0.8750000   1.1250000), wk =   0.0000000
        k(   17) = (  -0.3750000   0.3750000   0.3750000), wk =   0.0312500
        k(   18) = (  -0.3750000   0.3750000   1.3750000), wk =   0.0000000
        k(   19) = (   0.3750000  -0.3750000   1.1250000), wk =   0.0312500
        k(   20) = (   0.3750000  -0.3750000   2.1250000), wk =   0.0000000
        k(   21) = (   0.3750000  -0.1250000  -0.3750000), wk =   0.0625000
        k(   22) = (   0.3750000  -0.1250000   0.6250000), wk =   0.0000000
        k(   23) = (  -0.3750000   0.6250000   0.3750000), wk =   0.0625000
        k(   24) = (  -0.3750000   0.6250000   1.3750000), wk =   0.0000000
        k(   25) = (  -0.1250000   0.3750000   0.1250000), wk =   0.0625000
        k(   26) = (  -0.1250000   0.3750000   1.1250000), wk =   0.0000000
        k(   27) = (   0.1250000  -0.1250000   0.6250000), wk =   0.0312500
        k(   28) = (   0.1250000  -0.1250000   1.6250000), wk =   0.0000000
        k(   29) = (  -0.1250000   0.8750000   0.6250000), wk =   0.0625000
        k(   30) = (  -0.1250000   0.8750000   1.6250000), wk =   0.0000000
        k(   31) = (   0.8750000   0.6250000  -0.1250000), wk =   0.0625000
        k(   32) = (   0.8750000   0.6250000   0.8750000), wk =   0.0000000
        k(   33) = (   0.1250000   0.6250000   0.3750000), wk =   0.0625000
        k(   34) = (   0.1250000   0.6250000   1.3750000), wk =   0.0000000
        k(   35) = (   0.6250000   0.3750000   0.1250000), wk =   0.0625000
        k(   36) = (   0.6250000   0.3750000   1.1250000), wk =   0.0000000
        k(   37) = (   0.1250000  -0.1250000  -0.8750000), wk =   0.0312500
        k(   38) = (   0.1250000  -0.1250000   0.1250000), wk =   0.0000000
        k(   39) = (  -0.3750000   1.1250000   0.3750000), wk =   0.0625000
        k(   40) = (  -0.3750000   1.1250000   1.3750000), wk =   0.0000000
        k(   41) = (  -0.1250000   0.1250000   0.1250000), wk =   0.0312500
        k(   42) = (  -0.1250000   0.1250000   1.1250000), wk =   0.0000000
        k(   43) = (  -0.3750000   0.3750000  -0.1250000), wk =   0.0312500
        k(   44) = (  -0.3750000   0.3750000   0.8750000), wk =   0.0000000
        k(   45) = (   0.3750000  -0.3750000   0.6250000), wk =   0.0312500
        k(   46) = (   0.3750000  -0.3750000   1.6250000), wk =   0.0000000
        k(   47) = (   0.1250000  -0.1250000   0.3750000), wk =   0.0312500
        k(   48) = (   0.1250000  -0.1250000   1.3750000), wk =   0.0000000
        k(   49) = (  -0.1250000   0.6250000   0.1250000), wk =   0.0625000
        k(   50) = (  -0.1250000   0.6250000   1.1250000), wk =   0.0000000
        k(   51) = (   0.6250000  -0.1250000   0.8750000), wk =   0.0625000
        k(   52) = (   0.6250000  -0.1250000   1.8750000), wk =   0.0000000
        k(   53) = (   0.3750000   0.1250000   0.6250000), wk =   0.0625000
        k(   54) = (   0.3750000   0.1250000   1.6250000), wk =   0.0000000
        k(   55) = (  -0.1250000  -0.8750000   0.1250000), wk =   0.0625000
        k(   56) = (  -0.1250000  -0.8750000   1.1250000), wk =   0.0000000
        k(   57) = (  -0.3750000   0.3750000   0.3750000), wk =   0.0312500
        k(   58) = (  -0.3750000   0.3750000   1.3750000), wk =   0.0000000
        k(   59) = (   0.3750000  -0.3750000   1.1250000), wk =   0.0312500
        k(   60) = (   0.3750000  -0.3750000   2.1250000), wk =   0.0000000
        k(   61) = (   0.3750000  -0.1250000  -0.3750000), wk =   0.0625000
        k(   62) = (   0.3750000  -0.1250000   0.6250000), wk =   0.0000000
        k(   63) = (  -0.3750000   0.6250000   0.3750000), wk =   0.0625000
        k(   64) = (  -0.3750000   0.6250000   1.3750000), wk =   0.0000000
        k(   65) = (  -0.1250000   0.3750000   0.1250000), wk =   0.0625000
        k(   66) = (  -0.1250000   0.3750000   1.1250000), wk =   0.0000000
        k(   67) = (   0.1250000  -0.1250000   0.6250000), wk =   0.0312500
        k(   68) = (   0.1250000  -0.1250000   1.6250000), wk =   0.0000000
        k(   69) = (  -0.1250000   0.8750000   0.6250000), wk =   0.0625000
        k(   70) = (  -0.1250000   0.8750000   1.6250000), wk =   0.0000000
        k(   71) = (   0.8750000   0.6250000  -0.1250000), wk =   0.0625000
        k(   72) = (   0.8750000   0.6250000   0.8750000), wk =   0.0000000
        k(   73) = (   0.1250000   0.6250000   0.3750000), wk =   0.0625000
        k(   74) = (   0.1250000   0.6250000   1.3750000), wk =   0.0000000
        k(   75) = (   0.6250000   0.3750000   0.1250000), wk =   0.0625000
        k(   76) = (   0.6250000   0.3750000   1.1250000), wk =   0.0000000
        k(   77) = (   0.1250000  -0.1250000  -0.8750000), wk =   0.0312500
        k(   78) = (   0.1250000  -0.1250000   0.1250000), wk =   0.0000000
        k(   79) = (  -0.3750000   1.1250000   0.3750000), wk =   0.0625000
        k(   80) = (  -0.3750000   1.1250000   1.3750000), wk =   0.0000000

     PseudoPot. # 1 for Ni read from file:
     ./Ni.pbe-nd-rrkjus.UPF
     MD5 check sum: 8081f0a005c9a5470caab1a58e82ecb2
     Pseudo is Ultrasoft + core correction, Zval = 10.0
     Generated by new atomic code, or converted to UPF format
     Using radial grid of 1203 points,  6 beta functions with: 
                l(1) =   0
                l(2) =   0
                l(3) =   1
                l(4) =   1
                l(5) =   2
                l(6) =   2
     Q(r) pseudized with 0 coefficients 


     Mode symmetry, D_4h(4/mmm) point group:


     Atomic displacements:
     There are   2 irreducible representations

     Representation     1      1 modes -A_2u X_4' M_4'  To be done

     Representation     2      2 modes -E_u  X_5' M_5'  To be done



     Alpha used in Ewald sum =   2.8000
     PHONON       :     2.02s CPU         2.21s WALL



     Representation #  1 mode #   1

     Self-consistent Calculation

      iter #   1 total cpu time :     2.4 secs   av.it.:   5.0
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  4.698E-04

      iter #   2 total cpu time :     2.6 secs   av.it.:   6.8
      thresh= 2.168E-03 alpha_mix =  0.700 |ddv_scf|^2 =  2.689E-04

      iter #   3 total cpu time :     2.8 secs   av.it.:   6.2
      thresh= 1.640E-03 alpha_mix =  0.700 |ddv_scf|^2 =  3.460E-08

      iter #   4 total cpu time :     3.0 secs   av.it.:   6.6
      thresh= 1.860E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.151E-10

      iter #   5 total cpu time :     3.2 secs   av.it.:   5.7
      thresh= 1.775E-06 alpha_mix =  0.700 |ddv_scf|^2 =  1.020E-11

      iter #   6 total cpu time :     3.4 secs   av.it.:   6.1
      thresh= 3.194E-07 alpha_mix =  0.700 |ddv_scf|^2 =  5.191E-14

      iter #   7 total cpu time :     3.6 secs   av.it.:   6.6
      thresh= 2.278E-08 alpha_mix =  0.700 |ddv_scf|^2 =  9.422E-16

     End of self-consistent calculation

     Convergence has been achieved 


     Representation #  2 modes #   2  3

     Self-consistent Calculation

      iter #   1 total cpu time :     3.9 secs   av.it.:   4.4
      thresh= 1.000E-02 alpha_mix =  0.700 |ddv_scf|^2 =  1.426E-05

      iter #   2 total cpu time :     4.3 secs   av.it.:   7.7
      thresh= 3.776E-04 alpha_mix =  0.700 |ddv_scf|^2 =  5.396E-07

      iter #   3 total cpu time :     4.7 secs   av.it.:   7.5
      thresh= 7.346E-05 alpha_mix =  0.700 |ddv_scf|^2 =  3.168E-09

      iter #   4 total cpu time :     5.1 secs   av.it.:   7.0
      thresh= 5.628E-06 alpha_mix =  0.700 |ddv_scf|^2 =  5.804E-12

      iter #   5 total cpu time :     5.4 secs   av.it.:   7.0
      thresh= 2.409E-07 alpha_mix =  0.700 |ddv_scf|^2 =  1.502E-13

      iter #   6 total cpu time :     5.8 secs   av.it.:   7.3
      thresh= 3.876E-08 alpha_mix =  0.700 |ddv_scf|^2 =  8.996E-15

     End of self-consistent calculation

     Convergence has been achieved 

     Number of q in the star =    3
     List of q in the star:
          1   0.000000000   0.000000000   1.000000000
          2   0.000000000   1.000000000   0.000000000
          3   1.000000000   0.000000000   0.000000000

     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   1.000000000 ) 

 **************************************************************************
     freq (    1) =       6.621817 [THz] =     220.880032 [cm-1]
     freq (    2) =       6.621817 [THz] =     220.880032 [cm-1]
     freq (    3) =       8.965758 [THz] =     299.065511 [cm-1]
 **************************************************************************

     Mode symmetry, D_4h(4/mmm) point group:

     freq (  1 -  2) =        220.9  [cm-1]   --> E_u  X_5' M_5'     
     freq (  3 -  3) =        299.1  [cm-1]   --> A_2u X_4' M_4'     

     init_run     :      0.19s CPU      0.21s WALL (       1 calls)
     electrons    :      0.83s CPU      0.92s WALL (       1 calls)

     Called by init_run:
     wfcinit      :      0.00s CPU      0.00s WALL (       1 calls)
     potinit      :      0.01s CPU      0.02s WALL (       1 calls)

     Called by electrons:
     c_bands      :      0.83s CPU      0.91s WALL (       1 calls)
     v_of_rho     :      0.02s CPU      0.03s WALL (       2 calls)
     newd         :      0.02s CPU      0.02s WALL (       2 calls)

     Called by c_bands:
     init_us_2    :      0.04s CPU      0.03s WALL (     840 calls)
     cegterg      :      0.76s CPU      0.85s WALL (      80 calls)

     Called by sum_band:

     Called by *egterg:
     h_psi        :      1.36s CPU      1.55s WALL (    7271 calls)
     s_psi        :      0.34s CPU      0.30s WALL (   14013 calls)
     g_psi        :      0.00s CPU      0.00s WALL (    1129 calls)
     cdiaghg      :      0.32s CPU      0.39s WALL (    1209 calls)

     Called by h_psi:
     h_psi:pot    :      1.34s CPU      1.53s WALL (    7271 calls)
     h_psi:calbec :      0.13s CPU      0.16s WALL (    7271 calls)
     vloc_psi     :      1.08s CPU      1.20s WALL (    7271 calls)
     add_vuspsi   :      0.11s CPU      0.15s WALL (    7271 calls)

     General routines
     calbec       :      0.30s CPU      0.36s WALL (   15933 calls)
     fft          :      0.06s CPU      0.06s WALL (     712 calls)
     ffts         :      0.00s CPU      0.00s WALL (     215 calls)
     fftw         :      1.04s CPU      1.17s WALL (   90372 calls)
     interpolate  :      0.01s CPU      0.01s WALL (      86 calls)
     davcio       :      0.03s CPU      0.03s WALL (    3756 calls)

     Parallel routines
     fft_scatter  :      0.38s CPU      0.50s WALL (   91299 calls)

     PHONON       :     5.34s CPU         5.84s WALL

     INITIALIZATION: 
     phq_setup    :      0.06s CPU      0.06s WALL (       1 calls)
     phq_init     :      0.45s CPU      0.47s WALL (       1 calls)

     phq_init     :      0.45s CPU      0.47s WALL (       1 calls)
     set_drhoc    :      0.14s CPU      0.15s WALL (       3 calls)
     init_vloc    :      0.01s CPU      0.01s WALL (       2 calls)
     init_us_1    :      0.32s CPU      0.39s WALL (       2 calls)
     newd         :      0.02s CPU      0.02s WALL (       2 calls)
     dvanqq       :      0.07s CPU      0.07s WALL (       1 calls)
     drho         :      0.13s CPU      0.14s WALL (       1 calls)

     DYNAMICAL MATRIX:
     dynmat0      :      0.14s CPU      0.14s WALL (       1 calls)
     phqscf       :      3.32s CPU      3.62s WALL (       1 calls)
     dynmatrix    :      0.00s CPU      0.00s WALL (       1 calls)

     phqscf       :      3.32s CPU      3.62s WALL (       1 calls)
     solve_linter :      3.29s CPU      3.58s WALL (       2 calls)
     drhodv       :      0.03s CPU      0.03s WALL (       2 calls)

     dynmat0      :      0.14s CPU      0.14s WALL (       1 calls)
     dynmat_us    :      0.03s CPU      0.03s WALL (       1 calls)
     d2ionq       :      0.00s CPU      0.00s WALL (       1 calls)
     dynmatcc     :      0.11s CPU      0.11s WALL (       1 calls)

     dynmat_us    :      0.03s CPU      0.03s WALL (       1 calls)
     addusdynmat  :      0.00s CPU      0.00s WALL (       1 calls)

     phqscf       :      3.32s CPU      3.62s WALL (       1 calls)
     solve_linter :      3.29s CPU      3.58s WALL (       2 calls)

     solve_linter :      3.29s CPU      3.58s WALL (       2 calls)
     dvqpsi_us    :      0.08s CPU      0.11s WALL (     120 calls)
     ortho        :      0.07s CPU      0.08s WALL (     760 calls)
     cgsolve      :      1.68s CPU      1.84s WALL (     760 calls)
     incdrhoscf   :      0.14s CPU      0.16s WALL (     760 calls)
     addusddens   :      0.24s CPU      0.25s WALL (      15 calls)
     vpsifft      :      0.10s CPU      0.11s WALL (     640 calls)
     dv_of_drho   :      0.09s CPU      0.10s WALL (      19 calls)
     mix_pot      :      0.01s CPU      0.02s WALL (      13 calls)
     psymdvscf    :      0.55s CPU      0.56s WALL (      13 calls)
     newdq        :      0.26s CPU      0.27s WALL (      13 calls)
     adddvscf     :      0.01s CPU      0.04s WALL (     640 calls)
     drhodvus     :      0.00s CPU      0.00s WALL (       2 calls)

     dvqpsi_us    :      0.08s CPU      0.11s WALL (     120 calls)
     dvqpsi_us_on :      0.06s CPU      0.07s WALL (     120 calls)

     cgsolve      :      1.68s CPU      1.84s WALL (     760 calls)
     ch_psi       :      1.60s CPU      1.75s WALL (    5982 calls)

     ch_psi       :      1.60s CPU      1.75s WALL (    5982 calls)
     h_psi        :      1.36s CPU      1.55s WALL (    7271 calls)
     last         :      0.33s CPU      0.35s WALL (    5982 calls)

     h_psi        :      1.36s CPU      1.55s WALL (    7271 calls)
     add_vuspsi   :      0.11s CPU      0.15s WALL (    7271 calls)

     incdrhoscf   :      0.14s CPU      0.16s WALL (     760 calls)
     addusdbec    :      0.04s CPU      0.04s WALL (     880 calls)

     drhodvus     :      0.00s CPU      0.00s WALL (       2 calls)

      General routines
     calbec       :      0.30s CPU      0.36s WALL (   15933 calls)
     fft          :      0.06s CPU      0.06s WALL (     712 calls)
     ffts         :      0.00s CPU      0.00s WALL (     215 calls)
     fftw         :      1.04s CPU      1.17s WALL (   90372 calls)
     davcio       :      0.03s CPU      0.03s WALL (    3756 calls)
     write_rec    :      0.02s CPU      0.03s WALL (      15 calls)


     PHONON       :     5.34s CPU         5.84s WALL


   This run was terminated on:   0: 1:40   7Dec2016            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
"""


def test_read_espresso_ph_all():
    for i in range(1, 66, 5):
        to_read = globals()[f"test_{i}"]
        fd = StringIO(to_read)
        read_espresso_ph(fd)


def test_read_espresso_ph_1():
    fd = StringIO(test_1)
    results = read_espresso_ph(fd)

    assert len(results) == 8
    assert (0, 0, 0) in results
    assert np.unique(results[(0, 0, 0)]["freqs"]).shape[0] == 1
    assert np.unique(results[(0, 0, 0)]["freqs"])[0] == 0.173268
    assert len(results[(0, 0, 0)]["eqpoints"]) == 1
    assert results[(0, 0, 0)]["atoms"].symbols == ["Al"]

    assert (0.75, -0.25, 0.75) in results
    assert np.unique(results[(0.75, -0.25, 0.75)]["freqs"]).shape[0] == 3
    assert np.unique(results[(0.75, -0.25, 0.75)]["freqs"])[2] == 8.791383
    assert len(results[(0.75, -0.25, 0.75)]["eqpoints"]) == 24
    assert results[(0.75, -0.25, 0.75)]["atoms"].symbols == ["Al"]
