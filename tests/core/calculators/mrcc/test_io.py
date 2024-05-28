from __future__ import annotations

import io
from io import StringIO
from pathlib import Path

import pytest
from ase.atoms import Atoms
from ase.calculators.calculator import compare_atoms

from quacc.calculators.mrcc.io import (
    read_energy,
    read_geom_mrccinp,
    read_mrcc_output,
    read_mrcc_outputs,
    write_mrcc,
)


def test_read_geom_mrccinp():
    reference_inputfile = """calc=LNO-CCSD(T)
basis=cc-pVDZ
charge=0
mult=1
geom=xyz
3

H   1.0 0.0 0.0
H   2.0 0.0 0.0
O   3.0 0.0 0.0
"""
    atoms = read_geom_mrccinp(io.StringIO(reference_inputfile))
    atoms_ref = Atoms(
        "H2O", positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    )
    assert not compare_atoms(atoms, atoms_ref, tol=1e-7)


def test_write_mrcc(tmpdir):
    atoms = Atoms("H2O", positions=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    atoms.set_tags([71, 71, 0])

    params = {
        "charge": 0,
        "mult": 1,
        "mrccinput": {
            "calc": "LNO-CCSD(T)",
            "basis": "cc-pVDZ",
            "symm": "off",
            "localcc": "on",
            "lcorthr": "normal",
            "ccprog": "ccsd",
            "ccsdalg": "dfdirect",
            "dfbasis_cor": "cc-pVDZ-RI",
        },
        "mrccblocks": """scfmaxit=1000
usedisk=0""",
    }
    write_mrcc(Path(f"{tmpdir}/MINP"), atoms, params)

    with open(f"{tmpdir}/MINP") as fd:
        generated_inputfile = fd.readlines()

    reference_inputfile = [
        "calc=LNO-CCSD(T)\n",
        "basis=cc-pVDZ\n",
        "symm=off\n",
        "localcc=on\n",
        "lcorthr=normal\n",
        "ccprog=ccsd\n",
        "ccsdalg=dfdirect\n",
        "dfbasis_cor=cc-pVDZ-RI\n",
        "scfmaxit=1000\n",
        "usedisk=0 \n",
        "charge=0\n",
        "mult=1\n",
        "geom=xyz\n",
        "3\n",
        "\n",
        "H   1.0 0.0 0.0\n",
        "H   2.0 0.0 0.0\n",
        "O   3.0 0.0 0.0\n",
        "\n",
        "ghost=serialno\n",
        "1,2",
    ]

    assert generated_inputfile == reference_inputfile
    # Test when geom line is present in mrccblocks
    params = {
        "charge": 0,
        "mult": 1,
        "mrccinput": {
            "calc": "LNO-CCSD(T)",
            "basis": "cc-pVDZ",
            "symm": "off",
            "localcc": "on",
            "lcorthr": "normal",
            "ccprog": "ccsd",
            "ccsdalg": "dfdirect",
            "dfbasis_cor": "cc-pVDZ-RI",
        },
        "mrccblocks": """geom=xyz
3

H   2.0 0.0 0.0
H   2.0 0.0 0.0
O   3.0 0.0 0.0
""",
    }
    write_mrcc(Path(f"{tmpdir}/MINP"), atoms, params)

    with open(f"{tmpdir}/MINP") as fd:
        generated_inputfile = fd.readlines()

    reference_inputfile = [
        "calc=LNO-CCSD(T)\n",
        "basis=cc-pVDZ\n",
        "symm=off\n",
        "localcc=on\n",
        "lcorthr=normal\n",
        "ccprog=ccsd\n",
        "ccsdalg=dfdirect\n",
        "dfbasis_cor=cc-pVDZ-RI\n",
        "geom=xyz\n",
        "3\n",
        "\n",
        "H   2.0 0.0 0.0\n",
        "H   2.0 0.0 0.0\n",
        "O   3.0 0.0 0.0\n",
        " \n",
    ]

    assert generated_inputfile == reference_inputfile


def test_read_mrcc_output():
    reference_dft_outputfile = """...............................................................................

 SUCCESS...
 THE SCF ITERATION HAS CONVERGED!

                   A1  B1  B2  A2
 FINAL ALPHA OCC:   3   1   1   0
 FINAL BETA  OCC:   3   1   1   0

 ***FINAL KOHN-SHAM ENERGY:           -75.8491211955561653 [AU]

 RETURNING FROM SCF ALGORITHM
 ======================================================================

 ************************ 2024-05-28 16:10:52 *************************
                      Normal termination of mrcc.
 **********************************************************************
 """

    reference_cwft_outputfile = """...............................................................................
 ======================================================================

 SUCCESS...
 THE SCF ITERATION HAS CONVERGED!

                   A
 FINAL ALPHA OCC:   5
 FINAL BETA  OCC:   5

 ***FINAL HARTREE-FOCK ENERGY:        -75.7413954285433988 [AU]
...............................................................................
 LMP2 correlation energy [au]:                             -0.205131614733
 Total LMP2 energy [au]:                                  -75.946527043277

 Total MP2 correction for dropped NAFs [au]:                0.000018777971
 Total MP2 correction [au]:                                 0.000018777971

 CPU time for CCSD calculations [min]:                      0.124
 Total CPU time for CCSD [min]:                             0.414
 Wall time for CCSD calculations [min]:                     0.010
 Total wall time for CCSD [min]:                            0.050

 CCSD correlation energy [au]:                             -0.229368662969
 Total CCSD energy [au]:                                  -75.970764091513
 CCSD correlation energy + 0.5 MP2 corrections [au]:       -0.229359273984
 Total LNO-CCSD energy with MP2 corrections [au]:         -75.970754702527

 CPU time for (T) corrections [min]:                        0.013
 Total CPU time for CCSD(T) [min]:                          0.427
 Wall time for (T) corrections [min]:                       0.001
 Total wall time for CCSD(T) [min]:                         0.051

 Total (T) correlation energy contribution [au]:           -0.006949371048
 Total (T) correlation energy+0.5 MP2 correction [au]:     -0.006939982063
 CCSD(T) correlation energy [au]:                          -0.236318034018
 Total CCSD(T) energy [au]:                               -75.977713462561
 CCSD(T) correlation energy + MP2 corrections [au]:        -0.236299256047
 Total LNO-CCSD(T) energy with MP2 corrections [au]:      -75.977694684591

 ======================================================================
 ======================================================================

 ************************ 2024-05-27 13:39:19 *************************
                      Normal termination of mrcc.
 **********************************************************************
 """

    fd = StringIO(reference_dft_outputfile)
    generated_dft_output = read_energy(fd.readlines())

    fd = StringIO(reference_cwft_outputfile)
    generated_cwft_output = read_energy(fd.readlines())

    reference_dft_output = (-2063.959716461294, None, None, None)
    reference_cwft_output = (
        -2061.028349030339,
        -5.5819155543014425,
        -6.241183742647235,
        -6.430030273565713,
    )

    assert generated_dft_output[0] == pytest.approx(reference_dft_output[0])

    assert generated_cwft_output[0] == pytest.approx(reference_cwft_output[0])
    assert generated_cwft_output[1] == pytest.approx(reference_cwft_output[1])
    assert generated_cwft_output[2] == pytest.approx(reference_cwft_output[2])
    assert generated_cwft_output[3] == pytest.approx(reference_cwft_output[3])


def test_read_mrcc_output(tmpdir):
    reference_dft_outputfile = """...............................................................................

 SUCCESS...
 THE SCF ITERATION HAS CONVERGED!

                   A1  B1  B2  A2
 FINAL ALPHA OCC:   3   1   1   0
 FINAL BETA  OCC:   3   1   1   0

 ***FINAL KOHN-SHAM ENERGY:           -75.8491211955561653 [AU]

 RETURNING FROM SCF ALGORITHM
 ======================================================================

 ************************ 2024-05-28 16:10:52 *************************
                      Normal termination of mrcc.
 **********************************************************************
 """

    reference_cwft_outputfile = """...............................................................................
 ======================================================================

 SUCCESS...
 THE SCF ITERATION HAS CONVERGED!

                   A
 FINAL ALPHA OCC:   5
 FINAL BETA  OCC:   5

 ***FINAL HARTREE-FOCK ENERGY:        -75.7413954285433988 [AU]
...............................................................................
 LMP2 correlation energy [au]:                             -0.205131614733
 Total LMP2 energy [au]:                                  -75.946527043277

 Total MP2 correction for dropped NAFs [au]:                0.000018777971
 Total MP2 correction [au]:                                 0.000018777971

 CPU time for CCSD calculations [min]:                      0.124
 Total CPU time for CCSD [min]:                             0.414
 Wall time for CCSD calculations [min]:                     0.010
 Total wall time for CCSD [min]:                            0.050

 CCSD correlation energy [au]:                             -0.229368662969
 Total CCSD energy [au]:                                  -75.970764091513
 CCSD correlation energy + 0.5 MP2 corrections [au]:       -0.229359273984
 Total LNO-CCSD energy with MP2 corrections [au]:         -75.970754702527

 CPU time for (T) corrections [min]:                        0.013
 Total CPU time for CCSD(T) [min]:                          0.427
 Wall time for (T) corrections [min]:                       0.001
 Total wall time for CCSD(T) [min]:                         0.051

 Total (T) correlation energy contribution [au]:           -0.006949371048
 Total (T) correlation energy+0.5 MP2 correction [au]:     -0.006939982063
 CCSD(T) correlation energy [au]:                          -0.236318034018
 Total CCSD(T) energy [au]:                               -75.977713462561
 CCSD(T) correlation energy + MP2 corrections [au]:        -0.236299256047
 Total LNO-CCSD(T) energy with MP2 corrections [au]:      -75.977694684591

 ======================================================================
 ======================================================================

 ************************ 2024-05-27 13:39:19 *************************
                      Normal termination of mrcc.
 **********************************************************************
 """
    generated_cwft_mrcc_output = read_mrcc_output(StringIO(reference_cwft_outputfile))
    generated_dft_mrcc_output = read_mrcc_output(StringIO(reference_dft_outputfile))

    with open(f"{tmpdir}/mrcc_cwft.out", "w") as fd:
        fd.write(reference_cwft_outputfile)

    with open(f"{tmpdir}/mrcc_dft.out", "w") as fd:
        fd.write(reference_dft_outputfile)

    generated_cwft_mrcc_outputs = read_mrcc_outputs(".", f"{tmpdir}/mrcc_cwft.out")
    generated_dft_mrcc_outputs = read_mrcc_outputs(".", f"{tmpdir}/mrcc_dft.out")

    reference_dft_output = {
        "scf_energy": -2063.959716461294,
        "energy": -2063.959716461294,
    }

    reference_cwft_output = {
        "scf_energy": -2061.028349030339,
        "energy": -2067.4583788197597,
        "mp2_corr_energy": -5.5819155543014425,
        "ccsd_corr_energy": -6.2411833455514785,
        "ccsd(t)_corr_energy": -6.4300297894207326,
    }

    for key in reference_dft_output:
        assert generated_dft_mrcc_output[key] == pytest.approx(
            reference_dft_output[key]
        )
        assert generated_dft_mrcc_outputs[key] == pytest.approx(
            reference_dft_output[key]
        )

    for key in reference_cwft_output:
        assert generated_cwft_mrcc_output[key] == pytest.approx(
            reference_cwft_output[key]
        )
        assert generated_cwft_mrcc_outputs[key] == pytest.approx(
            reference_cwft_output[key]
        )
