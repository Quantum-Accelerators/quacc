from pathlib import Path

import numpy as np
import pytest
from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.newtonnet.core import freq_job, relax_job, static_job
from quacc.recipes.newtonnet.ts import irc_job, quasi_irc_job, ts_job

DEFAULT_SETTINGS = SETTINGS.copy()

pytest.importorskip("sella")
pytest.importorskip("newtonnet")


def setup_module():
    current_file_path = Path(__file__).parent.resolve()

    SETTINGS.NEWTONNET_CONFIG_PATH = current_file_path / "config0.yml"
    SETTINGS.NEWTONNET_MODEL_PATH = current_file_path / "best_model_state.tar"
    SETTINGS.CHECK_CONVERGENCE = False


def teardown_module():
    SETTINGS.NEWTONNET_CONFIG_PATH = DEFAULT_SETTINGS.NEWTONNET_CONFIG_PATH
    SETTINGS.NEWTONNET_MODEL_PATH = DEFAULT_SETTINGS.NEWTONNET_MODEL_PATH
    SETTINGS.CHECK_CONVERGENCE = DEFAULT_SETTINGS.CHECK_CONVERGENCE


@pytest.fixture(params=[
    {},
    {'disagreement': 'values'},
    {'method': None},
    {'method': None, 'disagreement': 'values'},
    {'method': 'autograd'},
    {'method': 'autograd', 'disagreement': 'values'},
    {'method': 'fwd_diff', 'grad_precision': 1e-5},
    {'method': 'fwd_diff', 'grad_precision': 1e-5, 'disagreement': 'values'},
    {'method': 'cnt_diff', 'grad_precision': 1e-5},
    {'method': 'cnt_diff', 'grad_precision': 1e-5, 'disagreement': 'values'},
])
def calc_swaps(request):
    return request.param


@pytest.fixture(params=[
    {'optimizer': FIRE},
    {'optimizer': BFGS},
    {'optimizer': BFGSLineSearch},
    {'optimizer': LBFGS},
    {'optimizer': LBFGSLineSearch},
    {'optimizer': GPMin},
    {'optimizer': Sella},
    {'optimizer': IRC},
])
def opt_swaps(request):
    return request.param


def test_static_job(tmpdir, calc_swaps):
    tmpdir.chdir()
    atoms = molecule("H2O")
    output = static_job(atoms, calc_swaps=calc_swaps)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["results"]["energy"] == pytest.approx(-9.532746332918718)


def test_relax_job(tmpdir, calc_swaps, opt_swaps):
    tmpdir.chdir()

    atoms = molecule("H2O")
    output = relax_job(atoms, calc_swaps=calc_swaps, opt_swaps=opt_swaps)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    if opt_swaps['optimizer'] == FIRE:
        assert output["results"]["energy"] == pytest.approx(-9.533049846183454)
    elif opt_swaps['optimizer'] == BFGS:
        assert output["results"]["energy"] == pytest.approx(-9.53305023894913)
    elif opt_swaps['optimizer'] == BFGSLineSearch:
        assert output["results"]["energy"] == pytest.approx(-9.53305023894913)
    elif opt_swaps['optimizer'] == LBFGS:
        assert output["results"]["energy"] == pytest.approx(-9.53305023894913)
    elif opt_swaps['optimizer'] == LBFGSLineSearch:
        assert output["results"]["energy"] == pytest.approx(-9.53305023894913)
    elif opt_swaps['optimizer'] == GPMin:
        assert output["results"]["energy"] == pytest.approx(-9.53305023894913)
    elif opt_swaps['optimizer'] == Sella:
        assert output["results"]["energy"] == pytest.approx(-9.53305023894913)
    elif opt_swaps['optimizer'] == IRC:
        assert output["results"]["energy"] == pytest.approx(-9.53305023894913)

    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


@pytest.mark.parametrize("atoms, "
                         "temperature, "
                         "pressure, "
                         "expected_energy, "
                         "expected_vib_freqs, "
                         "expected_imag_freqs", [
    (molecule("H2O"), 273.15, 1, -9.532746332918718,
     [
         1754.7438570677637,
         3905.8863167962554,
         3928.543281726882
     ],
     []),
    (molecule("CH3"), 1000, 20, -14.338570161038465,
     [949.4267389579275,
      1551.1215815612938,
      1669.107390832925,
      2661.1458956565873,
      2921.245283581416,
      3099.0363742761155],
     [-761.5004719152678]),
])
def test_freq_job(tmpdir,
                  atoms,
                  temperature,
                  pressure,
                  expected_energy,
                  expected_vib_freqs,
                  expected_imag_freqs):
    tmpdir.chdir()
    output = freq_job(atoms,
                      temperature=temperature,
                      pressure=pressure)

    assert output["atoms"] == atoms
    assert len(output["vib"]["results"]["vib_freqs"]) == len(expected_vib_freqs)
    assert output["thermo"]["atoms"] == atoms
    assert output["thermo"]["results"]["energy"] == pytest.approx(expected_energy)


@pytest.fixture
def atoms():
    symbols = ['C', 'O', 'C', 'O', 'H', 'H', 'H', 'H']
    positions = [
        [0.87589219, -0.13926227, 0.10160247],
        [-0.07821864, -1.06257641, 0.39840308],
        [0.38522759, -2.67766458, -0.24094265],
        [1.25339054, -2.63650946, -0.93076096],
        [0.46675900, 0.87142023, 0.48210765],
        [1.06360062, -0.03733460, -0.95079003],
        [1.80266010, -0.30212484, 0.61902953],
        [0.92072327, -2.76503751, 1.36753343]
    ]
    atoms = Atoms(symbols=symbols, positions=positions)
    return atoms


@pytest.mark.parametrize("use_custom_hessian, opt_swaps", [
    (True, None),
    (True, {"max_steps": 4}),
    (False, None),
    (False, {"max_steps": 4})
])
def test_ts_job(tmpdir, atoms, use_custom_hessian, opt_swaps):
    tmpdir.chdir()

    output = ts_job(atoms, use_custom_hessian=use_custom_hessian, opt_swaps=opt_swaps)

    assert isinstance(output, dict)
    assert "freq_job" in output
    assert "thermo" in output["freq_job"]

    if use_custom_hessian:
        if opt_swaps is None:
            assert output["results"]["energy"] == pytest.approx(-30.91756742496904)
            assert output["freq_job"]["vib"]["results"]["vib_energies"][0] == pytest.approx(0.012395738446304833)
        elif "max_steps" in opt_swaps and opt_swaps["max_steps"] == 4:
            assert output["results"]["energy"] == pytest.approx(-29.83567172077247)
            assert output["freq_job"]["vib"]["results"]["vib_energies"][0] == pytest.approx(0.024218215485424726)
    else:
        if opt_swaps is None:
            assert output["results"]["energy"] == pytest.approx(-30.917585898315)
            assert output["freq_job"]["vib"]["results"]["imag_vib_freqs"][0] == pytest.approx(-1676.5015763982956)
        elif "max_steps" in opt_swaps and opt_swaps["max_steps"] == 4:
            assert output["results"]["energy"] == pytest.approx(-29.831788105843444)
            assert output["freq_job"]["vib"]["results"]["imag_vib_freqs"][0] == pytest.approx(-1460.3321973013199)


@pytest.mark.parametrize("opt_swaps, freq_job_kwargs", [
    (None, None),
    ({"max_steps": 500}, None),
    ({"fmax": 0.001}, None),
    (None, {"temperature": 500.0, "pressure": 10.0}),
    ({"run_kwargs": {"direction": "reverse"}}, None),
])
def test_irc_job(tmpdir, atoms, opt_swaps, freq_job_kwargs):
    tmpdir.chdir()
    output = irc_job(atoms, opt_swaps=opt_swaps, freq_job_kwargs=freq_job_kwargs)

    assert isinstance(output, dict)
    print('\n\n\n\nopt_swaps\n\n\n\n', opt_swaps)
    if opt_swaps is None and freq_job_kwargs is None:
        assert output["results"]["energy"] == pytest.approx(-34.449834)
        assert output["freq_job"]["thermo"]["results"]["energy"] == pytest.approx(
            -34.449834
        )
    elif opt_swaps is not None:
        if "max_steps" in opt_swaps and opt_swaps["max_steps"] == 500 and freq_job_kwargs is None:
            assert output["results"]["energy"] == pytest.approx(-34.449834)
            assert output["freq_job"]["thermo"]["results"]["energy"] == pytest.approx(
                -34.449834
            )
    elif opt_swaps is None and\
            ("temperature" in freq_job_kwargs and
             freq_job_kwargs["temperature"] == 500 and
             "pressure" in freq_job_kwargs and freq_job_kwargs["pressure"] == 10.0):
        assert output["results"]["energy"] == pytest.approx(-34.449834)
        assert output["freq_job"]["thermo"]["results"]["energy"] == pytest.approx(
            -34.449834
        )
    elif opt_swaps is not None:
        if "run_kwargs" in opt_swaps and\
                "direction" in opt_swaps["run_kwargs"] and\
                opt_swaps["run_kwargs"]["direction"] == "reverse" and\
                freq_job_kwargs is None:
            assert output["results"]["energy"] == pytest.approx(-34.450195)
            assert output["freq_job"]["thermo"]["results"]["energy"] == pytest.approx(
                -34.450195
            )


@pytest.mark.parametrize("atoms, direction, freq_job_kwargs, irc_job_kwargs, expected_energy", [
    (molecule("H2O"), None, None, None, -9.533050993135458),
    (molecule("H2O"), "reverse", None, None, -9.533050993135458),
    (molecule("H2O"), None, {"temperature": 500.0, "pressure": 10.0}, None, -9.533050993135458),
    (molecule("H2O"), None, None, {"calc_swaps": {"run_kwargs": {"direction": "reverse"}}}, -9.533050993135458),
])
def test_quasi_irc_job(tmpdir,
                       atoms,
                       direction,
                       freq_job_kwargs,
                       irc_job_kwargs,
                       expected_energy):
    tmpdir.chdir()
    output = quasi_irc_job(atoms,
                           direction=direction,
                           freq_job_kwargs=freq_job_kwargs,
                           irc_job_kwargs=irc_job_kwargs)

    assert isinstance(output, dict)
    assert "irc_job" in output
    assert output["irc_job"]["results"]["energy"] == pytest.approx(expected_energy)
    assert output["results"]["energy"] == pytest.approx(expected_energy)
    assert output["freq_job"]["thermo"]["results"]["energy"] == pytest.approx(expected_energy)
