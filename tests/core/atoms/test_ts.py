from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture(scope="module", autouse=True)
def set_seed():
    np.random.seed(42)  # noqa: NPY002


from importlib.util import find_spec

from ase.atoms import Atoms
from ase.build import molecule

from quacc.atoms.ts import geodesic_interpolate_wrapper

has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))


@pytest.fixture()
def setup_test_environment(tmp_path):
    reactant = Atoms(
        symbols="CCHHCHH",
        positions=[
            [1.4835950817281542, -1.0145410211301968, -0.13209027203235943],
            [0.8409564131524673, 0.018549610257914483, -0.07338809662321308],
            [-0.6399757891931867, 0.01763740851518944, 0.0581573443268891],
            [-1.0005576455546672, 1.0430257532387608, 0.22197240310602892],
            [1.402180736662139, 0.944112416574632, -0.12179540364365492],
            [-1.1216961389434357, -0.3883639833876232, -0.8769102842015071],
            [-0.9645026578514683, -0.6204201840686793, 0.9240543090678239],
        ],
    )

    product = Atoms(
        symbols="CCHHCHH",
        positions=[
            [1.348003553501624, 0.4819311116778978, 0.2752537177143993],
            [0.2386618286631742, -0.3433222966734429, 0.37705518940917926],
            [-0.9741307940518336, 0.07686022294949588, 0.08710778043683955],
            [-1.8314843503320921, -0.5547344604780035, 0.1639037492534953],
            [0.3801391040059668, -1.3793340533058087, 0.71035902765307],
            [1.9296265384257907, 0.622088341468767, 1.0901733942191298],
            [-1.090815880212625, 1.0965111343610956, -0.23791518420660265],
        ],
    )
    return reactant, product


@pytest.mark.skipif(
    not has_geodesic_interpolate,
    reason="geodesic_interpolate function is not available",
)
def test_geodesic_interpolate_wrapper(setup_test_environment):
    reactant, product = setup_test_environment

    # Execute the geodesic_interpolate_wrapper function
    smoother_path = geodesic_interpolate_wrapper(
        reactant,
        product,
        n_images=20,
        redistribute_tol=5e-4,
        smoother_tol=1e-4,
        max_iterations=10,
        max_micro_iterations=10,
        morse_scaling=1.5,
        geometry_friction=1e-2,
        distance_cutoff=2.5,
    )
    assert smoother_path[1].positions[0][0] == pytest.approx(1.378384900, abs=1e-5)
    assert smoother_path[5].positions[0][2] == pytest.approx(-0.512075394, abs=1e-5)


@pytest.mark.skipif(
    not has_geodesic_interpolate,
    reason="geodesic_interpolate function is not available",
)
def test_geodesic_interpolate_wrapper_large_system(setup_test_environment):
    # Test with large system to trigger sweeping updates
    smoother_path = geodesic_interpolate_wrapper(molecule("C60"), molecule("C60"))
    assert len(smoother_path) == 10
