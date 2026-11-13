from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING, Any

import pytest

from quacc.recipes.espresso.phonons import _PHONON_ASR_DEFAULTS, phonon_dos_flow

if TYPE_CHECKING:
    from pathlib import Path

    from quacc.calculators.espresso.espresso import EspressoTemplate


@pytest.mark.parametrize("asr", ["simple", "all", "no"])
def test_phonon_dos_asr_defaults_and_overrides(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, asr: str
) -> None:
    """Keep defaults stable across overridden and repeated flow invocations."""
    calls: dict[str, Any] = {}
    original_defaults = deepcopy(_PHONON_ASR_DEFAULTS)

    def run_and_summarize(
        *, template: EspressoTemplate, calc_swaps: dict[str, Any], **kwargs: Any
    ) -> dict[str, Any]:
        calls[template.binary] = deepcopy(calc_swaps)
        return {"dir_name": tmp_path / template.binary}

    monkeypatch.setattr(
        "quacc.recipes.espresso.phonons.run_and_summarize", run_and_summarize
    )
    phonon_dos_flow()
    default_calls = deepcopy(calls)
    assert calls["q2r"]["input_data"]["input"]["zasr"] == "crystal"
    assert calls["matdyn"]["input_data"]["input"] == {
        "asr": "crystal",
        "dos": True,
        "nk1": 32,
        "nk2": 32,
        "nk3": 32,
    }
    assert calls["ph"]["input_data"]["inputph"] == {
        "tr2_ph": 1e-12,
        "alpha_mix(1)": 0.1,
        "verbosity": "high",
        "ldisp": True,
        "nq1": 4,
        "nq2": 4,
        "nq3": 4,
    }
    assert original_defaults == _PHONON_ASR_DEFAULTS

    job_params = {"matdyn_job": {"input_data": {"input": {"asr": asr}}}}
    original_params = deepcopy(job_params)
    phonon_dos_flow(job_params=job_params)
    assert calls["matdyn"]["input_data"]["input"]["asr"] == asr
    assert calls["q2r"] == default_calls["q2r"]
    assert calls["ph"] == default_calls["ph"]
    assert job_params == original_params
    assert original_defaults == _PHONON_ASR_DEFAULTS

    phonon_dos_flow()
    assert calls == default_calls
    assert original_defaults == _PHONON_ASR_DEFAULTS
