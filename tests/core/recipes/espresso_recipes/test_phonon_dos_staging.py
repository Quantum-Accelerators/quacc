from __future__ import annotations

import gzip
from copy import deepcopy
from typing import TYPE_CHECKING, Any

import pytest

from quacc.recipes.espresso.phonons import phonon_dos_flow
from quacc.utils.files import copy_decompress_files

if TYPE_CHECKING:
    from pathlib import Path

    from quacc.calculators.espresso.espresso import EspressoTemplate
    from quacc.types import CopyFiles


@pytest.fixture(
    params=[("", False), ("", True), (".xml", False), (".xml", True)],
    ids=["text", "text-gzip", "xml", "xml-gzip"],
)
def staged_jobs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, request: pytest.FixtureRequest
) -> dict[str, Any]:
    """Replace QE execution while exercising the real file-transfer mechanism."""
    suffix, compressed = request.param
    calls: dict[str, dict[str, Any]] = {}

    def write_file(directory: Path, name: str, content: str) -> None:
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if compressed:
            with gzip.open(path.with_name(path.name + ".gz"), "wt") as handle:
                handle.write(content)
        else:
            path.write_text(content)

    def run_and_summarize(
        *,
        template: EspressoTemplate,
        copy_files: CopyFiles | None = None,
        calc_swaps: dict[str, Any] | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        binary = template.binary
        directory = tmp_path / binary
        directory.mkdir()
        for spec in copy_files or []:
            copy_decompress_files(spec["source"], spec["filenames"], directory)
        calls[binary] = {
            "copy_files": deepcopy(copy_files),
            "prev_outdir": template.outdir,
            "calc_swaps": deepcopy(calc_swaps),
        }
        if binary == "ph":
            for name in ("matdyn0", "matdyn1"):
                write_file(directory, name + suffix, name)
            write_file(directory, "pwscf.save/wfc1.dat", "wavefunctions")
        elif binary == "q2r":
            assert (directory / ("matdyn0" + suffix)).read_text() == "matdyn0"
            assert (directory / ("matdyn1" + suffix)).read_text() == "matdyn1"
            write_file(directory, "q2r.fc" + suffix, "force constants")
        elif binary == "matdyn":
            assert (directory / ("q2r.fc" + suffix)).read_text() == "force constants"
            write_file(directory, "matdyn.dos", "0 0 0\n1 3 3\n")
        write_file(directory, f"{binary}.out", f"{binary} output")
        return {"dir_name": directory, "results": {}}

    monkeypatch.setattr(
        "quacc.recipes.espresso.phonons.run_and_summarize", run_and_summarize
    )
    return {
        "root": tmp_path,
        "suffix": suffix,
        "compressed": compressed,
        "calls": calls,
    }


def test_phonon_dos_stages_only_required_files(staged_jobs: dict[str, Any]) -> None:
    root = staged_jobs["root"]
    scf = root / "scf"
    scf.mkdir()
    (scf / "custom.dat").write_text("user-supplied input")
    copy_files = [{"source": scf, "filenames": "custom.dat"}]
    job_params = {"phonon_job": {"input_data": {"inputph": {"nq1": 2}}}}
    original_params = deepcopy(job_params)

    results = phonon_dos_flow(
        copy_files=copy_files,
        prev_outdir=scf,
        job_params=job_params,
        job_decorators={"all": None},
    )

    assert set(results) == {"phonon_job", "q2r_job", "matdyn_job"}
    for name, binary in [
        ("phonon_job", "ph"),
        ("q2r_job", "q2r"),
        ("matdyn_job", "matdyn"),
    ]:
        assert results[name]["dir_name"] == root / binary
    assert staged_jobs["calls"]["ph"]["copy_files"] == copy_files
    assert staged_jobs["calls"]["ph"]["prev_outdir"] == scf
    assert staged_jobs["calls"]["ph"]["calc_swaps"]["input_data"]["inputph"]["nq1"] == 2
    assert job_params == original_params
    assert (root / "ph/custom.dat").read_text() == "user-supplied input"
    assert not (root / "q2r/custom.dat").exists()
    assert not (root / "q2r/ph.out").exists()
    assert not (root / "q2r/pwscf.save").exists()
    assert not list((root / "matdyn").glob("matdyn[0-9]*"))
    assert not (root / "matdyn/q2r.out").exists()
    # Staging must leave the original calculation's restart data intact.
    compression = ".gz" if staged_jobs["compressed"] else ""
    assert (root / f"ph/pwscf.save/wfc1.dat{compression}").exists()
