"""Public optional-backend contracts, checked without executing calculations."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path
    from typing import assert_type

    from ase.atoms import Atoms

    from quacc.recipes.orca.fairchem import omol_static_job
    from quacc.recipes.tblite.core import static_job as tblite_static_job
    from quacc.recipes.torchsim._base import (
        AutobatcherDict,
        TrajectoryReporterDetails,
        TrajectoryReporterDict,
    )
    from quacc.types import BaderSchema, RunSchema

    def check_optional_api(
        atoms: Atoms,
        path: Path,
        reporter: TrajectoryReporterDetails,
        bader: BaderSchema,
    ) -> None:
        assert_type(omol_static_job(atoms), RunSchema)
        assert_type(tblite_static_job(atoms), RunSchema)
        assert_type(reporter["filenames"], list[str] | None)
        assert_type(bader["atomic_volume"], list[float])
        autobatcher: AutobatcherDict = {
            "memory_scales_with": "n_edges",
            "max_atoms_to_try": 100,
            "max_iterations": None,
        }
        assert_type(autobatcher, AutobatcherDict)
        trajectory: TrajectoryReporterDict = {"filenames": [path], "state_frequency": 1}
        assert_type(trajectory, TrajectoryReporterDict)
