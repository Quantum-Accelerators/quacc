"""Static API contracts checked by mypy; no calculations are executed."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path
    from typing import assert_type

    from ase.atoms import Atoms
    from ase.calculators.emt import EMT
    from ase.md.md import MolecularDynamics
    from ase.md.verlet import VelocityVerlet
    from ase.optimize.optimize import Optimizer
    from emmet.core.eos import EOSDoc
    from pymatgen.core.structure import Structure

    from quacc.atoms.slabs import flip_atoms
    from quacc.recipes.common.eos import eos_flow, eos_subflow
    from quacc.recipes.common.slabs import bulk_to_slabs_subflow, slab_to_ads_subflow
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow
    from quacc.recipes.vasp.qmof import QMOFRelaxSchema
    from quacc.recipes.vasp.slabs import relax_job, static_job
    from quacc.runners.ase import Runner
    from quacc.schemas.ase import Summarize
    from quacc.types import (
        DynSchema,
        MDParams,
        NebResultSchema,
        OptParams,
        OptSchema,
        RunSchema,
        TrajectoryLog,
        VaspSchema,
        VibThermoSchema,
    )
    from quacc.utils.files import find_recent_logfile

    def check_api(
        atoms: Atoms,
        adsorbate: Atoms,
        dynamics: DynSchema,
        qmof: QMOFRelaxSchema,
        thermo: VibThermoSchema,
        vasp_result: VaspSchema,
    ) -> None:
        assert_type(flip_atoms(atoms), Atoms)
        assert_type(flip_atoms(atoms, return_struct=True), Structure)
        assert_type(find_recent_logfile(".", ".log"), Path | None)
        assert_type(Runner(atoms, EMT()).run_calc(), Atoms)
        assert_type(Runner([atoms, atoms.copy()], EMT()).atoms, list[Atoms])
        assert_type(Runner([atoms, atoms.copy()], EMT()).run_neb(), Optimizer)
        neb: NebResultSchema = Summarize().neb(Optimizer(atoms), n_images=2)
        assert_type(neb["ts_atoms"], Atoms)

        md_params: MDParams = {"dynamics": VelocityVerlet}
        opt_params: OptParams = {"optimizer": VelocityVerlet, "relax_cell": False}
        assert_type(md_params["dynamics"], type[MolecularDynamics])
        assert_type(opt_params["relax_cell"], bool)

        assert_type(dynamics["trajectory_log"], list[TrajectoryLog])
        assert_type(qmof["double_relax"], list[VaspSchema])
        assert_type(thermo["results"]["entropy"], float)
        assert_type(thermo["results"]["vib_freqs"], list[float])
        assert_type(bulk_to_slabs_flow(atoms), dict[str, list[RunSchema | OptSchema]])
        assert_type(
            bulk_to_slabs_subflow(atoms, relax_job, static_job),
            dict[str, list[VaspSchema]],
        )
        assert_type(
            slab_to_ads_subflow(atoms, adsorbate, relax_job, static_job),
            dict[str, list[VaspSchema]],
        )
        eos = eos_subflow(vasp_result, relax_job, static_job)
        assert_type(eos["deformed_results"], dict[str, list[VaspSchema]])
        assert_type(eos["eos_doc"], EOSDoc)
        eos = eos_flow(atoms, relax_job, static_job)
        assert_type(eos["deformed_results"], dict[str, list[VaspSchema]])
