from htase.calculators.vasp import SmartVasp
from htase.schemas.vasp import summarize
from htase.util.atoms import make_slabs_from_bulk
from ase.io.jsonio import encode, decode
from jobflow import job, Flow, Response

NCORE = 4
KPAR = 4


@job
def run_slab_job(atoms_json, slab=True, static=False):
    atoms = decode(atoms_json)

    if static:
        updates = {
            "ismear": -5,
            "nsw": 0,
            "lwave": True,
            "lcharg": True,
            "laechg": True,
            "lvhar": True,
        }
    else:
        updates = {}

    if slab:
        atoms = SmartVasp(
            atoms, preset="SlabRelaxSet", ncore=NCORE, kpar=KPAR, **updates
        )
    else:
        atoms = SmartVasp(
            atoms,
            preset="SlabRelaxSet",
            isif=3,
            ncore=NCORE,
            kpar=KPAR * 4,
            auto_dipole=False,
            auto_kpts={"length_density": [50, 50, 50]},
            **updates,
        )

    atoms.get_potential_energy()
    results = summarize.get_results(atoms)

    return results


@job
def run_dos(atoms_json):
    atoms = decode(atoms_json)
    atoms = SmartVasp(
        atoms,
        preset="SlabRelaxSet",
        nsw=0,
        ismear=-5,
        isym=2,
        kpar=4,
        nedos=5000,
    )
    atoms.get_potential_energy()
    results = summarize.get_results(atoms=atoms)

    return results


@job
def bulk_to_slab_job(atoms_json, max_slabs=None, **slabgen_kwargs):
    atoms = decode(atoms_json)

    slabs = make_slabs_from_bulk(atoms, max_slabs=max_slabs, **slabgen_kwargs)
    jobs = []
    outputs = []
    for slab in slabs:
        relax_job = run_slab_job(encode(slab))
        jobs.append(relax_job)
        outputs.append(relax_job.output)

        static_job = run_slab_job(relax_job.output["atoms"], static=True)
        jobs.append(static_job)
        outputs.append(static_job.output)

        dos_job = run_dos(static_job.output["atoms"])
        jobs.append(dos_job)
        outputs.append(dos_job.output)

    return Response(replace=Flow(jobs, output=outputs))
