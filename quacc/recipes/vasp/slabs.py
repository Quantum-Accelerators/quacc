"""Recipes for slabs"""
import warnings
from dataclasses import dataclass
from typing import Any, Dict, List

from ase.atoms import Atoms
from jobflow import Flow, Maker, Response, job

from quacc.calculators.vasp import SmartVasp
from quacc.recipes.vasp.core import RelaxJob, StaticJob
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc
from quacc.util.slabs import (
    get_surface_energy,
    make_adsorbate_structures,
    make_max_slabs_from_bulk,
)


@dataclass
class SlabStaticJob(Maker):
    """
    Class to carry out a single-point calculation on a slab.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "VASP-SlabStatic"
    preset: str = None
    swaps: Dict[str, Any] = None

    @job
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        swaps = self.swaps or {}
        defaults = {
            "auto_dipole": True,
            "ismear": -5,
            "isym": 2,
            "laechg": True,
            "lcharg": True,
            "lvhar": True,
            "lwave": True,
            "nedos": 5001,
            "nsw": 0,
        }
        flags = merge_dicts(defaults, swaps)

        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class SlabRelaxJob(Maker):
    """
    Class to relax a slab.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "VASP-SlabRelax"
    preset: str = None
    swaps: Dict[str, Any] = None

    @job
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict:
            Summary of the run.
        """
        swaps = self.swaps or {}
        defaults = {
            "auto_dipole": True,
            "ediffg": -0.02,
            "isif": 2,
            "ibrion": 2,
            "isym": 0,
            "lcharg": False,
            "lwave": False,
            "nsw": 200,
        }
        flags = merge_dicts(defaults, swaps)

        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class BulkToSlabsJob(Maker):
    """
    Class to convert a bulk structure to slabs,
    along with the relaxations and statics for the slabs.

    Parameters
    ----------
    name
        Name of the job.
    slab_relax_job
        Maker to use for the relaxation of the slab.
    slab_static_job
        Default to use for the static calculation of the slab.
    """

    name: str = "VASP-BulkToSlabs"
    slab_relax_job: Maker = SlabRelaxJob()
    slab_static_job: Maker = SlabStaticJob()

    @job
    def make(self, atoms: Atoms, max_slabs: int = None, **slabgen_kwargs) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        max_slabs
            Maximum number of slabs to make. None implies no upper limit.
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs.
        """
        slabgen_kwargs = slabgen_kwargs or {}

        # Generate all the slab
        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slabgen_kwargs)

        # Generate the jobs for each slab
        jobs = []
        outputs = []
        all_atoms = []
        for slab in slabs:
            relax_job = self.slab_relax_job.make(slab)
            static_job = self.slab_static_job.make(relax_job.output["atoms"])

            jobs += [relax_job, static_job]
            outputs.append(static_job.output)
            all_atoms.append(static_job.output["atoms"])

        if len(slabs) == 0:
            return Response(
                output={"input_bulk": atoms, "generated_slabs": None},
                stop_children=True,
            )

        return Response(
            output={"input_bulk": atoms, "generated_slabs": slabs},
            replace=Flow(
                jobs,
                output={"all_atoms": all_atoms, "all_outputs": outputs},
                name=self.name,
            ),
        )


@dataclass
class SlabToAdsorbatesJob(Maker):
    """
    Class to convert a slab (or slabs) to one with an adsorbate present,
    along with the relaxations and statics for the slab-adsorbate systems.
    Multiple slab-adsorbate systems will be generated, one for each unique
    binding site.

    Parameters
    ----------
    name
        Name of the job.
    slab_ads_relax_job
        Maker to use for the relaxation of the slab-adsorbate system.
    slab_ads_static_job
        Maker to use for the static calculation of the slab-adsorbate system.
    """

    name: str = "VASP-SlabToAdsorbates"
    slab_ads_relax_job: Maker = SlabRelaxJob()
    slab_ads_static_job: Maker = SlabStaticJob()

    @job
    def make(
        self,
        slabs: Atoms | List[Atoms],
        adsorbates: Atoms | List[Atoms],
        **make_ads_kwargs
    ) -> Response:
        """
        Make the run.

        Parameters
        ----------
        slabs
            .Atoms object for the slab structure. Also takes a list of Atoms objects
            for the creation of a series of slabs with adsorbates.
        adsorbates
            .Atoms object for the adsorbate. Also takes a list of Atoms objects for
            the consideration of multiple different adsorbates.
        **make_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs with adsorbates.
        """
        make_ads_kwargs = make_ads_kwargs or {}

        if isinstance(slabs, Atoms):
            slabs_list = [slabs]
        else:
            slabs_list = slabs
        if isinstance(adsorbates, Atoms):
            adsorbates_list = [adsorbates]
        else:
            adsorbates_list = adsorbates

        all_ads_slabs = {}
        jobs = []
        outputs = []
        all_atoms = []
        for slab in slabs_list:
            for adsorbate in adsorbates_list:

                # Make slab-adsorbate systems
                ads_slabs = make_adsorbate_structures(
                    slab, adsorbate, **make_ads_kwargs
                )
                all_ads_slabs[adsorbate.get_chemical_formula()] = ads_slabs

                # Make a relaxation+static job for each slab-adsorbate ysstem
                for ads_slab in ads_slabs:
                    relax_job = self.slab_ads_relax_job.make(ads_slab)
                    static_job = self.slab_ads_static_job.make(
                        relax_job.output["atoms"]
                    )

                    jobs += [relax_job, static_job]
                    outputs.append(static_job.output)
                    all_atoms.append(static_job.output["atoms"])

        if len(jobs) == 0:
            return Response(
                output={"input_slabs": slabs, "generated_slab_ads": None},
                stop_children=True,
            )

        return Response(
            output={"input_slabs": slabs, "generated_slab_ads": all_ads_slabs},
            replace=Flow(
                jobs,
                output={"all_atoms": all_atoms, "all_outputs": outputs},
                name=self.name,
            ),
        )


@dataclass
class BulkToAdsorbatesFlow(Maker):
    """
    Flow consisting of:
    1. Bulk relaxation (optional)
    2. Bulk static (optional)
    3. Slab generation
    4. Selection of the most stable slab (optional)
    5. Addition of adsorbates to the slabs or most stable slab
    6. Slab relaxation(s) with adsorbates at multiple binding sites
    7. Slab static(s) with adsorbates at multiple binding sites

    Parameters
    ----------
    name
        Name of the job.
    bulk_relax_job
        Maker to use for the relaxation of the bulk structure.
    bulk_static_job
        Default to use for the static calculation of the bulk structure.
    bulk_to_slabs_job
        Maker to use for the BulkToSlabsJob.
    slab_to_adsorbate_job
        Maker to use for the SlabToAdsorbatesJob.
    """

    name: str = "VASP-BulkToAdsorbates"
    bulk_relax_job: Maker | None = RelaxJob()
    bulk_static_job: Maker | None = StaticJob()
    bulk_to_slabs_job: Maker = BulkToSlabsJob()
    slab_to_adsorbates_job: Maker = SlabToAdsorbatesJob()

    def make(
        self,
        atoms: Atoms,
        adsorbate: Atoms | List[Atoms],
        n_stable_slabs: int = None,
        max_slabs: int = None,
        slabgen_kwargs: Dict[str, Any] = None,
        make_ads_kwargs: Dict[str, Any] = None,
    ) -> Flow:
        """
        Make the Flow.

        Parameters
        ----------
        atoms
            .Atoms object for the structure.
        adsorbate
            .Atoms object for the adsorbate. Can also take a list of adsorbates.
        n_stable_slabs
            If set to an integer, the n most stable slabs will be considered
            for adsorption. If None, then all will be considered.
        max_slabs
            Maximum number of slabs to make. None implies no upper limit.
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()
        make_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()

        Returns
        -------
        Flow
            The Flow for this process.
        """
        jobs = []
        slabgen_kwargs = slabgen_kwargs or {}
        make_ads_kwargs = make_ads_kwargs or {}

        if n_stable_slabs and not self.bulk_static_job:
            raise ValueError("Cannot use n_stable_slabs without a bulk_static_job")

        if self.bulk_relax_job:
            bulk_relax_job = self.bulk_relax_job.make(atoms)
            atoms = bulk_relax_job.output["atoms"]
            jobs.append(bulk_relax_job)

        if self.bulk_static_job:
            bulk_static_job = self.bulk_static_job.make(atoms)
            atoms = bulk_static_job.output["atoms"]
            jobs.append(bulk_static_job)

        bulk_to_slabs_job = self.bulk_to_slabs_job.make(
            atoms, max_slabs=max_slabs, **slabgen_kwargs
        )
        jobs.append(bulk_to_slabs_job)

        if n_stable_slabs:
            find_stable_slab_job = _get_slab_stability(
                bulk_static_job.output,
                bulk_to_slabs_job.output["all_outputs"],
                n_stable_slabs=n_stable_slabs,
            )
            slab_to_adsorbates_job = self.slab_to_adsorbates_job.make(
                find_stable_slab_job.output["stable_slabs"]["all_atoms"],
                adsorbate,
                **make_ads_kwargs
            )
            jobs += [find_stable_slab_job, slab_to_adsorbates_job]
        else:
            slab_to_adsorbates_job = self.slab_to_adsorbates_job.make(
                bulk_to_slabs_job.output["all_atoms"], adsorbate, **make_ads_kwargs
            )
            jobs.append(slab_to_adsorbates_job)

        return Flow(jobs, output=slab_to_adsorbates_job.output, name=self.name)


@job
def _get_slab_stability(
    bulk_summary: Dict[str, Any],
    slab_summaries: Dict[str, Any],
    n_stable_slabs: int = 1,
) -> Dict[str, Any]:
    """
    A job that determine the most stable surface slabs (based on surface energy) for
    a given bulk summary and list of slab summaries.

    Parameters
    ----------
    bulk_summary
        Output of a VASP job corresponding to the bulk structure.
    slab_summaries
        List of outputs of VASP jobs corresponding to the slab structures.
    n_stable_slabs
        The n most stable slabs are returned.

    Returns
    -------
    Dict
        VASP output summaries for the stable and unstable slabs formatted as
        {
            {"stable_slabs": {"all_atoms": [.Atoms, .Atoms, ...], "all_outputs": [...],
            {"unstable_slabs": {"all_atoms": [.Atoms, .Atoms, ...], "all_outputs": [...]
        }
    """
    bulk = bulk_summary["atoms"]
    bulk_energy = bulk_summary["output"]["energy"]
    surface_energies = []

    if n_stable_slabs > len(slab_summaries):
        warnings.warn(
            "n_stable_slabs is larger than the number of slabs. Setting n_stable_slabs to the number of slabs."
        )
        n_stable_slabs = len(slab_summaries)

    # Iterate through each slab summary and determine the most stable slab
    for slab_summary in slab_summaries:
        slab = slab_summary["atoms"]
        slab_energy = slab_summary["output"]["energy"]

        # Calculate the surface energy
        surface_energy = get_surface_energy(bulk, slab, bulk_energy, slab_energy)

        # Insert the surface energy into the slab summary
        slab_summary["avg_surface_energy"] = surface_energy

        # Store the slab energy in a convenient dict
        surface_energies.append(surface_energy)

    slab_summaries_sorted = [
        slab_summary
        for _, slab_summary in sorted(zip(surface_energies, slab_summaries))
    ]
    stable_slab_summaries = slab_summaries_sorted[0:n_stable_slabs]
    unstable_slab_summaries = slab_summaries_sorted[n_stable_slabs:]

    output = {
        "stable_slabs": {
            "all_atoms": [summary["atoms"] for summary in stable_slab_summaries],
            "all_outputs": stable_slab_summaries,
        },
        "unstable_slabs": {
            "all_atoms": [summary["atoms"] for summary in unstable_slab_summaries]
            if unstable_slab_summaries
            else [],
            "all_outputs": unstable_slab_summaries,
        },
    }

    return output
