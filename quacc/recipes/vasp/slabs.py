from typing import Any, Dict
from ase.io.jsonio import encode, decode
from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.slabs import make_max_slabs_from_bulk, make_adsorbate_structures
from jobflow import job, Flow, Response, Maker
from dataclasses import dataclass


@dataclass
class SlabRelaxMaker(Maker):
    """
    Class to relax a slab.

    Parameters
    ----------
    name:
        Name of the job.
    preset:
        Preset to use.
    ncore:
        VASP NCORE parameter.
    kpar:
        VASP KPAR parameter.
    """

    name: str = "SlabRelax"
    preset: None | str = None
    ncore: int = 1
    kpar: int = 1

    @job
    def make(self, atoms_json: str, **kwargs) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms_json:
            Encoded .Atoms object

        Returns
        -------
        Dict:
            Summary of the run.
        """
        atoms = decode(atoms_json)
        flags = {
            "auto_dipole": True,
            "ediff": 1e-5,
            "ediffg": -0.02,
            "isif": 2,
            "ibrion": 2,
            "ismear": 0,
            "isym": 0,
            "kpar": self.kpar,
            "lcharg": False,
            "lwave": False,
            "ncore": self.ncore,
            "nsw": 200,
            "sigma": 0.05,
        }
        for k, v in kwargs.items():
            flags[k] = v

        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms.get_potential_energy()
        summary = summarize_run(atoms)

        return summary


@dataclass
class SlabStaticMaker(Maker):
    """
    Class to carry out a single-point calculation on a slab.

    Parameters
    ----------
    name:
        Name of the job.
    preset:
        Preset to use.
    ncore:
        VASP NCORE parameter.
    kpar:
        VASP KPAR parameter.
    """

    name: str = "SlabStatic"
    preset: None | str = None
    ncore: int = 1
    kpar: int = 1

    @job
    def make(self, atoms_json: str, **kwargs) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms_json:
            Encoded .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        atoms = decode(atoms_json)
        flags = {
            "auto_dipole": True,
            "ediff": 1e-6,
            "ismear": -5,
            "isym": 2,
            "kpar": self.kpar,
            "laechg": True,
            "lcharg": True,
            "lvhar": True,
            "lwave": True,
            "ncore": self.ncore,
            "nedos": 5001,
            "nsw": 0,
            "sigma": 0.05,
        }
        for k, v in kwargs.items():
            flags[k] = v

        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms.get_potential_energy()
        summary = summarize_run(atoms)

        return summary


@dataclass
class BulkToSlabMaker(Maker):
    """
    Class to convert a bulk structure to a slab,
    along with the relaxations and statics for the slabs.

    Parameters
    ----------
    name:
        Name of the job.
    preset:
        Preset to use.
    ncore:
        VASP NCORE parameter.
    kpar:
        VASP KPAR parameter.
    """

    name: str = "BulkToSlab"
    preset: None | str = None
    ncore: int = 1
    kpar: int = 1

    @job
    def make(
        self, atoms_json: str, max_slabs: None | int = None, **slab_kwargs
    ) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms_json
            Encoded .Atoms object
        max_slabs
            Maximum number of slabs to make
        **slab_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs.
        """
        atoms = decode(atoms_json)
        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slab_kwargs)
        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = SlabRelaxMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(encode(slab))
            jobs.append(relax_job)
            outputs.append(relax_job.output)

            static_job = SlabStaticMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))


@dataclass
class SlabToAdsSlabMaker(Maker):
    """
    Class to convert a slab structure to one with adsorbates present,
    along with the relaxations and statics for the slab-adsorbate systems.

    Parameters
    ----------
    name:
        Name of the job.
    preset:
        Preset to use.
    ncore:
        VASP NCORE parameter.
    kpar:
        VASP KPAR parameter.
    """

    name: str = "SlabToAdsSlab"
    preset: None | str = None
    npar: int = 1
    kpar: int = 1

    @job
    def make(self, atoms_json: str, adsorbate_json: str, **slab_ads_kwargs) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms_json
            Encoded .Atoms object for the structure.
        adsorbate_json
            Encoded .Atoms object for the adsorbate.
        **slab_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs with adsorbates.
        """
        atoms = decode(atoms_json)
        adsorbate = decode(adsorbate_json)

        slabs = make_adsorbate_structures(atoms, adsorbate, **slab_ads_kwargs)
        if slabs is None:
            return None

        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = SlabRelaxMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(encode(slab))
            jobs.append(relax_job)
            outputs.append(relax_job.output)

            static_job = SlabStaticMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))
