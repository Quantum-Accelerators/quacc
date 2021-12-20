from ase.calculators.vasp import Vasp
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath
import warnings
import numpy as np
import os


def SmartVasp(
    atoms, force_gamma=True, incar_copilot=True, verbose=True, **kwargs,
):
    """
    This is a wrapper around the VASP calculator that adjusts INCAR parameters on-the-fly.
    Also supports several automatic k-point generation schemes from Pymatgen.

    Parameters
    ----------
    atoms : ase.Atoms
        The Atoms object to be used for the calculation.
    force_gamma : bool
        If True, the KPOINTS will be set to gamma-centered.
        Defaults to True.
    incar_copilot : bool
        If True, the INCAR parameters will be adjusted if they go against the VASP manual.
        Defaults to True.
    verbose : bool
        If True, warnings will be raised when INCAR parameters are changed.
        Defaults to True.
    **kwargs :
        Additional arguments to be passed to the VASP calculator, e.g. xc='PBE', encut=520.
    """

    if incar_copilot:

        # Move final magmoms to initial magmoms if present and any
        # are > 0.02 in magnitude (unless the user has specified some)
        if np.all(atoms.get_initial_magnetic_moments() == 0):
            try:
                mags = atoms.get_magnetic_moments()
            except:
                mags = None
            if mags and np.any(np.abs(mags > 0.02)):
                atoms.set_initial_magnetic_moments(mags)

    # Initialize calculator
    calc = Vasp(atoms=atoms, **kwargs)

    # Shortcuts for pymatgen k-point generation schemes.
    # Options include: line_density (for band structures),
    # reciprocal_density (by volume), grid_density (by number of atoms),
    # max_mixed_density (max of reciprocal_density volume or atoms), and length_density (good for slabs).
    # These are formatted as {"line_density": float}, {"reciprocal_density": float},
    # {"grid_density": float}, {"vol_kkpa_density": [float, float]}, and
    # {"length_density": [float, float, float]}.
    if isinstance(calc.kpts, dict):
        struct = AseAtomsAdaptor().get_structure(atoms)

        if "line_density" in calc.kpts:
            kpath = HighSymmKpath(struct, path_type="latimer_munro")
            kpts, _ = kpath.get_kpoints(
                line_density=calc.kpts["line_density"], coords_are_cartesian=True
            )
            calc.set(kpts=kpts, reciprocal=True)

        else:
            if "max_mixed_density" in calc.kpts:
                pmg_kpts1 = Kpoints.automatic_density_by_vol(
                    struct, calc.kpts["max_mixed_density"][0], force_gamma
                )
                pmg_kpts2 = Kpoints.automatic_density(
                    struct, calc.kpts["max_mixed_density"][1], force_gamma
                )
                if np.product(pmg_kpts1.kpts[0]) >= np.product(pmg_kpts2.kpts[0]):
                    pmg_kpts = pmg_kpts1
                else:
                    pmg_kpts = pmg_kpts2
            elif "reciprocal_density" in calc.kpts:
                pmg_kpts = Kpoints.automatic_density_by_vol(
                    struct, calc.kpts["reciprocal_density"], force_gamma
                )
            elif "grid_density" in calc.kpts:
                pmg_kpts = Kpoints.automatic_density(
                    struct, calc.kpts["grid_density"], force_gamma
                )
            elif "length_density" in calc.kpts:
                pmg_kpts = Kpoints.automatic_density_by_length(
                    struct, calc.kpts["length_density"], force_gamma
                )
            else:
                raise ValueError(f"Unsupported k-point generation scheme: {kpts}.")

            kpts = pmg_kpts.kpts[0]
            if pmg_kpts.style.name.lower() == "gamma":
                gamma = True
            else:
                gamma = False
            calc.set(kpts, gamma)

    # Handle INCAR swaps as needed
    if incar_copilot:

        if any(atoms.get_atomic_numbers() > 56):
            if verbose:
                warnings.warn(
                    "Copilot: Setting LMAXMIX = 6 because you have an f element."
                )
            calc.set(lmaxmix=6)
        elif any(atoms.get_atomic_numbers() > 20):
            if verbose:
                warnings.warn(
                    "Copilot: Setting LMAXMIX = 4 because you have a d element"
                )
            calc.set(lmaxmix=4)

        if (
            calc.asdict()["inputs"].get("luse_vdw", False) is True
            or calc.asdict()["inputs"].get("lhfcalc", False) is True
            or calc.asdict()["inputs"].get("ldau", False) is True
            or calc.asdict()["inputs"].get("metagga", None) is not None
        ) and calc.asdict()["inputs"].get("lasph", False) is False:
            if verbose:
                warnings.warn(
                    "Copilot: Setting LASPH = True because you have a +U, vdW, meta-GGA, or hybrid calculation."
                )
            calc.set(lasph=True)

        if (
            calc.asdict()["inputs"].get("lasph", False) is True
            and calc.asdict()["inputs"].get("lmaxtau", 6) < 8
        ):
            if verbose:
                warnings.warn(
                    "Copilot: Setting LMAXTAU = 8 because you have LASPH = True and an f element."
                )
            calc.set(lmaxtau=8)

        if (
            calc.asdict()["inputs"].get("lhfcalc", False) is True
            or calc.asdict()["inputs"].get("metagga", None) is not None
        ) and calc.asdict()["inputs"].get("algo", "Normal") != "All":
            if verbose:
                warnings.warn(
                    "Copilot: Setting ALGO = All because you have a meta-GGA or hybrid calculation."
                )
            calc.set(algo="All")

        if (
            calc.asdict()["inputs"].get("nedos", 301) > 301
            and calc.asdict()["inputs"].get("ismear", 1) != -5
            and calc.asdict()["inputs"].get("nsw", 0) == 0
            and np.product(calc.asdict()["inputs"].get("kpts", (1, 1, 1))) >= 4
        ):
            if verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = -5 because you have a static DOS calculation."
                )
            calc.set(ismear=-5)

        if (
            np.product(calc.asdict()["inputs"].get("kpts", (1, 1, 1))) < 4
            and calc.asdict()["inputs"].get("ismear", 1) == -5
        ):
            if verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = 0 because you don't have enough k-points for ISMEAR = -5."
                )
            calc.set(ismear=0)

        if (
            calc.asdict()["inputs"].get("kspacing", 0) > 0.5
            and calc.asdict()["inputs"].get("ismear", 1) == -5
        ):
            if verbose:
                warnings.warn(
                    "Copilot: KSPACING might be too large for ISMEAR = -5. Custodian should save you if needed, but you can also change ISMEAR to 0."
                )
            pass  # let Custodian deal with it

        if (
            calc.asdict()["inputs"].get("nsw", 0) > 0
            and calc.asdict()["inputs"].get("laechg", False) is True
        ):
            if verbose:
                warnings.warn(
                    "Copilot: Setting LAECHG = False because you have NSW > 0. LAECHG is not compatible with NSW > 0."
                )
            calc.set(laechg=False)

        if (
            calc.asdict()["inputs"].get("ldauprint", 0) == 0
            and calc.asdict()["inputs"].get("ldau", False) is True
        ):
            if verbose:
                warnings.warn("Copilot: Setting LDAUPRINT = 1 because LDAU = True.")
            calc.set(laechg=False)

        if (
            calc.asdict()["inputs"].get("lreal", False) is True
            and calc.asdict()["inputs"].get("nsw", 0) <= 1
        ):
            if verbose:
                warnings.warn(
                    "Copilot: Setting LREAL = False because you are running a static calculation. LREAL != False can be bad for energies."
                )
            calc.set(lreal=False)

        if calc.asdict()["inputs"].get("lorbit", None) is None and (
            calc.asdict()["inputs"].get("ispin", 1) == 2
            or np.any(atoms.get_initial_magnetic_moments() != 0)
        ):
            if verbose:
                warnings.warn(
                    "Copilot: Setting LORBIT = 11 because you have a spin-polarized calculation."
                )
            calc.set(lorbit=11)
