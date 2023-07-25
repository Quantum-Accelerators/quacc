"""
A wrapper around ASE's Vasp calculator that makes it better suited for high-throughput DFT.
"""
from __future__ import annotations

import inspect
import os
import warnings
from typing import Literal

import numpy as np
from ase import Atoms
from ase.calculators.vasp import Vasp as Vasp_
from ase.calculators.vasp import setups as ase_setups
from ase.constraints import FixAtoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

from quacc import SETTINGS
from quacc.custodian import vasp as custodian_vasp
from quacc.util.atoms import check_is_metal, set_magmoms
from quacc.util.files import load_yaml_calc


class Vasp(Vasp_):
    """
    This is a wrapper around the ASE Vasp calculator that adjusts INCAR parameters on-the-fly,
    allows for ASE to run VASP via Custodian, and supports several automatic k-point generation schemes
    from Pymatgen.

    Parameters
    ----------
    input_atoms
        The input Atoms object to be used for the calculation.
    preset
        The name of a YAML file containing a list of INCAR parameters to use as a "preset" for the calculator.
        quacc will automatically look in the `VASP_PRESET_DIR` (default: quacc/presets/vasp) for the file, such
        that preset="BulkSet" is supported, for instance. The .yaml extension is not necessary. Any user-suppplied
        calculator **kwargs will override any corresponding preset values.
    use_custodian
        Whether to use Custodian to run VASP.
        Default is True in settings.
    incar_copilot
        If True, the INCAR parameters will be adjusted if they go against the VASP manual.
        Default is True in settings.
    copy_magmoms
        If True, any pre-existing `atoms.get_magnetic_moments()` will be set in `atoms.set_initial_magnetic_moments()`.
        Set this to False if you want to use a preset's magnetic moments every time.
    preset_mag_default
        Default magmom value for sites without one explicitly specified in the preset. Only used if a preset is
        specified with an elemental_mags_dict key-value pair.
        Default is 1.0 in settings.
    mag_cutoff
        Set all initial magmoms to 0 if all have a magnitude below this value.
        Default is 0.05 in settings.
    verbose
        If True, warnings will be raised when INCAR parameters are automatically changed.
        Default is True in settings.
    **kwargs
        Additional arguments to be passed to the VASP calculator, e.g. `xc='PBE'`, `encut=520`. Takes all valid
        ASE calculator arguments, in addition to those custom to quacc.

    Returns
    -------
    Atoms
        The ASE Atoms object with attached VASP calculator.
    """

    def __init__(
        self,
        input_atoms: Atoms,
        preset: None | str = None,
        use_custodian: bool | None = None,
        incar_copilot: bool | None = None,
        copy_magmoms: bool | None = None,
        preset_mag_default: float | None = None,
        mag_cutoff: None | float = None,
        verbose: bool | None = None,
        **kwargs,
    ):
        # Set defaults
        use_custodian = (
            SETTINGS.VASP_USE_CUSTODIAN if use_custodian is None else use_custodian
        )
        incar_copilot = (
            SETTINGS.VASP_INCAR_COPILOT if incar_copilot is None else incar_copilot
        )
        copy_magmoms = (
            SETTINGS.VASP_COPY_MAGMOMS if copy_magmoms is None else copy_magmoms
        )
        preset_mag_default = (
            SETTINGS.VASP_PRESET_MAG_DEFAULT
            if preset_mag_default is None
            else preset_mag_default
        )
        mag_cutoff = SETTINGS.VASP_MAG_CUTOFF if mag_cutoff is None else mag_cutoff
        verbose = SETTINGS.VASP_VERBOSE if verbose is None else verbose

        # Assign variables to self
        self.input_atoms = input_atoms
        self.preset = preset
        self.use_custodian = use_custodian
        self.incar_copilot = incar_copilot
        self.copy_magmoms = copy_magmoms
        self.preset_mag_default = preset_mag_default
        self.mag_cutoff = mag_cutoff
        self.verbose = verbose
        self.kwargs = kwargs

        # Check constraints
        if (
            use_custodian
            and input_atoms.constraints
            and not all(isinstance(c, FixAtoms) for c in input_atoms.constraints)
        ):
            raise ValueError(
                "Atoms object has a constraint that is not compatible with Custodian. Set use_custodian = False."
            )

        # Get VASP executable command, if necessary, and specify child environment
        # variables
        command = self._manage_environment()

        # Get user-defined preset parameters for the calculator
        if preset:
            calc_preset = load_yaml_calc(
                os.path.join(SETTINGS.VASP_PRESET_DIR, preset)
            )["inputs"]
        else:
            calc_preset = {}

        # Collect all the calculator parameters and prioritize the kwargs
        # in the case of duplicates.
        self.user_calc_params = calc_preset | kwargs
        none_keys = [k for k, v in self.user_calc_params.items() if v is None]
        for none_key in none_keys:
            del self.user_calc_params[none_key]

        # Allow the user to use setups='mysetups.yaml' to load in a custom setups
        # from a YAML file
        if (
            isinstance(self.user_calc_params.get("setups"), str)
            and self.user_calc_params["setups"] not in ase_setups.setups_defaults
        ):
            self.user_calc_params["setups"] = load_yaml_calc(
                os.path.join(SETTINGS.VASP_PRESET_DIR, self.user_calc_params["setups"])
            )["inputs"]["setups"]

        # If the preset has auto_kpts but the user explicitly requests kpts, then
        # we should honor that.
        if kwargs.get("kpts") and calc_preset.get("auto_kpts"):
            del self.user_calc_params["auto_kpts"]

        # Handle special arguments in the user calc parameters that
        # ASE does not natively support
        if self.user_calc_params.get("elemental_magmoms"):
            elemental_mags_dict = self.user_calc_params["elemental_magmoms"]
        else:
            elemental_mags_dict = None
        if self.user_calc_params.get("auto_kpts"):
            auto_kpts = self.user_calc_params["auto_kpts"]
        else:
            auto_kpts = None
        if self.user_calc_params.get("auto_dipole"):
            auto_dipole = self.user_calc_params["auto_dipole"]
        else:
            auto_dipole = None
        self.user_calc_params.pop("elemental_magmoms", None)
        self.user_calc_params.pop("auto_kpts", None)
        self.user_calc_params.pop("auto_dipole", None)

        # Make automatic k-point mesh
        if auto_kpts:
            kpts, gamma, reciprocal = self._convert_auto_kpts(auto_kpts)
            self.user_calc_params["kpts"] = kpts
            if reciprocal and self.user_calc_params.get("reciprocal") is None:
                self.user_calc_params["reciprocal"] = reciprocal
            if self.user_calc_params.get("gamma") is None:
                self.user_calc_params["gamma"] = gamma

        # Add dipole corrections if requested
        if auto_dipole:
            com = input_atoms.get_center_of_mass(scaled=True)
            if "dipol" not in self.user_calc_params:
                self.user_calc_params["dipol"] = com
            if "idipol" not in self.user_calc_params:
                self.user_calc_params["idipol"] = 3
            if "ldipol" not in self.user_calc_params:
                self.user_calc_params["ldipol"] = True

        # Set magnetic moments
        set_magmoms(
            input_atoms,
            elemental_mags_dict=elemental_mags_dict,
            copy_magmoms=copy_magmoms,
            elemental_mags_default=preset_mag_default,
            mag_cutoff=mag_cutoff,
        )

        # Handle INCAR swaps as needed
        if incar_copilot:
            self.user_calc_params = self._calc_swaps(auto_kpts=auto_kpts)

        # Remove unused INCAR flags
        self.user_calc_params = self._remove_unused_flags()

        # Instantiate the calculator!
        super().__init__(atoms=input_atoms, command=command, **self.user_calc_params)

    def _manage_environment(self) -> str:
        """
        Manage the environment for the VASP calculator.

        Returns
        -------
        str
            The command flag to pass to the Vasp calculator.
        """

        # Check ASE environment variables
        if "VASP_PP_PATH" not in os.environ:
            warnings.warn(
                "The VASP_PP_PATH environment variable must point to the library of VASP pseudopotentials. See the ASE Vasp calculator documentation for details.",
                UserWarning,
            )

        # Check if Custodian should be used and confirm environment variables are set
        if self.use_custodian:
            # Return the command flag
            run_vasp_custodian_file = os.path.abspath(inspect.getfile(custodian_vasp))
            return f"python {run_vasp_custodian_file}"

        if "ASE_VASP_COMMAND" not in os.environ and "VASP_SCRIPT" not in os.environ:
            warnings.warn(
                "ASE_VASP_COMMAND or VASP_SCRIPT must be set in the environment to run VASP. See the ASE Vasp calculator documentation for details.",
                UserWarning,
            )
        return None

    def _remove_unused_flags(self) -> dict:
        """
        Removes unused flags in the INCAR, like EDIFFG if you are doing NSW = 0.

        Returns
        -------
        Dict
            Adjusted user-specified calculation parameters
        """

        if self.user_calc_params.get("nsw", 0) == 0:
            # Turn off opt flags if NSW = 0
            opt_flags = ("ediffg", "ibrion", "isif", "potim", "iopt")
            for opt_flag in opt_flags:
                self.user_calc_params.pop(opt_flag, None)

        if not self.user_calc_params.get(
            "ldau", False
        ) and not self.user_calc_params.get("ldau_luj"):
            # Turn off +U flags if +U is not even used
            ldau_flags = (
                "ldau",
                "ldauu",
                "ldauj",
                "ldaul",
                "ldautype",
                "ldauprint",
                "ldau_luj",
            )
            for ldau_flag in ldau_flags:
                self.user_calc_params.pop(ldau_flag, None)

        return self.user_calc_params

    def _calc_swaps(
        self,
        auto_kpts: None
        | dict[Literal["line_density", "reciprocal_density", "grid_density"], float]
        | dict[Literal["max_mixed_density"], list[float, float]]
        | dict[Literal["length_density"], list[float, float, float]],
    ) -> dict:
        """
        Swaps out bad INCAR flags.

        Parameters
        ----------
        auto_kpts
            The automatic k-point scheme dictionary

        Returns
        -------
        dict
            Dictionary of new user-specified calculation parameters
        """
        is_metal = check_is_metal(self.input_atoms)
        calc = Vasp_(**self.user_calc_params)
        max_Z = max(el.Z for el in self.input_atoms)

        if (
            not calc.int_params["lmaxmix"] or calc.int_params["lmaxmix"] < 6
        ) and max_Z > 57:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LMAXMIX = 6 because you have f electrons.",
                    UserWarning,
                )
            calc.set(lmaxmix=6)
        elif (
            not calc.int_params["lmaxmix"] or calc.int_params["lmaxmix"] < 4
        ) and max_Z > 20:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LMAXMIX = 4 because you have d electrons.",
                    UserWarning,
                )
            calc.set(lmaxmix=4)

        if (
            calc.bool_params["luse_vdw"]
            or calc.bool_params["lhfcalc"]
            or calc.bool_params["ldau"]
            or calc.dict_params["ldau_luj"]
            or calc.string_params["metagga"]
        ) and not calc.bool_params["lasph"]:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LASPH = True because you have a +U, vdW, meta-GGA, or hybrid calculation.",
                    UserWarning,
                )
            calc.set(lasph=True)

        if calc.string_params["metagga"] and (
            not calc.string_params["algo"]
            or calc.string_params["algo"].lower() != "all"
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ALGO = All because you have a meta-GGA calculation.",
                    UserWarning,
                )
            calc.set(algo="all")

        if calc.bool_params["lhfcalc"] and (
            not calc.string_params["algo"]
            or calc.string_params["algo"].lower() not in ["all", "damped"]
        ):
            if is_metal:
                calc.set(algo="damped", time=0.5)
                if self.verbose:
                    warnings.warn(
                        "Copilot: Setting ALGO = Damped, TIME = 0.5 because you have a hybrid calculation with a metal.",
                        UserWarning,
                    )
            else:
                calc.set(algo="all")
                if self.verbose:
                    warnings.warn(
                        "Copilot: Setting ALGO = All because you have a hybrid calculation.",
                        UserWarning,
                    )

        if (
            is_metal
            and (calc.int_params["ismear"] and calc.int_params["ismear"] < 0)
            and (calc.int_params["nsw"] and calc.int_params["nsw"] > 0)
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: You are relaxing a likely metal. Setting ISMEAR = 1 and SIGMA = 0.1.",
                    UserWarning,
                )
            calc.set(ismear=1, sigma=0.1)

        if (
            calc.int_params["nedos"]
            and calc.int_params["ismear"] != -5
            and calc.int_params["nsw"] in (None, 0)
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = -5 because you have a static DOS calculation.",
                    UserWarning,
                )
            calc.set(ismear=-5)

        if (
            calc.int_params["ismear"] == -5
            and np.product(calc.kpts) < 4
            and calc.float_params["kspacing"] is None
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = 0 because you don't have enough k-points for ISMEAR = -5.",
                    UserWarning,
                )
            calc.set(ismear=0)

        if (
            auto_kpts
            and auto_kpts.get("line_density", None)
            and calc.int_params["ismear"] != 0
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = 0 and SIGMA = 0.01 because you are doing a line mode calculation.",
                    UserWarning,
                )
            calc.set(ismear=0, sigma=0.01)

        if calc.int_params["ismear"] == -5 and (
            not calc.float_params["sigma"] or calc.float_params["sigma"] > 0.05
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting SIGMA = 0.05 because ISMEAR = -5 was requested with SIGMA > 0.05.",
                    UserWarning,
                )
            calc.set(sigma=0.05)

        if (
            calc.float_params["kspacing"]
            and calc.float_params["kspacing"] > 0.5
            and calc.int_params["ismear"] == -5
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: KSPACING is likely too large for ISMEAR = -5. Setting ISMEAR = 0.",
                    UserWarning,
                )
            calc.set(ismear=0)

        if (
            calc.int_params["nsw"]
            and calc.int_params["nsw"] > 0
            and calc.bool_params["laechg"]
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LAECHG = False because you have NSW > 0. LAECHG is not compatible with NSW > 0.",
                    UserWarning,
                )
            calc.set(laechg=False)

        if calc.int_params["ldauprint"] in (None, 0) and (
            calc.bool_params["ldau"] or calc.dict_params["ldau_luj"]
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LDAUPRINT = 1 because LDAU = True.", UserWarning
                )
            calc.set(ldauprint=1)

        if calc.special_params["lreal"] and calc.int_params["nsw"] in (None, 0, 1):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LREAL = False because you are running a static calculation. LREAL != False can be bad for energies.",
                    UserWarning,
                )
            calc.set(lreal=False)

        if not calc.int_params["lorbit"] and (
            calc.int_params["ispin"] == 2
            or np.any(self.input_atoms.get_initial_magnetic_moments() != 0)
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LORBIT = 11 because you have a spin-polarized calculation.",
                    UserWarning,
                )
            calc.set(lorbit=11)

        if (
            (calc.int_params["ncore"] and calc.int_params["ncore"] > 1)
            or (calc.int_params["npar"] and calc.int_params["npar"] > 1)
        ) and (
            calc.bool_params["lhfcalc"] is True
            or calc.bool_params["lrpa"] is True
            or calc.bool_params["lepsilon"] is True
            or calc.int_params["ibrion"] in [5, 6, 7, 8]
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting NCORE = 1 because NCORE/NPAR is not compatible with this job type.",
                    UserWarning,
                )
            calc.set(ncore=1)
            calc.set(npar=None)

        if (
            (calc.int_params["ncore"] and calc.int_params["ncore"] > 1)
            or (calc.int_params["npar"] and calc.int_params["npar"] > 1)
        ) and len(self.input_atoms) <= 4:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting NCORE = 1 because you have a very small structure.",
                    UserWarning,
                )
            calc.set(ncore=1)
            calc.set(npar=None)

        if (
            calc.int_params["kpar"]
            and calc.int_params["kpar"] > np.prod(calc.kpts)
            and calc.float_params["kspacing"] is None
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting KPAR = 1 because you have too few k-points to parallelize.",
                    UserWarning,
                )
            calc.set(kpar=1)

        if (
            calc.int_params["nsw"]
            and calc.int_params["nsw"] > 0
            and calc.int_params["isym"]
            and calc.int_params["isym"] > 0
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISYM = 0 because you are running a relaxation.",
                    UserWarning,
                )
            calc.set(isym=0)

        if calc.bool_params["lhfcalc"] is True and calc.int_params["isym"] in (1, 2):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISYM = 3 because you are running a hybrid calculation.",
                    UserWarning,
                )
            calc.set(isym=3)

        if calc.bool_params["lsorbit"]:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISYM = -1 because you are running an SOC calculation.",
                    UserWarning,
                )
            calc.set(isym=-1)

        if (
            (calc.int_params["ncore"] and calc.int_params["ncore"] > 1)
            or (calc.int_params["npar"] and calc.int_params["npar"] > 1)
        ) and (calc.bool_params["lelf"] is True):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting NPAR = 1 because NCORE/NPAR is not compatible with this job type.",
                    UserWarning,
                )
            calc.set(npar=1)
            calc.set(ncore=None)

        if not calc.string_params["efermi"]:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting EFERMI = MIDGAP per the VASP manual.", UserWarning
                )
            calc.set(efermi="midgap")

        if calc.bool_params["luse_vdw"] and "ASE_VASP_VDW" not in os.environ:
            warnings.warn(
                "ASE_VASP_VDW was not set, yet you requested a vdW functional.",
                UserWarning,
            )

        return calc.parameters

    def _convert_auto_kpts(
        self,
        auto_kpts: None
        | dict[Literal["line_density", "reciprocal_density", "grid_density"], float]
        | dict[Literal["max_mixed_density"], list[float, float]]
        | dict[Literal["length_density"], list[float, float, float]],
        force_gamma: bool = True,
    ) -> tuple[list[tuple[int, int, int]], None | bool, None | bool]:
        """
        Shortcuts for pymatgen k-point generation schemes.

        Parameters
        ----------
        auto_kpts
            Dictionary describing the automatic k-point scheme
        force_gamma
            Whether a gamma-centered mesh should be returned

        Returns
        -------
        List[int, int, int]
            List of k-points for use with the ASE Vasp calculator
        Optional[bool]
            The gamma command for use with the ASE Vasp calculator
        Optional[bool]
            The reciprocal command for use with the ASE Vasp calculator
        """
        struct = AseAtomsAdaptor.get_structure(self.input_atoms)

        if auto_kpts.get("line_density", None):
            # TODO: Support methods other than latimer-munro
            kpath = HighSymmKpath(struct, path_type="latimer_munro")
            kpts, _ = kpath.get_kpoints(
                line_density=auto_kpts["line_density"], coords_are_cartesian=True
            )
            kpts = np.stack(kpts)
            reciprocal = True
            gamma = None

        else:
            reciprocal = None
            if auto_kpts.get("max_mixed_density", None):
                if len(auto_kpts["max_mixed_density"]) != 2:
                    raise ValueError("Must specify two values for max_mixed_density.")

                if (
                    auto_kpts["max_mixed_density"][0]
                    > auto_kpts["max_mixed_density"][1]
                ):
                    warnings.warn(
                        "Warning: It is not usual that kppvol > kppa. Please make sure you have chosen the right k-point densities.",
                        UserWarning,
                    )
                pmg_kpts1 = Kpoints.automatic_density_by_vol(
                    struct, auto_kpts["max_mixed_density"][0], force_gamma=force_gamma
                )
                pmg_kpts2 = Kpoints.automatic_density(
                    struct, auto_kpts["max_mixed_density"][1], force_gamma=force_gamma
                )
                if np.product(pmg_kpts1.kpts[0]) >= np.product(pmg_kpts2.kpts[0]):
                    pmg_kpts = pmg_kpts1
                else:
                    pmg_kpts = pmg_kpts2
            elif auto_kpts.get("reciprocal_density", None):
                pmg_kpts = Kpoints.automatic_density_by_vol(
                    struct, auto_kpts["reciprocal_density"], force_gamma=force_gamma
                )
            elif auto_kpts.get("grid_density", None):
                pmg_kpts = Kpoints.automatic_density(
                    struct, auto_kpts["grid_density"], force_gamma=force_gamma
                )
            elif auto_kpts.get("length_density", None):
                if len(auto_kpts["length_density"]) != 3:
                    raise ValueError("Must specify three values for length_density.")
                pmg_kpts = Kpoints.automatic_density_by_lengths(
                    struct, auto_kpts["length_density"], force_gamma=force_gamma
                )
            else:
                raise ValueError(f"Unsupported k-point generation scheme: {auto_kpts}.")

            kpts = pmg_kpts.kpts[0]
            gamma = pmg_kpts.style.name.lower() == "gamma"

        return kpts, gamma, reciprocal
