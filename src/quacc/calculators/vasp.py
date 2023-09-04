"""
A wrapper around ASE's Vasp calculator that makes it better suited for
high-throughput DFT.
"""
from __future__ import annotations

import inspect
import os
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.calculators.vasp import Vasp as Vasp_
from ase.calculators.vasp import setups as ase_setups
from ase.constraints import FixAtoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

from quacc import SETTINGS
from quacc.custodian import vasp as custodian_vasp
from quacc.utils.atoms import check_is_metal, set_magmoms
from quacc.utils.files import load_yaml_calc

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms


class Vasp(Vasp_):
    """
    This is a wrapper around the ASE Vasp calculator that adjusts INCAR
    parameters on-the-fly, allows for ASE to run VASP via Custodian, and
    supports several automatic k-point generation schemes from Pymatgen.

    Parameters
    ----------
    input_atoms
        The input Atoms object to be used for the calculation.
    preset
        The name of a YAML file containing a list of INCAR parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `VASP_PRESET_DIR` (default: quacc/presets/vasp) for the file, such that
        preset="BulkSet" is supported, for instance. The .yaml extension is not
        necessary. Any user-supplied calculator **kwargs will override any
        corresponding preset values.
    use_custodian
        Whether to use Custodian to run VASP. Default is True in settings.
    incar_copilot
        If True, the INCAR parameters will be adjusted if they go against the
        VASP manual. Default is True in settings.
    copy_magmoms
        If True, any pre-existing `atoms.get_magnetic_moments()` will be set in
        `atoms.set_initial_magnetic_moments()`. Set this to False if you want to
        use a preset's magnetic moments every time.
    preset_mag_default
        Default magmom value for sites without one explicitly specified in the
        preset. Only used if a preset is specified with an elemental_mags_dict
        key-value pair. Default is 1.0 in settings.
    mag_cutoff
        Set all initial magmoms to 0 if all have a magnitude below this value.
        Default is 0.05 in settings.
    elemental_magmoms
        A dictionary of elemental initial magnetic moments to pass to
        `quacc.utils.atoms.set_magmoms`, e.g. `{"Fe": 5, "Ni": 4}`.
    auto_kpts
        An automatic k-point generation scheme from Pymatgen. Options include:

        - {"line_density": float}. This will call
          `pymatgen.symmetry.bandstructure.HighSymmKpath`
            with `path_type="latimer_munro"`. The `line_density` value will be
            set in the `.get_kpoints` attribute.

        - {"kppvol": float}. This will call
          `pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_vol`
            with the given value for `kppvol`.

        - {"kppa": float}. This will call
          `pymatgen.io.vasp.inputs.Kpoints.automatic_density`
            with the given value for `kppa`.

        - {"length_densities": [float, float, float]}. This will call
          `pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_lengths`
            with the given value for `length_densities`.

        If multiple options are specified, the most dense k-point scheme will be
        chosen.
    auto_dipole
        If True, will automatically set dipole moment correction parameters
        based on the center of mass (in the c dimension by default).
    verbose
        If True, warnings will be raised when INCAR parameters are automatically
        changed. Default is True in settings.
    **kwargs
        Additional arguments to be passed to the VASP calculator, e.g.
        `xc='PBE'`, `encut=520`. Takes all valid ASE calculator arguments.

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
        elemental_magmoms: dict | None = None,
        auto_kpts: dict[Literal["line_density", "kppvol", "kppa"], float]
        | dict[Literal["length_densities"], list[float]] = None,
        auto_dipole: bool | None = None,
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
        self.elemental_magmoms = elemental_magmoms
        self.auto_kpts = auto_kpts
        self.auto_dipole = auto_dipole
        self.verbose = verbose
        self.kwargs = kwargs

        # Check constraints
        if (
            use_custodian
            and input_atoms.constraints
            and not all(isinstance(c, FixAtoms) for c in input_atoms.constraints)
        ):
            msg = "Atoms object has a constraint that is not compatible with Custodian. Set use_custodian = False."
            raise ValueError(msg)

        # Get VASP executable command, if necessary, and specify child
        # environment variables
        command = self._manage_environment()

        # Get user-defined preset parameters for the calculator
        if preset:
            calc_preset = load_vasp_yaml_calc(Path(SETTINGS.VASP_PRESET_DIR, preset))[
                "inputs"
            ]
        else:
            calc_preset = {}

        # Collect all the calculator parameters and prioritize the kwargs in the
        # case of duplicates.
        self.user_calc_params = calc_preset | kwargs
        none_keys = [k for k, v in self.user_calc_params.items() if v is None]
        for none_key in none_keys:
            del self.user_calc_params[none_key]

        # Allow the user to use setups='mysetups.yaml' to load in a custom
        # setups from a YAML file
        if (
            isinstance(self.user_calc_params.get("setups"), str)
            and self.user_calc_params["setups"] not in ase_setups.setups_defaults
        ):
            self.user_calc_params["setups"] = load_vasp_yaml_calc(
                Path(SETTINGS.VASP_PRESET_DIR, self.user_calc_params["setups"])
            )["inputs"]["setups"]

        # Handle special arguments in the user calc parameters that ASE does not
        # natively support
        if (
            self.user_calc_params.get("elemental_magmoms")
            and self.elemental_magmoms is None
        ):
            self.elemental_magmoms = self.user_calc_params["elemental_magmoms"]
        if self.user_calc_params.get("auto_kpts") and self.auto_kpts is None:
            self.auto_kpts = self.user_calc_params["auto_kpts"]
        if self.user_calc_params.get("auto_dipole") and self.auto_dipole is None:
            self.auto_dipole = self.user_calc_params["auto_dipole"]
        for k in {"elemental_magmoms", "auto_kpts", "auto_dipole"}:
            self.user_calc_params.pop(k, None)

        # Make automatic k-point mesh
        if self.auto_kpts and not self.user_calc_params.get("kpts"):
            self._convert_auto_kpts()

        # Add dipole corrections if requested
        if self.auto_dipole:
            self._set_auto_dipole()

        # Set magnetic moments
        set_magmoms(
            input_atoms,
            elemental_mags_dict=self.elemental_magmoms,
            copy_magmoms=copy_magmoms,
            elemental_mags_default=preset_mag_default,
            mag_cutoff=mag_cutoff,
        )

        # Handle INCAR swaps as needed
        if incar_copilot:
            self._calc_swaps()

        # Remove unused INCAR flags
        self._remove_unused_flags()

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

        # Check if Custodian should be used and confirm environment variables
        # are set
        if self.use_custodian:
            # Return the command flag
            run_vasp_custodian_file = Path(inspect.getfile(custodian_vasp)).resolve()
            return f"python {run_vasp_custodian_file}"

        if "ASE_VASP_COMMAND" not in os.environ and "VASP_SCRIPT" not in os.environ:
            warnings.warn(
                "ASE_VASP_COMMAND or VASP_SCRIPT must be set in the environment to run VASP. See the ASE Vasp calculator documentation for details.",
                UserWarning,
            )
        return None

    def _remove_unused_flags(self) -> None:
        """
        Removes unused flags in the INCAR, like EDIFFG if you are doing NSW = 0.

        Returns
        -------
        None
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

    def _set_auto_dipole(self) -> None:
        """
        Sets flags related to the auto_dipole kwarg.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        com = self.input_atoms.get_center_of_mass(scaled=True)
        if "dipol" not in self.user_calc_params:
            self.user_calc_params["dipol"] = com
        if "idipol" not in self.user_calc_params:
            self.user_calc_params["idipol"] = 3
        if "ldipol" not in self.user_calc_params:
            self.user_calc_params["ldipol"] = True

    def _calc_swaps(self) -> None:
        """
        Swaps out bad INCAR flags.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        is_metal = check_is_metal(self.input_atoms)
        calc = Vasp_(**self.user_calc_params)
        max_Z = self.input_atoms.get_atomic_numbers().max()

        if (
            not calc.int_params["lmaxmix"] or calc.int_params["lmaxmix"] < 6
        ) and max_Z > 56:
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
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ALGO = All because you have a hybrid calculation.",
                    UserWarning,
                )
            calc.set(algo="all")

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
            calc.int_params["ismear"] != -5
            and calc.int_params["nsw"] in (None, 0)
            and (
                np.prod(calc.kpts) >= 4
                or (
                    calc.float_params["kspacing"]
                    and calc.float_params["kspacing"] <= 0.5
                )
            )
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = -5 because you have a static calculation.",
                    UserWarning,
                )
            calc.set(ismear=-5)

        if (
            calc.int_params["ismear"] == -5
            and np.prod(calc.kpts) < 4
            and calc.float_params["kspacing"] is None
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = 0 because you don't have enough k-points for ISMEAR = -5.",
                    UserWarning,
                )
            calc.set(ismear=0)

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
            self.auto_kpts
            and self.auto_kpts.get("line_density", None)
            and calc.int_params["ismear"] != 0
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting ISMEAR = 0 and SIGMA = 0.01 because you are doing a line mode calculation.",
                    UserWarning,
                )
            calc.set(ismear=0, sigma=0.01)

        if calc.int_params["ismear"] == 0 and (
            not calc.float_params["sigma"] or calc.float_params["sigma"] > 0.05
        ):
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting SIGMA = 0.05 because ISMEAR = 0 was requested with SIGMA > 0.05.",
                    UserWarning,
                )
            calc.set(sigma=0.05)

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

        if calc.special_params["lreal"] and len(self.input_atoms) < 30:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting LREAL = False because you have a small system (< 30 atoms/cell).",
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
            calc.set(ncore=1, npar=None)

        if (
            (calc.int_params["ncore"] and calc.int_params["ncore"] > 1)
            or (calc.int_params["npar"] and calc.int_params["npar"] > 1)
        ) and len(self.input_atoms) <= 4:
            if self.verbose:
                warnings.warn(
                    "Copilot: Setting NCORE = 1 because you have a very small structure.",
                    UserWarning,
                )
            calc.set(ncore=1, npar=None)

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
            calc.set(npar=1, ncore=None)

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

        self.user_calc_params = calc.parameters

    def _convert_auto_kpts(
        self,
    ) -> None:
        """
        Shortcuts for pymatgen k-point generation schemes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        struct = AseAtomsAdaptor.get_structure(self.input_atoms)

        if self.auto_kpts.get("line_density", None):
            # TODO: Support methods other than latimer-munro
            kpath = HighSymmKpath(
                struct,
                path_type="latimer_munro",
                has_magmoms=np.any(struct.site_properties.get("magmom", None)),
            )
            kpts, _ = kpath.get_kpoints(
                line_density=self.auto_kpts["line_density"], coords_are_cartesian=True
            )
            kpts = np.stack(kpts)
            reciprocal = True
            gamma = None

        else:
            reciprocal = None
            force_gamma = self.user_calc_params.get("gamma", False)
            max_pmg_kpts = None
            for k, v in self.auto_kpts.items():
                if k == "kppvol":
                    pmg_kpts = Kpoints.automatic_density_by_vol(
                        struct,
                        v,
                        force_gamma=force_gamma,
                    )
                elif k == "kppa":
                    pmg_kpts = Kpoints.automatic_density(
                        struct,
                        v,
                        force_gamma=force_gamma,
                    )
                elif k == "length_densities":
                    pmg_kpts = Kpoints.automatic_density_by_lengths(
                        struct,
                        v,
                        force_gamma=force_gamma,
                    )
                else:
                    msg = f"Unsupported k-point generation scheme: {self.auto_kpts}."
                    raise ValueError(msg)

                max_pmg_kpts = (
                    pmg_kpts
                    if (
                        not max_pmg_kpts
                        or np.prod(pmg_kpts.kpts[0]) >= np.prod(max_pmg_kpts.kpts[0])
                    )
                    else max_pmg_kpts
                )

            kpts = max_pmg_kpts.kpts[0]
            gamma = max_pmg_kpts.style.name.lower() == "gamma"

        self.user_calc_params["kpts"] = kpts
        if reciprocal and self.user_calc_params.get("reciprocal") is None:
            self.user_calc_params["reciprocal"] = reciprocal
        if self.user_calc_params.get("gamma") is None:
            self.user_calc_params["gamma"] = gamma


def load_vasp_yaml_calc(yaml_path: str | Path) -> dict:
    """
    Loads a YAML file containing calculator settings. Used for VASP calculations
    and can read quacc-formatted YAMLs that are of the following format:
    ```yaml
    inputs:
      xc: pbe
      algo: all
      setups:
        Cu: Cu_pv
      elemental_magmoms:
        Fe: 5
        Cu: 1
    ```
    where `inputs` is a dictionary of ASE-style input parameters, `setups` is a
    dictionary of ASE-style pseudopotentials, and and `elemental_magmoms` is a
    dictionary of element-wise initial magmoms.

    Parameters
    ----------
    yaml_path
        Path to the YAML file. This function will look in the `VASP_PRESET_DIR`
        (default: quacc/presets/vasp) for the file, thereby assuming that
        `yaml_path` is a relative path within that folder.
    Returns
    -------

    dict
        The calculator configuration (i.e. settings).
    """

    config = load_yaml_calc(yaml_path)

    # Allow for either "Cu_pv" and "_pv" style setups
    if "inputs" in config:
        config["inputs"] = {
            k.lower(): v.lower() if isinstance(v, str) else v
            for k, v in config["inputs"].items()
        }
        for k, v in config["inputs"].get("setups", {}).items():
            if k in v:
                config["inputs"]["setups"][k] = v.split(k)[-1]

    return config
