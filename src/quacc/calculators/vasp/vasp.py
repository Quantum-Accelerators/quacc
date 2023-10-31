"""A wrapper around ASE's Vasp calculator that makes it better suited for high-throughput DFT."""
from __future__ import annotations

import inspect
import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.calculators.vasp import Vasp as Vasp_
from ase.calculators.vasp import setups as ase_setups
from ase.constraints import FixAtoms

from quacc.calculators.vasp import custodian
from quacc.calculators.vasp.io import load_vasp_yaml_calc
from quacc.calculators.vasp.params import (
    convert_auto_kpts,
    get_param_swaps,
    remove_unused_flags,
    set_auto_dipole,
)
from quacc.runners.prep import set_magmoms

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms


class Vasp(Vasp_):
    """
    This is a wrapper around the ASE Vasp calculator that adjusts INCAR parameters
    on-the-fly, allows for ASE to run VASP via Custodian, and supports several automatic
    k-point generation schemes from Pymatgen.
    """

    def __init__(
        self,
        input_atoms: Atoms,
        preset: None | str = None,
        use_custodian: bool | None = None,
        incar_copilot: Literal["off", "on", "aggressive"] | None = None,
        copy_magmoms: bool | None = None,
        preset_mag_default: float | None = None,
        mag_cutoff: None | float = None,
        elemental_magmoms: dict[str, float] | None = None,
        auto_kpts: dict[Literal["line_density", "kppvol", "kppa"], float]
        | dict[Literal["length_densities"], list[float]]
        | None = None,
        auto_dipole: bool | None = None,
        **kwargs,
    ) -> None:
        """
        Initialize the VASP calculator.

        Parameters
        ----------
        input_atoms
            The input Atoms object to be used for the calculation.
        preset
            The name of a YAML file containing a list of INCAR parameters to use as
            a "preset" for the calculator. quacc will automatically look in the
            `VASP_PRESET_DIR` (default: quacc/calculators/vasp/presets) for the
            file, such that preset="BulkSet" is supported, for instance. The .yaml
            extension is not necessary. Any user-supplied calculator **kwargs will
            override any corresponding preset values.
        use_custodian
            Whether to use Custodian to run VASP. Default is True in settings.
        incar_copilot
            Controls VASP co-pilot mode for automated INCAR parameter handling.
            Options include:
            off: Do not use co-pilot mode. INCAR parameters will be unmodified.
            on: Use co-pilot mode. This will only modify INCAR flags not already set by the user.
            aggressive: Use co-pilot mode in agressive mode. This will modify INCAR flags even if they are already set by the user.
        copilot_override
            If False, INCAR swaps enabled by the INCAR co-pilot will not override
            the user's chosen value. It will only override values that aren't set.
            Default is False in settings.
        copy_magmoms
            If True, any pre-existing `atoms.get_magnetic_moments()` will be set in
            `atoms.set_initial_magnetic_moments()`. Set this to False if you want to
            use a preset's magnetic moments every time. Default is True in settings.
        preset_mag_default
            Default magmom value for sites without one explicitly specified in the
            preset. Only used if a preset is specified with an elemental_mags_dict
            key-value pair. Default is 1.0 in settings.
        mag_cutoff
            Set all initial magmoms to 0 if all have a magnitude below this value.
            Default is 0.05 in settings.
        elemental_magmoms
            A dictionary of elemental initial magnetic moments to pass to
            [quacc.runners.prep.set_magmoms][], e.g. `{"Fe": 5, "Ni": 4}`.
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
        **kwargs
            Additional arguments to be passed to the VASP calculator, e.g.
            `xc='PBE'`, `encut=520`. Takes all valid ASE calculator arguments.

        Returns
        -------
        None
        """
        from quacc import SETTINGS

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
        self.kwargs = kwargs

        # Initialize for later
        self.user_calc_params = {}

        # Cleanup parameters
        self._cleanup_params()

        # Get VASP executable command, if necessary, and specify child
        # environment variables
        command = self._manage_environment()

        # Instantiate the calculator!
        super().__init__(
            atoms=self.input_atoms, command=command, **self.user_calc_params
        )

    def _manage_environment(self) -> str:
        """
        Manage the environment for the VASP calculator.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The command flag to pass to the Vasp calculator.
        """
        from quacc import SETTINGS

        # Set the VASP pseudopotential directory
        if SETTINGS.VASP_PP_PATH:
            os.environ["VASP_PP_PATH"] = str(SETTINGS.VASP_PP_PATH)

        # Set the ASE_VASP_VDW environmentvariable
        if SETTINGS.VASP_VDW:
            os.environ["ASE_VASP_VDW"] = str(SETTINGS.VASP_VDW)
        if self.user_calc_params.get("luse_vdw") and "ASE_VASP_VDW" not in os.environ:
            raise OSError(
                "VASP_VDW setting was not provided, yet you requested a vdW functional."
            )

        # Return Custodian executable command
        if self.use_custodian:
            run_vasp_custodian_file = Path(inspect.getfile(custodian)).resolve()
            return f"python {run_vasp_custodian_file}"

        # Return vanilla ASE command
        vasp_cmd = (
            SETTINGS.VASP_GAMMA_CMD
            if np.prod(self.user_calc_params.get("kpts", [1, 1, 1])) == 1
            else SETTINGS.VASP_CMD
        )

        return f"{SETTINGS.VASP_PARALLEL_CMD} {vasp_cmd}"

    def _cleanup_params(self) -> None:
        """
        Clean up various calculator attributes and parameters.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        from quacc import SETTINGS

        # Check constraints
        if (
            self.use_custodian
            and self.input_atoms.constraints
            and not all(isinstance(c, FixAtoms) for c in self.input_atoms.constraints)
        ):
            msg = "Atoms object has a constraint that is not compatible with Custodian."
            raise ValueError(msg)

        # Get user-defined preset parameters for the calculator
        if self.preset:
            calc_preset = load_vasp_yaml_calc(SETTINGS.VASP_PRESET_DIR / self.preset)[
                "inputs"
            ]
        else:
            calc_preset = {}

        # Collect all the calculator parameters and prioritize the kwargs in the
        # case of duplicates.
        self.user_calc_params = calc_preset | self.kwargs

        # Allow the user to use setups='mysetups.yaml' to load in a custom
        # setups from a YAML file
        if (
            isinstance(self.user_calc_params.get("setups"), (str, Path))
            and self.user_calc_params["setups"] not in ase_setups.setups_defaults
        ):
            self.user_calc_params["setups"] = load_vasp_yaml_calc(
                SETTINGS.VASP_PRESET_DIR / self.user_calc_params["setups"]
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
        for k in ["elemental_magmoms", "auto_kpts", "auto_dipole"]:
            self.user_calc_params.pop(k, None)

        # Make automatic k-point mesh
        if self.auto_kpts and not self.user_calc_params.get("kpts"):
            self.user_calc_params = convert_auto_kpts(
                self.user_calc_params, self.auto_kpts, self.input_atoms
            )

        # Add dipole corrections if requested
        if self.auto_dipole:
            self.user_calc_params = set_auto_dipole(
                self.user_calc_params, self.input_atoms
            )

        # Set magnetic moments
        self.input_atoms = set_magmoms(
            self.input_atoms,
            elemental_mags_dict=self.elemental_magmoms,
            copy_magmoms=self.copy_magmoms,
            elemental_mags_default=self.preset_mag_default,
            mag_cutoff=self.mag_cutoff,
        )

        # Handle INCAR swaps
        self.user_calc_params = get_param_swaps(
            self.user_calc_params,
            self.auto_kpts,
            self.input_atoms,
            self.incar_copilot,
        )

        # Remove unused INCAR flags
        self.user_calc_params = remove_unused_flags(self.user_calc_params)
