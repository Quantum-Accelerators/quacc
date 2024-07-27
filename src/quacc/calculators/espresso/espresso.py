"""Custom Espresso calculator and template."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atoms
from ase.calculators.espresso import EspressoProfile
from ase.calculators.espresso import EspressoTemplate as EspressoTemplate_
from ase.calculators.genericfileio import GenericFileIOCalculator
from ase.io import read, write
from ase.io.espresso import (
    Namelist,
    read_espresso_ph,
    write_espresso_ph,
    write_fortran_namelist,
)
from ase.io.espresso_namelist.keys import ALL_KEYS

from quacc import get_settings
from quacc.calculators.espresso.utils import (
    espresso_prepare_dir,
    get_pseudopotential_info,
    remove_conflicting_kpts_kspacing,
)
from quacc.utils.dicts import Remove, recursive_dict_merge, remove_dict_entries
from quacc.utils.files import load_yaml_calc, safe_decompress_dir

if TYPE_CHECKING:
    from typing import Any

LOGGER = logging.getLogger(__name__)


class EspressoTemplate(EspressoTemplate_):
    """
    A wrapper around the ASE Espresso template that allows for the use of
    other binaries such as pw.x, ph.x, cp.x, etc.
    """

    def __init__(
        self,
        binary: str = "pw",
        test_run: bool = False,
        autorestart: bool = False,
        outdir: str | Path | None = None,
    ) -> None:
        """
        Initialize the Espresso template.

        Parameters
        ----------
        binary
            The name of the espresso binary to use. This is used to set the
            input/output file names. By default we fall back to "pw".
        test_run
            If True, a test run is performed to check that the calculation
            input_data is correct or to generate some files/info if needed.
        autorestart
            If True, the calculation will automatically switch to 'restart'
            if this calculator performs more than one run. (ASE-relax/MD/NEB)
        outdir
            The directory that will be used as `outdir` in the input_data. If
            None, the directory will be set to the current working directory.

        Returns
        -------
        None
        """
        super().__init__()

        self.inputname = f"{binary}.in"
        self.outputname = f"{binary}.out"
        self.errorname = f"{binary}.err"
        self.binary = binary
        self._ase_known_binary = self.binary in ALL_KEYS
        self.test_run = test_run
        self.nruns = 0
        self.autorestart = autorestart
        self.outdir = outdir

    def write_input(
        self,
        profile: EspressoProfile,
        directory: Path | str,
        atoms: Atoms,
        parameters: dict[str, Any],
        properties: Any,
    ) -> None:
        """
        The function that should be used instead of the one in ASE EspressoTemplate to
        write the input file. It calls a customly defined write function.

        Parameters
        ----------
        profile
            The profile to use.
        directory
            The directory in which to write the input file.
        atoms
            The atoms object to use.
        parameters
            The parameters to use.
        properties
            Special ASE properties

        Returns
        -------
        None
        """
        directory = Path(directory)
        self._output_handler(parameters, directory)
        parameters = self._sanity_checks(parameters)

        if self.outdir:
            safe_decompress_dir(self.outdir)

        if self.test_run:
            self._test_run(parameters, directory)

        if self.binary == "pw":
            if self.autorestart and self.nruns > 0:
                parameters["input_data"]["electrons"]["startingpot"] = "file"
                parameters["input_data"]["electrons"]["startingwfc"] = "file"
            write(
                directory / self.inputname,
                atoms,
                format="espresso-in",
                pseudo_dir=str(profile.pseudo_dir),
                properties=properties,
                **parameters,
            )
        elif self.binary in ["ph", "phcg"]:
            with Path(directory, self.inputname).open(mode="w") as fd:
                write_espresso_ph(fd=fd, properties=properties, **parameters)
        else:
            with Path(directory, self.inputname).open(mode="w") as fd:
                write_fortran_namelist(
                    fd,
                    binary=self.binary if self._ase_known_binary else None,
                    properties=properties,
                    **parameters,
                )

    def execute(self, *args: Any, **kwargs: Any) -> None:
        super().execute(*args, **kwargs)
        self.nruns += 1

    @staticmethod
    def _search_keyword(parameters: dict[str, Any], key_to_search: str) -> str | None:
        """
        Function that searches for a keyword in the input_data.

        Parameters
        ----------
        parameters
            input_data, to search for the keyword

        Returns
        -------
        str
            The value of the keyword
        """
        input_data = parameters.get("input_data", {})

        for section in input_data:
            for key in input_data[section]:
                if key == key_to_search:
                    return input_data[section][key]
        return None

    @staticmethod
    def _test_run(parameters: dict[str, Any], directory: Path) -> None:
        """
        Almost all QE binaries will do a test run if a file named <prefix>.EXIT is
        present in the working directory. This function will create this file.

        Parameters
        ----------
        parameters
            input_data, which are needed to know the prefix
        directory
            The directory in which to write the EXIT file.

        Returns
        -------
        None
        """
        prefix = EspressoTemplate._search_keyword(parameters, "prefix") or "pwscf"

        Path(directory, f"{prefix}.EXIT").touch()

    def read_results(self, directory: os.PathLike) -> dict[str, Any]:
        """
        The function that should be used instead of the one in ASE EspressoTemplate to
        read the output file. It calls a customly defined read function. It also adds
        the "energy" key to the results dictionnary if it is not present. This is needed
        if the calculation is not made with pw.x.

        Parameters
        ----------
        directory
            The directory in which to read the output file.

        Returns
        -------
        dict
            The results dictionnary
        """
        results = {}
        if self.binary == "pw":
            atoms = read(Path(directory) / self.outputname, format="espresso-out")
            results = dict(atoms.calc.properties())
        elif self.binary in ["ph", "phcg"]:
            with Path(directory, self.outputname).open() as fd:
                results = read_espresso_ph(fd)
        elif self.binary == "dos":
            with Path(directory, "pwscf.dos").open() as fd:
                lines = fd.readlines()
                match = re.search(r"-?\d+\.?\d*", lines[0])
                fermi = float(match.group(0)) if match else None
                dos = np.loadtxt(lines[1:])
            results = {"dos_results": {"dos": dos, "fermi": fermi}}
        elif self.binary == "projwfc":
            with Path(directory, "pwscf.pdos_tot").open() as fd:
                lines = np.loadtxt(fd.readlines())
                energy = lines[1:, 0]
                dos = lines[1:, 1]
                pdos = lines[1:, 2]
            results = {"projwfc_results": {"energy": energy, "dos": dos, "pdos": pdos}}
        elif self.binary == "matdyn":
            fldos = Path(directory, "matdyn.dos")
            if fldos.exists():
                phonon_dos = np.loadtxt(fldos)
                results = {"matdyn_results": {"phonon_dos": phonon_dos}}

        if "energy" not in results:
            results["energy"] = None

        return results

    def _output_handler(
        self, parameters: dict[str, Any], directory: Path | str
    ) -> dict[str, Any]:
        """
        Function that handles the various output of espresso binaries. It will force the
        output directory and other output files to be set or deleted if needed.

        It will also prevent the user from setting environment variables that change the
        output directories.

        Parameters
        ----------
        parameters
            User-supplied kwargs
        directory
            The `directory` kwarg from the calculator.

        Returns
        -------
        dict[str, Any]
            The merged kwargs
        """
        os.environ.pop("ESPRESSO_TMPDIR", None)
        os.environ.pop("ESPRESSO_FILDVSCF_DIR", None)
        os.environ.pop("ESPRESSO_FILDRHO_DIR", None)

        espresso_outdir = Path(self.outdir or directory).expanduser().resolve()
        outkeys = espresso_prepare_dir(espresso_outdir, self.binary)

        input_data = parameters.get("input_data", {})
        input_data = recursive_dict_merge(input_data, outkeys, verbose=True)

        parameters["input_data"] = input_data

        return parameters

    def _sanity_checks(self, parameters: dict[str, Any]) -> dict[str, Any]:
        """
        Function that performs sanity checks on the input_data. It is meant
        to catch common mistakes that are not caught by the espresso binaries.

        Parameters
        ----------
        parameters
            The parameters dictionary which is assumed to already be in
            the nested format.

        Returns
        -------
        dict
            The modified dictionary parameters.
        """
        input_data = parameters.get("input_data", {})

        if self.binary == "pw":
            system = input_data.get("system", {})

            occupations = system.get("occupations", "fixed")
            smearing = system.get("smearing", None)
            degauss = system.get("degauss", None)

            if occupations == "fixed" and (smearing is not None or degauss is not None):
                LOGGER.warning(
                    "The occupations are set to 'fixed' but smearing or degauss is also set. This will be ignored."
                )
                system["smearing"] = Remove
                system["degauss"] = Remove

            parameters["input_data"]["system"] = system

        elif self.binary in ["ph", "phcg"]:
            input_ph = input_data.get("inputph", {})
            qpts = parameters.get("qpts", (0, 0, 0))

            qplot = input_ph.get("qplot", False)
            lqdir = input_ph.get("lqdir", False)
            recover = input_ph.get("recover", False)
            ldisp = input_ph.get("ldisp", False)

            is_grid = input_ph.get("start_q") or input_ph.get("start_irr")
            # Temporary patch for https://gitlab.com/QEF/q-e/-/issues/644
            if qplot and lqdir and recover and is_grid:
                prefix = input_ph.get("prefix", "pwscf")
                outdir = input_ph.get("outdir", ".")

                Path(outdir, "_ph0", f"{prefix}.q_1").mkdir(parents=True, exist_ok=True)
            if not (ldisp or qplot):
                if np.array(qpts).shape == (1, 4):
                    LOGGER.warning(
                        "qpts is a 2D array despite ldisp and qplot being set to False. Converting to 1D array"
                    )
                    qpts = tuple(qpts[0])
                if lqdir and is_grid and qpts != (0, 0, 0):
                    LOGGER.warning(
                        "lqdir is set to True but ldisp and qplot are set to False. The band structure will still be computed at each step. Setting lqdir to False"
                    )
                    input_ph["lqdir"] = False

            parameters["input_data"]["inputph"] = input_ph
            parameters["qpts"] = qpts

        return remove_dict_entries(parameters, remove_trigger=Remove)


class Espresso(GenericFileIOCalculator):
    """
    A wrapper around the ASE Espresso calculator that adjusts input_data
    parameters and allows for the use of presets.
    Templates are used to set the binary and input/output file names.
    """

    def __init__(
        self,
        input_atoms: Atoms | None = None,
        preset: str | None = None,
        template: EspressoTemplate | None = None,
        **kwargs,
    ) -> None:
        """
        Initialize the Espresso calculator.

        Parameters
        ----------
        input_atoms
            The input Atoms object to be used for the calculation.
        preset
            The name of a YAML file containing a list of parameters to use as
            a "preset" for the calculator. quacc will automatically look in the
            `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets),
            The .yaml extension is not necessary. Any user-supplied calculator
            **kwargs will override any corresponding preset values.
        template
            ASE calculator templace which can be used to specify which espresso
            binary will be used in the calculation. This is taken care of by recipe
            in most cases.
        **kwargs
            Additional arguments to be passed to the Espresso calculator. Takes all valid
            ASE calculator arguments, such as `input_data` and `kpts`. Refer to
            [ase.calculators.espresso.Espresso][] for details. Note that the full input
            must be described; use `{"system":{"ecutwfc": 60}}` and not the `{"ecutwfc": 60}`
            short-hand.

        Returns
        -------
        None
        """
        self.input_atoms = input_atoms or Atoms()
        self.preset = preset
        self.kwargs = kwargs
        self.user_calc_params = {}
        self._settings = get_settings()
        template = template or EspressoTemplate("pw")
        self._binary = template.binary
        full_path = Path(
            self._settings.ESPRESSO_BIN_DIR,
            self._settings.ESPRESSO_BINARIES[self._binary],
        )
        self._bin_path = str(full_path)

        if template._ase_known_binary:
            self._cleanup_params()
        else:
            LOGGER.warning(
                f"the binary you requested, `{self._binary}`, is not supported by ASE. This means that presets and usual checks will not be carried out, your `input_data` must be provided in nested format."
            )

            self.kwargs["input_data"] = Namelist(self.kwargs.get("input_data"))
            self.user_calc_params = self.kwargs

        self._pseudo_path = (
            self.user_calc_params.get("input_data", {})
            .get("control", {})
            .get("pseudo_dir", str(self._settings.ESPRESSO_PSEUDO))
        )

        profile = EspressoProfile(
            f"{self._settings.ESPRESSO_PARALLEL_CMD[0]} {self._bin_path} {self._settings.ESPRESSO_PARALLEL_CMD[1]}",
            self._pseudo_path,
        )

        super().__init__(
            template=template,
            profile=profile,
            directory=".",
            parameters=self.user_calc_params,
        )

    def _cleanup_params(self) -> None:
        """
        Function that handles the kwargs. It will merge the user-supplied kwargs with
        the preset values, using the former as priority.

        Returns
        -------
        None
        """
        if self.kwargs.get("directory"):
            raise NotImplementedError("quacc does not support the directory argument.")

        self.kwargs["input_data"] = Namelist(self.kwargs.get("input_data"))
        self.kwargs["input_data"].to_nested(binary=self._binary, **self.kwargs)

        if self.preset:
            calc_preset = load_yaml_calc(
                self._settings.ESPRESSO_PRESET_DIR / f"{self.preset}"
            )
            calc_preset["input_data"] = Namelist(calc_preset.get("input_data"))
            calc_preset["input_data"].to_nested(binary=self._binary, **calc_preset)
            if "pseudopotentials" in calc_preset:
                ecutwfc, ecutrho, pseudopotentials = get_pseudopotential_info(
                    calc_preset["pseudopotentials"], self.input_atoms
                )
                calc_preset.pop("pseudopotentials", None)
                calc_preset = remove_conflicting_kpts_kspacing(calc_preset, self.kwargs)
                self.user_calc_params = recursive_dict_merge(
                    calc_preset,
                    {
                        "input_data": {
                            "system": {"ecutwfc": ecutwfc, "ecutrho": ecutrho}
                        },
                        "pseudopotentials": pseudopotentials,
                    },
                    self.kwargs,
                )
            else:
                self.user_calc_params = recursive_dict_merge(calc_preset, self.kwargs)
        else:
            self.user_calc_params = self.kwargs

        if self.user_calc_params.get("kpts") is not None and self.user_calc_params.get(
            "kspacing"
        ):
            raise ValueError("Cannot specify both kpts and kspacing.")
