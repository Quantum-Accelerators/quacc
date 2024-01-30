"""Custom Espresso calculator and template."""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase import Atoms
from ase.calculators.espresso import Espresso as Espresso_
from ase.calculators.espresso import EspressoProfile
from ase.calculators.espresso import EspressoTemplate as EspressoTemplate_
from ase.io import read, write
from ase.io.espresso import (
    Namelist,
    read_espresso_ph,
    write_espresso_ph,
    write_fortran_namelist,
)

from quacc import SETTINGS
from quacc.calculators.espresso.utils import get_pseudopotential_info, sanity_checks
from quacc.utils.dicts import recursive_dict_merge
from quacc.utils.files import load_yaml_calc

if TYPE_CHECKING:
    from typing import Any


class EspressoTemplate(EspressoTemplate_):
    """This is a wrapper around the ASE Espresso template that allows for the use of
    other binaries such as pw.x, ph.x, cp.x, etc."""

    def __init__(
        self, binary: str = "pw", test_run: bool = False, autorestart: bool = False
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

        Returns
        -------
        None
        """
        super().__init__()

        self.inputname = f"{binary}.in"
        self.outputname = f"{binary}.out"

        self.binary = binary

        self.outdirs = {
            "outdir": os.environ.get("ESPRESSO_TMPDIR", "."),
            "wfcdir": os.environ.get("ESPRESSO_TMPDIR", "."),
        }

        self.outfiles = {"fildos": "pwscf.dos"}

        self.test_run = test_run

        self.nruns = 0
        self.autorestart = autorestart

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
        parameters = sanity_checks(parameters, binary=self.binary)

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
        elif self.binary == "ph":
            with Path.open(directory / self.inputname, "w") as fd:
                write_espresso_ph(fd=fd, properties=properties, **parameters)
        else:
            with Path.open(directory / self.inputname, "w") as fd:
                write_fortran_namelist(
                    fd, binary=self.binary, properties=properties, **parameters
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
    def _test_run(parameters: dict[str, Any], directory: Path) -> dict[str, Any]:
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

    def read_results(self, directory: Path | str) -> dict[str, Any]:
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

        if self.binary == "pw":
            atoms = read(directory / self.outputname, format="espresso-out")
            results = dict(atoms.calc.properties())
        elif self.binary == "ph":
            with Path.open(directory / self.outputname, "r") as fd:
                results = read_espresso_ph(fd)
        elif self.binary == "dos":
            fildos = self.outfiles["fildos"]
            with Path(fildos).open("r") as fd:
                lines = fd.readlines()
                fermi = float(re.search(r"-?\d+\.?\d*", lines[0])[0])
                dos = np.loadtxt(lines[1:])
            results = {fildos.name: {"dos": dos, "fermi": fermi}}
        else:
            results = {}

        if "energy" not in results:
            results["energy"] = None

        return results

    def _output_handler(
        self, parameters: dict[str, Any], directory: Path
    ) -> dict[str, Any]:
        """
        Function that handles the various output of espresso binaries. If they are
        relative, they are resolved against `directory`. In any other case the
        function will raise a ValueError. This is to avoid the user to use absolute
        paths that might lead to unexpected behaviour when using Quacc.

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
        input_data = parameters.get("input_data", {})

        all_out = {**self.outdirs, **self.outfiles}
        working_dir = Path(directory).expanduser().resolve()

        for key in all_out:
            path = Path(all_out[key]).expanduser().resolve()
            for section in input_data:
                if key in input_data[section]:
                    path = Path(input_data[section][key]).expanduser().resolve()
                    input_data[section][key] = path

            try:
                path.relative_to(working_dir)
            except ValueError as e:
                raise ValueError(
                    f"Cannot use {key}={path} because it is not a subpath of {working_dir}. When using Quacc please provide subpaths relative to the working directory."
                ) from e
            if key in self.outdirs:
                path.mkdir(parents=True, exist_ok=True)
                self.outdirs[key] = path
            elif key in self.outfiles:
                self.outfiles[key] = path

        parameters["input_data"] = input_data

        return parameters


class Espresso(Espresso_):
    """
    This is a wrapper around the ASE Espresso calculator that adjusts input_data
    parameters and allows for the use of presets.

    Templates are used to set the binary and input/output file names.
    """

    def __init__(
        self,
        input_atoms: Atoms | None = None,
        preset: str | None = None,
        parallel_info: dict[str, Any] | None = None,
        template: EspressoTemplate | None = None,
        profile: EspressoProfile | None = None,
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
        parallel_info
            parallel_info is a dictionary passed to the ASE Espresso calculator
            profile. It is used to specify prefixes for the command line arguments.
            See the ASE documentation for more details.
        template
            ASE calculator templace which can be used to specify which espresso
            binary will be used in the calculation. This is taken care of by recipe
            in most cases.
        profile
            ASE calculator profile which can be used to specify the location of
            the espresso binary and pseudopotential files. This is taken care of
            internally using quacc settings.
        **kwargs
            Additional arguments to be passed to the Espresso calculator. Takes all valid
            ASE calculator arguments, such as `input_data` and `kpts`. Refer to
            `ase.calculators.espresso.Espresso` for details. Note that the full input
            must be described; use `{"system":{"ecutwfc": 60}}` and not the `{"ecutwfc": 60}`
            short-hand.

        Returns
        -------
        None
        """
        self.input_atoms = input_atoms or Atoms()
        self.preset = preset
        self.parallel_info = parallel_info
        self.kwargs = kwargs
        self._user_calc_params = {}

        template = template or EspressoTemplate("pw")
        full_path = Path(
            SETTINGS.ESPRESSO_BIN_DIR, SETTINGS.ESPRESSO_BINARIES[template.binary]
        )
        self._bin_path = str(full_path)
        self._binary = template.binary
        self._cleanup_params()
        self._pseudo_path = (
            self._user_calc_params.get("input_data", {})
            .get("control", {})
            .get("pseudo_dir", str(SETTINGS.ESPRESSO_PSEUDO))
        )
        self.profile = profile or EspressoProfile(
            binary=self._bin_path,
            parallel_info=parallel_info,
            pseudo_dir=self._pseudo_path,
        )

        super().__init__(
            profile=self.profile,
            parallel_info=self.parallel_info,
            **self._user_calc_params,
        )

        self.template = template

    def _cleanup_params(self) -> None:
        """
        Function that handles the kwargs. It will merge the user-supplied kwargs with
        the preset values, using the former as priority.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.kwargs.get("directory"):
            raise ValueError("quacc does not support the directory argument.")

        self.kwargs["input_data"] = Namelist(self.kwargs.get("input_data"))
        self.kwargs["input_data"].to_nested(binary=self._binary, **self.kwargs)

        if self.preset:
            calc_preset = load_yaml_calc(
                SETTINGS.ESPRESSO_PRESET_DIR / f"{self.preset}"
            )
            calc_preset["input_data"] = Namelist(calc_preset.get("input_data"))
            calc_preset["input_data"].to_nested(binary=self._binary, **calc_preset)
            if "pseudopotentials" in calc_preset:
                ecutwfc, ecutrho, pseudopotentials = get_pseudopotential_info(
                    calc_preset["pseudopotentials"], self.input_atoms
                )
                calc_preset.pop("pseudopotentials", None)
                if "kpts" in self.kwargs:
                    calc_preset.pop("kspacing", None)
                if "kspacing" in self.kwargs:
                    calc_preset.pop("kpts", None)
                self._user_calc_params = recursive_dict_merge(
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
                self._user_calc_params = recursive_dict_merge(calc_preset, self.kwargs)
        else:
            self._user_calc_params = self.kwargs

        if self._user_calc_params.get("kpts") and self._user_calc_params.get(
            "kspacing"
        ):
            raise ValueError("Cannot specify both kpts and kspacing.")
