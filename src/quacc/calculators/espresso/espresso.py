from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from ase import Atoms
from ase.calculators.espresso import Espresso as Espresso_
from ase.calculators.espresso import EspressoProfile
from ase.calculators.espresso import EspressoTemplate as EspressoTemplate_
from ase.io.espresso import construct_namelist

from quacc import SETTINGS
from quacc.calculators.espresso.io import read, write
from quacc.calculators.espresso.keys import ALL_KEYS
from quacc.calculators.espresso.utils import parse_pw_preset
from quacc.utils.dicts import recursive_dict_merge
from quacc.utils.files import load_yaml_calc

if TYPE_CHECKING:
    from typing import Any


class EspressoTemplate(EspressoTemplate_):
    """
    This is a wrapper around the ASE Espresso template that allows for the use
    of other binaries such as pw.x, ph.x, cp.x, etc.
    """

    def __init__(self, binary: str = "pw") -> None:
        """
        Initialize the Espresso template.

        Parameters
        ----------
        binary
            The name of the espresso binary to use. This is used to set the
            input/output file names. By default we fall bacl on "pw".

        Returns
        -------
        None
        """
        super().__init__()

        self.inputname = f"{binary}.in"
        self.outputname = f"{binary}.out"

        self.binary = binary

        self.outdirs = {
            "outdir": os.environ.get("ESPRESSO_TMPDIR"),
            "wfcdir": os.environ.get("ESPRESSO_TMPDIR"),
        }

    def write_input(
        self,
        profile: EspressoProfile,
        directory: Path | str,
        atoms: Atoms,
        parameters: dict[str, Any],
        properties: Any,
    ) -> None:
        """
        The function that should be used instead of the one in ASE EspressoTemplate
        to write the input file. It calls a customly defined write function.

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
        self._outdir_handler(parameters, directory)

        write(
            directory / self.inputname,
            atoms,
            binary=self.binary,
            properties=properties,
            pseudo_dir=str(profile.pseudo_path),
            **parameters,
        )

    def read_results(self, directory: Path | str) -> dict[str, Any]:
        """
        The function that should be used instead of the one in ASE EspressoTemplate
        to read the output file. It calls a customly defined read function. It also
        adds the "energy" key to the results dictionnary if it is not present. This
        is needed if the calculation is not made with pw.x.

        Parameters
        ----------
        directory
            The directory in which to read the output file.

        Returns
        -------
        dict
            The results dictionnary
        """

        results = read(Path(directory) / self.outputname, binary=self.binary)
        if "energy" not in results:
            results["energy"] = None
        return results

    def _outdir_handler(
        self, parameters: dict[str, Any], directory: Path
    ) -> dict[str, Any]:
        """
        Function that handles the various outdir of espresso binaries. If they are relative,
        they are resolved against `directory`, which is the recommended approach.
        If the user-supplied paths are absolute, they are resolved and checked
        against `directory`, which is typically `os.getcwd()`. If they are not in `directory`,
        they will be ignored.

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

        for section in input_data:
            for d_key in self.outdirs.keys():
                if d_key in input_data[section]:
                    path = Path(input_data[section][d_key])
                    path = path.expanduser().resolve()
                    if directory.expanduser().resolve() not in path.parents:
                        self.outdirs[d_key] = path
                        continue
                    path.mkdir(parents=True, exist_ok=True)
                    input_data[section][d_key] = path

        self.outdirs = [path for path in self.outdirs.values() if path is not None]

        parameters["input_data"] = input_data

        return parameters


class Espresso(Espresso_):
    """
    This is a wrapper around the ASE Espresso calculator that adjusts input_data
    parameters and allows for the use of presets. Templates are used to set
    the binary and input/output file names.
    """

    def __init__(
        self,
        input_atoms: Atoms | None = None,
        preset: str | None = None,
        template: EspressoTemplate | None = None,
        profile: EspressoProfile | None = None,
        calc_defaults: dict[str, Any] | None = None,
        parallel_info: dict[str, Any] | None = None,
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
        profile
            ASE calculator profile which can be used to specify the location of
            the espresso binary and pseudopotential files. This is taken care of
            internally using quacc settings.
        calc_defaults
            A dictionary of default input_data parameters to pass to the Espresso
            calculator. These will be overridden by any user-supplied calculator
            **kwargs.
        parallel_info
            parallel_info is a dictionary passed to the ASE Espresso calculator
            profile. It is used to specify prefixes for the command line arguments.
            See the ASE documentation for more details.
        **kwargs
            Additional arguments to be passed to the Espresso calculator, e.g.
            `input_data`, `kpts`... Takes all valid ASE calculator arguments.

        Returns
        -------
        None
        """
        self.preset = preset
        self.input_atoms = input_atoms or Atoms()
        self.calc_defaults = calc_defaults or {}

        template = template or EspressoTemplate("pw")

        kwargs = self._kwargs_handler(template.binary, **kwargs)

        pseudo_path = (
            kwargs["input_data"]
            .get("control", {})
            .get("pseudo_dir", str(SETTINGS.ESPRESSO_PSEUDO))
        )

        bin_path = SETTINGS.ESPRESSO_BIN_PATHS[template.binary]
        profile = profile or EspressoProfile(
            binary=str(bin_path), parallel_info=parallel_info, pseudo_path=pseudo_path
        )

        if kwargs.get("directory"):
            raise ValueError("quacc does not support the directory argument.")

        super().__init__(
            profile=profile, directory=".", parallel_info=parallel_info, **kwargs
        )

        self.template = template

    def _kwargs_handler(self, binary: str, **kwargs) -> dict[str, Any]:
        """
        Function that handles the kwargs. It will merge the user-supplied
        kwargs with the defaults and preset values. Priority order is as follow:

        User-supplied kwargs > preset > defaults

        Parameters
        ----------
        binary
            The espresso binary used to construct the namelist
        **kwargs
            User-supplied kwargs

        Returns
        -------
        kwargs
            The merged kwargs
        """
        keys = ALL_KEYS[binary]
        kwargs["input_data"] = construct_namelist(kwargs.get("input_data"), keys=keys)
        self.calc_defaults["input_data"] = construct_namelist(
            self.calc_defaults.get("input_data"), keys=keys
        )

        kpts = kwargs.get("kpts")
        kspacing = kwargs.get("kspacing")

        if kpts and kspacing:
            raise ValueError("Cannot specify both kpts and kspacing.")

        if self.preset:
            config = load_yaml_calc(SETTINGS.ESPRESSO_PRESET_DIR / f"{self.preset}")
            preset = parse_pw_preset(config, self.input_atoms)
            kwargs = recursive_dict_merge(preset, kwargs)

        if kpts:
            kwargs.pop("kspacing", None)
        elif kspacing:
            kwargs.pop("kpts", None)

        kwargs = recursive_dict_merge(self.calc_defaults, kwargs)

        if kwargs.get("kpts") == "gamma":
            kwargs["kpts"] = None

        return kwargs
