from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase.calculators.espresso import Espresso as Espresso_
from ase.calculators.espresso import EspressoProfile
from ase.calculators.espresso import EspressoTemplate as EspressoTemplate_
from ase.io.espresso import construct_namelist

from quacc import SETTINGS
from quacc.calculators.espresso.io import read, write
from quacc.calculators.espresso.keys import ALL_KEYS
from quacc.calculators.espresso.utils import parse_pp_and_cutoff
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import load_yaml_calc

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms


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
        write(
            Path(directory) / self.inputname,
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
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults

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

        super().__init__(profile=profile, parallel_info=parallel_info, **kwargs)

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
            self.calc_defaults["input_data"], keys=keys
        )

        if self.preset:
            config = load_yaml_calc(SETTINGS.ESPRESSO_PRESET_DIR / f"{self.preset}")
            preset_pp = parse_pp_and_cutoff(config, self.input_atoms)
            kwargs = merge_dicts(preset_pp, kwargs)
        return merge_dicts(self.calc_defaults, kwargs)
