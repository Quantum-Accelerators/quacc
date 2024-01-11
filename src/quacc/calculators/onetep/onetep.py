from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms
from ase.calculators.onetep import Onetep as Onetep_
from ase.calculators.onetep import OnetepProfile

from quacc import SETTINGS
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any


class Onetep(Onetep_):
    def __init__(
        self,
        input_atoms: Atoms = None,
        preset: str | None = None,
        profile: OnetepProfile = None,
        calc_defaults: dict[str | Any] | None = None,
        parallel_info: dict[str | Any] | None = None,
        **kwargs,
    ):
        profile = OnetepProfile(str(SETTINGS.ONETEP_CMD), parallel_info=parallel_info)
        self.preset = preset
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults
        kwargs = recursive_dict_merge(self.calc_defaults, kwargs)

        if kwargs.get("directory"):
            raise ValueError("quacc does not support the directory argument.")

        super().__init__(
            profile=profile, directory=".", parallel_info=parallel_info, **kwargs
        )
