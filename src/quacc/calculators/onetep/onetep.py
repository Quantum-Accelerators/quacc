from __future__ import annotations

from typing import TYPE_CHECKING
from quacc import SETTINGS
from ase import Atoms
from ase.calculators.onetep import Onetep as Onetep_
from ase.calculators.onetep import OnetepProfile, OnetepTemplate
from quacc.utils.dicts import merge_dicts
if TYPE_CHECKING:
    from typing import Any


class Onetep(Onetep_):
    def __init__(
        self,
        input_atoms: Atoms = None,
        preset: str = None,
        template: OnetepTemplate = None,
        profile: OnetepProfile = None,
        calc_defaults: dict[str | Any] = None,
        parallel_info: dict[str | Any] = None,
        **kwargs,
    ):
        template = OnetepTemplate("onetep",append=True)
        profile = OnetepProfile(SETTINGS.ONETEP_CMD,parallel_info = parallel_info)
        self.preset = preset
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults
        kwargs = merge_dicts(self.calc_defaults, kwargs)
        self.template = template
        super().__init__(profile=profile, parallel_info=parallel_info, **kwargs)

