from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms
from ase.calculators.onetep import Onetep as Onetep_
from ase.calculators.onetep import OnetepProfile
from ase.calculators.onetep import OnetepTemplate
from ase.io.onetep import read_onetep_in,read_onetep_out,write_onetep_in,get_onetep_keywords

from quacc import SETTINGS
from quacc.utils.dicts import merge_dicts
from quacc.utils.files import load_yaml_calc

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
        self.preset = preset
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults
        self.template = template
        self.profile = profile
        ## KWARGS should be a dict and since we only got 1 binary we do not need the handler
        self.kwargs = kwargs
        ### Onetepprof does not have the pseudo_path
        #pseudo_path = kwargs["pseudo_path"].get("pseudo_dir", str(SETTINGS.ESPRESSO_PP_PATH))


        super().__init__(profile=profile, parallel_info=parallel_info, **kwargs)

        self.template = template

