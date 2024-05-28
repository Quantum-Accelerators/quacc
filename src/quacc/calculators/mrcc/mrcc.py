from __future__ import annotations

import re

from ase.calculators.genericfileio import (
    BaseProfile,
    CalculatorTemplate,
    GenericFileIOCalculator,
)

import quacc.calculators.mrcc.io as io


def get_version_from_mrcc_header(mrcc_header):
    match = re.search(r"Release date: (.*)$", mrcc_header, re.M)
    return match.group(1)


class MrccProfile(BaseProfile):
    def __init__(self, binary, **kwargs):
        """
        Parameters
        ----------
        binary : str
            Full path to the dmrcc binary, the path to the dmrcc_mpi binary must be specified if MRCC is to be run in parallel.
        """
        # Because MRCC handles its parallelization without being called with
        # mpirun/mpiexec/etc parallel should be set to False.
        # Whether or not it is run in parallel is controlled by mrccinput or mrccblocks
        super().__init__(parallel=False, parallel_info={})
        self.binary = binary

    def version(self):
        # XXX Allow MPI in argv; the version call should not be parallel.
        from ase.calculators.genericfileio import read_stdout

        stdout = read_stdout([self.binary, "does_not_exist"])
        return get_version_from_mrcc_header(stdout)

    def get_calculator_command(self, inputfile):
        return [self.binary, inputfile]


class MrccTemplate(CalculatorTemplate):
    _label = "mrcc"

    def __init__(self):
        super().__init__(
            "mrcc",
            implemented_properties=[
                "energy",
                "mp2_corr_energy",
                "ccsd_corr_energy",
                "ccsd(t)_corr_energy",
            ],
        )

        self.inputname = "MINP"
        self.outputname = f"{self._label}.out"
        self.errorname = f"{self._label}.err"

    def execute(self, directory, profile) -> None:
        profile.run(
            directory, self.inputname, self.outputname, errorfile=self.errorname
        )

    def write_input(self, profile, directory, atoms, parameters, properties):
        parameters = dict(parameters)

        mrccinput = {"calc": "PBE", "basis": "def2-SVP"}

        kw = {"charge": 0, "mult": 1, "mrccinput": mrccinput, "mrccblocks": ""}
        kw.update(parameters)

        io.write_mrcc(directory / self.inputname, atoms, kw)

    def read_results(self, directory):
        return io.read_mrcc_outputs(directory, directory / self.outputname)

    def load_profile(self, cfg, **kwargs):
        return MrccProfile.from_config(cfg, self.name, **kwargs)


class MRCC(GenericFileIOCalculator):
    """Class for doing MRCC calculations.

    Example:

      calc = MRCC(charge=0, mult=1, mrccinput={'calc':PBE, 'basis':'def2-SVP'},
        mrccblocks=' ')
    """

    def __init__(
        self,
        *,
        profile=None,
        directory=".",
        parallel_info=None,
        parallel=None,
        **kwargs,
    ):
        """Construct MRCC-calculator object.

        Parameters
        ==========
        charge: int

        mult: int

        mrccinput : dict[str,str]

        mrccblocks: str


        Examples
        ========
        Use default values:

        >>> from quacc.calculators.mrcc.mrcc import MRCC
        >>> MyMrccProfile = MrccProfile("/path/to/mrcc/dir/dmrcc")
        >>> h = Atoms(
        ...     "H",
        ...     calculator=MRCC(
        ...         profile=MyMrccProfile,
        ...         charge=0,
        ...         mult=1,
        ...         directory="water",
        ...         mrccinput={"calc": PBE, "basis": "def2-SVP"},
        ...         mrccblocks=" ",
        ...     ),
        ... )

        """

        assert (
            parallel is None
        ), "MRCC does not support keyword parallel - use mrccblocks or mrccinput"
        assert (
            parallel_info is None
        ), "MRCC does not support keyword parallel_info - use mrccblocks or mrccinput"

        super().__init__(
            template=MrccTemplate(),
            profile=profile,
            directory=directory,
            parameters=kwargs,
        )
