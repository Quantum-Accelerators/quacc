"""A Q-Chem calculator built on Pymatgen and Custodian functionality."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.calculator import FileIOCalculator

from quacc.calculators.qchem.io import read_qchem, write_qchem
from quacc.calculators.qchem.params import cleanup_attrs, make_qc_input
from quacc.calculators.qchem.qchem_custodian import run_custodian

if TYPE_CHECKING:
    from typing import Any, ClassVar, Literal, TypedDict

    from ase.atoms import Atoms
    from numpy.typing import NDArray

    class Results(TypedDict, total=False):
        """
        Type hint for the `results` attribute in [quacc.calculators.qchem.qchem.QChem][].
        """

        energy: float  # electronic energy in eV
        forces: NDArray  # forces in eV/A
        hessian: NDArray  # Hessian in eV/A^2/amu
        enthalpy: float  # total enthalpy in eV
        entropy: float  # total entropy in eV/K
        qc_output: dict[
            str, Any
        ]  # Output from `pymatgen.io.qchem.outputs.QCOutput.data`
        qc_input: dict[
            str, Any
        ]  # Input from `pymatgen.io.qchem.inputs.QCInput.as_dict()`
        custodian: dict[str, Any]  # custodian.json file metadata


class QChem(FileIOCalculator):
    """Custom Q-Chem calculator built on Pymatgen and Custodian."""

    implemented_properties: ClassVar[list[str]] = [
        "energy",
        "forces",
        "hessian",
        "enthalpy",
        "entropy",
        "qc_output",
        "qc_input",
        "custodian",
    ]
    results: ClassVar[Results] = {}

    def __init__(
        self,
        atoms: Atoms,
        charge: int = 0,
        spin_multiplicity: int = 1,
        rem: dict | None = None,
        opt: dict[str, list[str]] | None = None,
        pcm: dict | None = None,
        solvent: dict | None = None,
        smx: dict | None = None,
        scan: dict[str, list] | None = None,
        van_der_waals: dict[str, float] | None = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
        plots: dict | None = None,
        nbo: dict | None = None,
        geom_opt: dict | None = None,
        cdft: list[list[dict]] | None = None,
        almo_coupling: list[list[tuple[int, int]]] | None = None,
        svp: dict | None = None,
        pcm_nonels: dict | None = None,
        qchem_dict_set_params: dict[str, Any] | None = None,
        **fileiocalculator_kwargs,
    ) -> None:
        """
        Initialize the Q-Chem calculator. Most of the input parameters here are used to
        create a `pymatgen.io.qchem.inputs.QCInput` object. See the documentation for
        that class for more information.

        Parameters
        ----------
        atoms
            The Atoms object to be used for the calculation.
        charge
            The total charge of the molecular system.
        spin_multiplicity
            The spin multiplicity of the molecular system.
        rem
            A dictionary of all the input parameters for the rem section of
            QChem input file. e.g. rem = {'method': 'rimp2', 'basis': '6-31*G++'}
        opt
            A dictionary of opt sections, where each opt section is a key and
            the corresponding values are a list of strings. Strings must be
            formatted as instructed by the QChem manual. The different opt
            sections are: CONSTRAINT, FIXED, DUMMY, and CONNECT e.g. opt =
            {"CONSTRAINT": ["tors 2 3 4 5 25.0", "tors 2 5 7 9 80.0"], "FIXED":
            ["2 XY"]}
        pcm
            A dictionary of the PCM section, defining behavior for use of the
            polarizable continuum model. e.g. pcm = {"theory": "cpcm",
            "hpoints": 194}
        solvent
            A dictionary defining the solvent parameters used with PCM. e.g.
            solvent = {"dielectric": 78.39, "temperature": 298.15}
        smx
            A dictionary defining solvent parameters used with the SMD method, a
            solvent method that adds short-range terms to PCM. e.g. smx =
            {"solvent": "water"}
        scan
            A dictionary of scan variables. Because two constraints of the same
            type are allowed (for instance, two torsions or two bond stretches),
            each TYPE of variable (stre, bend, tors) should be its own key in
            the dict, rather than each variable. Note that the total number of
            variable (sum of lengths of all lists) CANNOT be more than two. e.g.
            scan = {"stre": ["3 6 1.5 1.9 0.1"], "tors": ["1 2 3 4 -180 180
            15"]}
        van_der_waals
            A dictionary of custom van der Waals radii to be used when
            constructing cavities for the PCM model or when computing, e.g.
            Mulliken charges. They keys are strs whose meaning depends on the
            value of vdw_mode, and the values are the custom radii in angstroms.
        vdw_mode
            Method of specifying custom van der Waals radii - 'atomic' or
            'sequential'. In 'atomic' mode (default), dict keys represent the
            atomic number associated with each radius (e.g., 12 = carbon). In
            'sequential' mode, dict keys represent the sequential position of a
            single specific atom in the input structure.
        plots
            A dictionary of all the input parameters for the plots section of
            the QChem input file.
        nbo
            A dictionary of all the input parameters for the nbo section of the
            QChem input file.
        geom_opt
            A dictionary of input parameters for the geom_opt section of the
            QChem input file. This section is required when using the libopt3
            geometry optimizer.
        cdft
            A list of lists of dictionaries, where each dictionary represents a
            charge constraint in the cdft section of the QChem input file. Each
            entry in the main list represents one state (allowing for
            multi-configuration calculations using constrained density
            functional theory - configuration interaction (CDFT-CI). Each state
            is represented by a list, which itself contains some number of
            constraints (dictionaries).

            1. For a single-state calculation with two constraints:
            cdft=[[
                {"value": 1.0, "coefficients": [1.0], "first_atoms": [1],
                "last_atoms": [2], "types": [None]}, {"value": 2.0,
                "coefficients": [1.0, -1.0], "first_atoms": [1, 17],
                "last_atoms": [3, 19],
                    "types": ["s"]}
            ]]

            Note that a type of None will default to a charge constraint (which
            can also be accessed by requesting a type of "c" or "charge".

            2. For a multi-reference calculation:
            cdft=[
                [
                    {"value": 1.0, "coefficients": [1.0], "first_atoms": [1],
                    "last_atoms": [27],
                        "types": ["c"]},
                    {"value": 0.0, "coefficients": [1.0], "first_atoms": [1],
                    "last_atoms": [27],
                        "types": ["s"]},
                ], [
                    {"value": 0.0, "coefficients": [1.0], "first_atoms": [1],
                    "last_atoms": [27],
                        "types": ["c"]},
                    {"value": -1.0, "coefficients": [1.0], "first_atoms": [1],
                    "last_atoms": [27],
                        "types": ["s"]},
                ]
            ]
        almo_coupling
            A list of lists of int 2-tuples used for calculations of
            diabatization and state coupling calculations
            relying on the absolutely localized molecular orbitals (ALMO)
            methodology. Each entry in the main list represents a single
            state (two states are included in an ALMO calculation). Within a
            single state, each 2-tuple represents the charge and spin
            multiplicity of a single fragment.
            e.g. almo=[[(1, 2), (0, 1)], [(0, 1), (1, 2)]]
        svp
            TODO.
        pcm_nonels
            TODO.
        qchem_dict_set_params
            Keyword arguments to be passed to `pymatgen.io.qchem.sets.QChemDictSet`,
            which will generate a `QCInput`. If `qchem_dict_set_params` is specified,
            the resulting `QCInput` will be merged with the `QCInput` generated from
            the `QChem` calculator kwargs, with the former taking priority. Accepts
            all arguments that `pymatgen.io.qchem.sets.QChemDictSet` accepts, except
            for `molecule`, which will always be generated from `atoms`. By default,
            `job_type`, `basis_set`, and `scf_algorithm` will be pulled from the `rem`
            kwarg if not specified in `qchem_dict_set_params`. `qchem_version` will
            default to 6 if not specified in `qchem_dict_set_params`.
        **fileiocalculator_kwargs
            Additional arguments to be passed to
            [ase.calculators.calculator.FileIOCalculator][].

        Returns
        -------
        None
        """

        # Assign variables to self
        self.atoms = atoms
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.rem = rem or {}
        self.opt = opt
        self.pcm = pcm
        self.solvent = solvent
        self.smx = smx
        self.scan = scan
        self.van_der_waals = van_der_waals
        self.vdw_mode = vdw_mode
        self.plots = plots
        self.nbo = nbo
        self.geom_opt = geom_opt
        self.cdft = cdft
        self.almo_coupling = almo_coupling
        self.svp = svp
        self.pcm_nonels = pcm_nonels
        self.qchem_dict_set_params = qchem_dict_set_params or {}
        self.fileiocalculator_kwargs = fileiocalculator_kwargs

        # Instantiate previous orbital coefficients
        self.prev_orbital_coeffs = None

        # Clean up parameters
        cleanup_attrs(self)

        # Set default params
        self._set_default_params()

        # Instantiate the calculator
        super().__init__(
            restart=None,
            label=None,
            command="",
            atoms=self.atoms,
            profile=None,
            **self.fileiocalculator_kwargs,
        )

    def write_input(
        self,
        atoms: Atoms,
        properties: list[str] | None = None,
        system_changes: list[str] | None = None,
    ) -> None:
        """
        Write the Q-Chem input files.

        Parameters
        ----------
        atoms
            The Atoms object to be used for the calculation.
        properties
            List of properties to calculate.
        system_changes
            List of changes to the system since last calculation.

        Returns
        -------
        None
        """
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        qc_input = make_qc_input(self, atoms)

        write_qchem(
            qc_input,
            directory=self.directory,
            prev_orbital_coeffs=self.prev_orbital_coeffs,
        )

    def execute(self) -> int:
        """
        Execute Q-Chem.

        Returns
        -------
        int
            The return code.
        """

        run_custodian(directory=self.directory)
        return 0

    def read_results(self) -> None:
        """
        Read the Q-Chem output files. Update the .results and .prev_orbital_coeffs
        attributes.

        Returns
        -------
        None
        """
        results, prev_orbital_coeffs = read_qchem(directory=self.directory)
        self.results = results
        self.prev_orbital_coeffs = prev_orbital_coeffs

    def _set_default_params(self) -> None:
        """
        Store the parameters that have been passed to the Q-Chem calculator in
        FileIOCalculator's self.default_parameters.

        Returns
        -------
        None
        """
        params = {
            "charge": self.charge,
            "spin_multiplicity": self.spin_multiplicity,
            "rem": self.rem,
            "opt": self.opt,
            "pcm": self.pcm,
            "solvent": self.solvent,
            "smx": self.smx,
            "scan": self.scan,
            "van_der_waals": self.van_der_waals,
            "vdw_mode": self.vdw_mode,
            "plots": self.plots,
            "nbo": self.nbo,
            "geom_opt": self.geom_opt,
            "cdft": self.cdft,
            "almo_coupling": self.almo_coupling,
            "svp": self.svp,
            "pcm_nonels": self.pcm_nonels,
            "qchem_dict_set_params": self.qchem_dict_set_params,
        }

        self.default_parameters = {k: v for k, v in params.items() if v is not None}
