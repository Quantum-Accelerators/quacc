"""A Q-Chem calculator built on Pymatgen and Custodian functionality"""
from __future__ import annotations

import inspect
from pathlib import Path
from typing import TYPE_CHECKING

from ase.calculators.calculator import FileIOCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.sets import QChemDictSet
from pymatgen.io.qchem.utils import lower_and_check_unique

from quacc.calculators.qchem import custodian
from quacc.calculators.qchem.io import read_qchem, write_qchem
from quacc.utils.dicts import remove_dict_nones

if TYPE_CHECKING:
    from typing import Any, ClassVar, Literal

    from ase import Atoms
    from pymatgen.core.structure import Molecule

    from quacc.calculators.qchem.io import Results


class QChem(FileIOCalculator):
    """
    Custom Q-Chem calculator built on Pymatgen and Custodian.
    """

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
        atoms: Atoms | list[Atoms] | Literal["read"],
        charge: int,
        spin_multiplicity: int,
        rem: dict,
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
        qchem_dict_set_kwargs: dict[str, Any] | None = None,
        **fileiocalculator_kwargs,
    ) -> None:
        """
        Initialize the Q-Chem calculator.

        Parameters
        ----------
        atoms
            The Atoms object to be used for the calculation.
        charge
            The total charge of the molecular system.
        spin_multiplicity
            The spin multiplicity of the molecular system.
        ... TODO
        **fileiocalculator_kwargs
            Additional arguments to be passed to
            `ase.calculators.calculator.FileIOCalculator`.

        Returns
        -------
        None
        """

        # Assign variables to self
        self.atoms = atoms
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.rem = rem
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
        self.qchem_dict_set_kwargs = qchem_dict_set_kwargs or {}
        self.fileiocalculator_kwargs = fileiocalculator_kwargs

        # Instantiate previous orbital coefficients
        self._prev_orbital_coeffs = None

        if "directory" in self.fileiocalculator_kwargs:
            raise NotImplementedError("The directory kwarg is not supported.")

        # Clean up parameters
        for attr in [
            "rem",
            "pcm",
            "solvent",
            "smx",
            "scan",
            "van_der_waals",
            "plots",
            "nbo",
            "geom_opt",
            "svp",
            "pcm_nonels",
        ]:
            attr_val = lower_and_check_unique(getattr(self, attr))
            setattr(self, attr, attr_val)

        self._molecule = self._get_molecule()
        self._set_default_params()

        # Get Q-Chem executable command
        self.command = self._manage_environment()

        # Instantiate the calculator
        FileIOCalculator.__init__(
            self,
            restart=None,
            ignore_bad_restart_file=FileIOCalculator._deprecated,
            label=None,
            atoms=self.atoms,
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

        # TODO: Merge both input sets if both are passed.
        if self.qchem_dict_set_kwargs:
            qc_input = QChemDictSet(self._molecule, **self.qchem_dict_set_kwargs)
        else:
            qc_input = QCInput(
                self._molecule,
                self.rem,
                opt=self.opt,
                pcm=self.pcm,
                solvent=self.solvent,
                smx=self.smx,
                scan=self.scan,
                van_der_waals=self.van_der_waals,
                vdw_mode=self.vdw_mode,
                plots=self.plots,
                nbo=self.nbo,
                geom_opt=self.geom_opt,
                cdft=self.cdft,
                almo_coupling=self.almo_coupling,
                svp=self.svp,
                pcm_nonels=self.pcm_nonels,
            )
        write_qchem(
            qc_input,
            prev_orbital_coeffs=self._prev_orbital_coeffs,
        )

    def read_results(self) -> None:
        """
        Read the Q-Chem output files. Update the .results and ._prev_orbital_coeffs
        attributes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        results, prev_orbital_coeffs = read_qchem()
        self.results = results
        self._prev_orbital_coeffs = prev_orbital_coeffs

    def _manage_environment(self) -> str:
        """
        Return the command to run the Q-Chem calculator via Custodian.

        Returns
        -------
        str
            The command flag to run Q-Chem with Custodian.
        """

        qchem_custodian_script = Path(inspect.getfile(custodian)).resolve()
        return f"python {qchem_custodian_script}"

    def _get_molecule(self) -> Molecule | list[Molecule] | Literal["read"]:
        """
        Clean up q-chem input parameters for the Q-Chem calculator.
        Modifies self.qchem_input_params in place.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # TODO: We should probably not be setting this here...
        if "scf_guess" not in self.rem:
            self.rem["scf_guess"] = "read"

        adaptor = AseAtomsAdaptor()

        if isinstance(self.atoms, Atoms):
            atoms_.charge = self.charge
            atoms_.spin_multiplicity = self.spin_multiplicity
            molecule = adaptor.get_molecule(self.atoms)
            return molecule
        if isinstance(self.atoms, list):
            molecule = []
            for atoms_ in self.atoms:
                atoms_.charge = self.charge
                atoms_.spin_multiplicity = self.spin_multiplicity
                molecule.append(adaptor.get_molecule(atoms_))
            return molecule
        if isinstance(self.atoms, str):
            return self.atoms

    def _set_default_params(self) -> None:
        """
        Store the parameters that have been passed to the Q-Chem
        calculator in FileIOCalculator's self.default_parameters.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.default_parameters = remove_dict_nones(
            {
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
            }
        )
