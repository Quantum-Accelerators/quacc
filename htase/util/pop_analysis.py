from pymatgen.command_line.bader_caller import bader_analysis_from_path
from pymatgen.command_line.bader_caller import chargemol_analysis_from_path
import os


def run_bader(path=None):
    """
    Runs a Bader partial charge and spin moment analysis using the VASP
    output files in the given path. This function requires that `bader` or
    `bader.exe` is located in your PATH environment variable. See
    http://theory.cm.utexas.edu/henkelman/code/bader for the bader code.

    Args:
        path (str): The path where the VASP output files are located.
        Must include CHGCAR, AECCAR0, AECCAR2, and POTCAR files. These
        files can be gzip'd or not -- it doesn't matter.
            Default: None (current working directory).

    Returns:
        Dictionary containing the Bader analysis summary.

    """

    if path is None:
        path = os.getcwd()

    # Make sure files are present.
    for f in ["CHGCAR", "AECCAR0", "AECCAR2", "POTCAR"]:
        if not os.path.exists(os.path.join(path, f)) and not os.path.exists(
            os.path.join(path, f"{f}.gz")
        ):
            raise ValueError("Could not find {} in {}".format(f, path))

    # Run Bader analysis
    bader_stats = bader_analysis_from_path(path)

    # Store the partial charge, which is much more useful than the
    # raw charge and is more intuitive than the charge transferred.
    # An atom with a positive partial charge is cationic, whereas
    # an atom with a negative partial charge is anionic.
    bader_stats["partial_charge"] = [-c for c in bader_stats["charge_transfer"]]

    # Some cleanup of the returned dictionary
    if "magmom" in bader_stats:
        bader_stats["spin_moment"] = bader_stats["magmom"]
    bader_stats.pop("charge", None)
    bader_stats.pop("charge_transfer", None)
    bader_stats.pop("reference_used", None)
    bader_stats.pop("magmom", None)

    return bader_stats


def run_chargemol(path=None):
    """
    Runs a Chargemol (i.e. DDEC6 + CM5) analysis using the VASP output files
    in the given path. This function requires that the chargemol executable,
    given by the name `Chargemol_09_26_2017_linux_parallel`,
    `Chargemol_09_26_2017_linux_serial`, or `chargemol` is in the system PATH
    environment variable. See https://sourceforge.net/projects/ddec/files for
    the Chargemol code.

    Args:
        path (str): The path where the VASP output files are located.
        Must include CHGCAR, AECCAR0, AECCAR2, and POTCAR files. These
        files can be gzip'd or not -- it doesn't matter.
            Default: None (current working directory).

    Returns:
        Dictionary containing the Chargemol analysis summary.

    """
    if path is None:
        path = os.getcwd()

    # Make sure files are present.
    for f in ["CHGCAR", "AECCAR0", "AECCAR2", "POTCAR"]:
        if not os.path.exists(os.path.join(path, f)) and not os.path.exists(
            os.path.join(path, f"{f}.gz")
        ):
            raise ValueError("Could not find {} in {}".format(f, path))

    # Run Chargemol analysis
    chargemol_stats = chargemol_analysis_from_path(path)

    return chargemol_stats
