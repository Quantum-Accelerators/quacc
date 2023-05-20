from pathlib import Path

from setuptools import find_packages, setup

module_dir = Path(__file__).resolve().parent

with open(module_dir / "README.md") as f:
    long_description = f.read()

if __name__ == "__main__":
    setup(
        name="quacc",
        description="A platform to enable high-throughput, database-driven quantum chemistry and computational materials science",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Andrew S. Rosen",
        author_email="asrosen93@gmail.com",
        url="https://github.com/arosen93/quacc",
        python_requires=">=3.8.0",
        version="0.1.0",
        packages=find_packages(),
        license="BSD-3",
        keywords="high-throughput automated workflow dft",
        data_files=["LICENSE.md"],
        zip_safe=False,
        install_requires=[
            "ase @ git+https://gitlab.com/ase/ase.git",  # waiting on >3.22.1
            "atomate2[cclib] @ git+https://github.com/materialsproject/atomate2.git",  # waiting on >0.0.10
            "covalent @ git+https://github.com/AgnostiqHQ/covalent.git@refs/pull/1641/head",  # waiting on PR
            "emmet-core>=0.55.1",
            "jobflow>=0.1.11",
            "maggma>=0.50.4",
            "pymatgen>=2023.5.10",
            "monty",
            "numpy",
            "pydantic",
        ],
        extras_require={
            "fireworks": ["fireworks>=2.0.3"],
            "vasp": ["custodian>=2023.5.7"],
            "tblite": ["tblite[ase]>=0.3.0"],
            "docs": [
                "autodoc_pydantic==1.8.0",
                "furo==2023.03.27",
                "ipython==8.13.2",
                "jsonschema[format]",
                "myst_parser==1.0.0",
                "numpydoc==1.5.0",
                "sphinx_design==0.4.1",
            ],
            "dev": ["black==23.3.0", "isort==5.12.0", "pytest==7.3.1"],
            "strict": [
                "ase @ git+https://gitlab.com/ase/ase.git",
                "atomate2[cclib] @ git+https://github.com/materialsproject/atomate2.git",
                "covalent @ git+https://github.com/AgnostiqHQ/covalent.git@refs/pull/1641/head",
                # "covalent-slurm-plugin==0.16.0rc0",
                "custodian==2023.5.12",
                "emmet-core==0.55.2",
                "jobflow==0.1.11",
                "maggma==0.50.4",
                "monty==2023.5.8",
                "pymatgen==2023.5.10",
                "tblite[ase]==0.3.0",
                "numpy",
                "pydantic",
            ],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
