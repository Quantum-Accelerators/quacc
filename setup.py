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
            "ase @ git+https://gitlab.com/ase/ase.git",
            "covalent @ git+https://github.com/arosen93/covalent.git@rosen-unique-workdir",
            "jobflow>=0.1.11",
            "pymatgen>=2023.5.10",
            "emmet-core>=0.55.1",
            "maggma>=0.50.4",
            "atomate2 @ git+https://github.com/materialsproject/atomate2.git",
            "cclib>=1.7.2",
            "monty",
            "numpy",
            "pydantic",
        ],
        extras_require={
            "fireworks": ["fireworks>=2.0.3"],
            "vasp": ["custodian>=2023.5.7"],
            "xtb": ["tblite[ase]>=0.3.0"],
            "docs": [
                "numpydoc==1.5.0",
                "ipython==8.13.2",
                "autodoc_pydantic==1.8.0",
                "myst_parser==1.0.0",
                "furo==2023.03.27",
                "jsonschema[format]",
                "sphinx_design==0.4.1",
            ],
            "dev": ["pytest==7.3.1", "black==23.3.0", "isort==5.12.0"],
            "strict": [
                "ase @ git+https://gitlab.com/ase/ase.git",
                "covalent @ git+https://github.com/arosen93/covalent.git@rosen-unique-workdir",
                # "covalent-slurm-plugin==0.16.0rc0",
                "monty==2023.5.8",
                "jobflow==0.1.11",
                "pymatgen==2023.5.10",
                "emmet-core==0.55.1",
                "maggma==0.50.4",
                "atomate2 @ git+https://github.com/materialsproject/atomate2.git",
                "custodian==2023.5.7",
                "cclib==1.7.2",
                "tblite[ase]==0.3.0",
                "numpy",
                "pydantic",
            ],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
