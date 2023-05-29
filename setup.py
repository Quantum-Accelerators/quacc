from pathlib import Path

from setuptools import find_packages, setup

module_dir = Path(__file__).resolve().parent

long_description = Path(module_dir / "README.md").read_text()

if __name__ == "__main__":
    setup(
        name="quacc",
        description="A platform to enable high-throughput, database-driven quantum chemistry and computational materials science",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Andrew S. Rosen",
        author_email="asrosen93@gmail.com",
        url="https://github.com/quantum-accelerators/quacc",
        python_requires=">=3.8.0, <3.10",
        version="0.1.0",
        packages=find_packages(),
        license="BSD-3",
        keywords="high-throughput automated workflow dft",
        classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "Operating System :: Unix",
            "Operating System :: MacOS",
        ],
        data_files=["LICENSE.md"],
        zip_safe=False,
        install_requires=[
            "ase @ git+https://github.com/Quantum-Accelerators/rASE.git",  # waiting on >3.22.1 for regular ase
            "atomate2[cclib] @ git+https://github.com/materialsproject/atomate2.git",  # waiting on >0.0.10
            "covalent @ git+https://github.com/AgnostiqHQ/covalent.git",  # waiting on > 0.222.0
            "custodian>=2023.5.12",
            "pymatgen>=2023.5.10",
            "monty>=2023.4.10",
        ],
        extras_require={
            "fireworks": ["fireworks>=2.0.3"],
            "tblite": ["tblite[ase]>=0.3.0"],
            "dev": ["black>=23.3.0", "isort>=5.12.0", "pytest>=7.3.1"],
            "docs": [
                "autodoc_pydantic==1.8.0",
                "furo==2023.5.20",
                "ipython==8.13.2",
                "jsonschema[format]==4.17.3",
                "myst_parser==1.0.0",
                "numpydoc==1.5.0",
                "sphinx_design==0.4.1",
            ],
            "strict": [
                "ase @ git+https://github.com/Quantum-Accelerators/rASE.git",  # waiting on >3.22.1 for regular ase
                "atomate2[cclib] @ git+https://github.com/materialsproject/atomate2.git",  # waiting on >0.0.10
                "black==23.3.0",
                "cclib==1.7.2",
                "covalent @ git+https://github.com/AgnostiqHQ/covalent.git",  # waiting on > 0.222.0
                "custodian==2023.5.12",
                "emmet-core==0.55.2",
                "fireworks==2.0.3",
                "isort==5.12.0",
                "jobflow==0.1.11",
                "maggma==0.51.1",
                "monty==2023.5.8",
                "numpy==1.24.3",
                "pydantic==1.10.8",
                "pymatgen==2023.5.10",
                "pytest==7.3.1",
                "tblite[ase]==0.3.0",
            ],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
