from pathlib import Path

from setuptools import find_packages, setup

module_dir = Path(__file__).resolve().parent

long_description = Path(module_dir / "README.md").read_text(encoding="utf8")

if __name__ == "__main__":
    setup(
        name="quacc",
        description="A platform to enable high-throughput, database-driven quantum chemistry and computational materials science",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Andrew S. Rosen",
        author_email="asrosen93@gmail.com",
        url="https://github.com/quantum-accelerators/quacc",
        python_requires=">=3.9.0, <3.11",
        version="0.1.0",
        packages=find_packages(),
        license="BSD-3",
        keywords="high-throughput automated workflow dft",
        classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "Operating System :: Unix",
            "Operating System :: MacOS",
        ],
        data_files=["LICENSE.md"],
        zip_safe=False,
        install_requires=[
            "ase @ https://gitlab.com/argon214/ase/-/archive/efermi/ase-efermi.zip",  # waiting on >3.22.1
            "atomate2[cclib] @ git+https://github.com/materialsproject/atomate2.git",  # waiting on >0.0.10
            "covalent>=0.224.0-rc.0",  # waiting on > 0.222.0
            "custodian",
            "emmet-core>=0.51.11",
            "monty>=2023.4.10",
            "numpy",
            "pydantic",
            "pymatgen>=2022.7.8",
        ],
        extras_require={
            "fireworks": ["fireworks"],
            "optimizers": ["sella"],
            "tblite": ["tblite[ase]"],
            "dev": ["black==23.3.0", "isort==5.12.0", "pytest==7.3.2"],
            "docs": [
                "autodoc_pydantic==1.8.0",
                "furo==2023.5.20",
                "ipython==8.14.0",
                "jsonschema[format]==4.17.3",
                "nbsphinx==0.9.2",
                "myst_parser==1.0.0",
                "numpydoc==1.5.0",
                "sphinx_design==0.4.1",
                "sphinx-copybutton==0.5.2",
            ],
            "strict": [
                "ase @ https://gitlab.com/argon214/ase/-/archive/efermi/ase-efermi.zip",  # waiting on >3.22.1
                "atomate2 @ git+https://github.com/materialsproject/atomate2.git",  # waiting on >0.0.10
                "cclib==1.7.2",
                "covalent==0.226.0rc0",  # waiting on > 0.222.0
                "custodian==2023.6.5",
                "emmet-core==0.55.5",
                "jobflow==0.1.11",
                "maggma==0.51.7",
                "monty==2023.5.8",
                "numpy==1.24.3",
                "pydantic==1.10.2",
                "pymatgen==2023.5.31",
            ],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
