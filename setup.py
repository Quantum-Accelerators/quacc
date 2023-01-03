from pathlib import Path

from setuptools import find_packages, setup

module_dir = Path(__file__).resolve().parent

with open(module_dir / "README.md") as f:
    long_description = f.read()

if __name__ == "__main__":
    setup(
        name="quacc",
        description="Enhance ASE for high-throughput DFT",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Andrew S. Rosen",
        author_email="asrosen93@gmail.com",
        url="https://github.com/arosen93/quacc",
        python_requires=">=3.8.0",
        version="0.0.3",
        packages=find_packages(),
        license="modified BSD",
        keywords="high-throughput automated workflow dft",
        data_files=["LICENSE.md"],
        zip_safe=False,
        install_requires=[
            "ase @ git+https://gitlab.com/argon214/ase.git@rosen-dftb",
            "pymatgen==2022.11.7",
            "custodian==2022.5.26",
            "jobflow==0.1.9",
            "atomate2==0.0.8",
            "monty==2022.9.9",
            "fireworks==2.0.3",
            "cclib==1.7.2",
        ],
        extras_require={
            "tests": ["pytest==7.2.0"],
            "codes": ["xtb==22.1"],
            "docs": [
                "sphinx==6.0.0",
                "furo==2022.12.7",
                "m2r2==0.3.3",
                "ipython==8.8.0",
                "nbsphinx==0.8.11",
                "nbsphinx-link==1.3.0",
                "autodoc_pydantic==1.8.0",
            ],
            "dev": ["pytest==7.2.0", "black==22.12.0", "isort==5.11.4"],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
