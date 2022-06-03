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
            "pymatgen==2022.5.26",
            "custodian==2022.5.26",
            "jobflow==0.1.8",
            "atomate2==0.0.7",
            "monty==2022.4.26",
            "fireworks==2.0.3",
            "cclib==1.7.1",
        ],
        extras_require={
            "tests": ["pytest==7.1.2"],
            "docs": [
                "sphinx==4.5.0",
                "furo==2022.6.4",
                "m2r2==0.3.2",
                "ipython==8.4.0",
                "nbsphinx==0.8.8",
                "nbsphinx-link==1.3.0",
                "autodoc_pydantic==1.7.1",
            ],
            "dev": ["pytest==7.1.2", "black==22.3.0", "isort==5.10.1"],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
