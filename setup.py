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
        python_requires=">=3.10.0",
        version="0.0.3",
        packages=find_packages(),
        license="modified BSD",
        keywords="high-throughput automated workflow dft vasp",
        data_files=["LICENSE.md"],
        zip_safe=False,
        install_requires=[
            "ase @ git+https://gitlab.com/ase/ase.git",
            "pymatgen==2022.3.7",
            "custodian==2022.2.13",
            "jobflow==0.1.6",
            "atomate2==0.0.6",
            "monty==2022.3.12",
            "FireWorks==2.0.2",
            "cclib==1.7.1",
        ],
        extras_require={
            "tests": ["pytest==7.1.0"],
            "docs": [
                "sphinx==4.4.0",
                "numpydoc==1.2",
                "m2r2==0.3.2",
                "ipython==8.1.1",
                "mistune==0.8.4",
                "pydata-sphinx-theme==0.8.0",
                "sphinx_panels==0.6.0",
                "autodoc_pydantic==1.6.1",
            ],
            "dev": ["pytest==7.1.0", "black==22.1.0"],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
