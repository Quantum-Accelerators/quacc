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
        version="0.0.6",
        packages=find_packages(),
        license="modified BSD",
        keywords="high-throughput automated workflow dft",
        data_files=["LICENSE.md"],
        zip_safe=False,
        install_requires=[
            "ase @ git+https://gitlab.com/ase/ase.git",
            "monty==2023.4.10",
            "pymatgen==2023.3.23",
            "emmet-core==0.51.13",
            "atomate2 @ git+https://github.com/arosen93/atomate2.git@molecule",
            "cclib==1.7.2",
            "numpy",
        ],
        extras_require={
            "covalent": [
                "covalent==0.209.1",
                "covalent-slurm-plugin==0.8.0",
                "covalent-ssh-plugin==0.17.0",
                "covalent-aws-plugins[all]==0.13.0"
            ],
            "jobflow": ["fireworks==2.0.2", "jobflow==0.1.11"],
            "vasp": ["custodian==2023.3.10"],
            "xtb": ["tblite[ase]==0.3.0", "xtb==22.1"],
            "docs": [
                "sphinx==6.1.3",
                "furo==2023.3.27",
                "m2r2==0.3.3.post2",
                "ipython==8.12.0",
                "nbsphinx==0.9.1",
                "nbsphinx-link==1.3.0",
                "autodoc_pydantic==1.8.0",
            ],
            "dev": ["pytest==7.3.1", "black==23.3.0", "isort==5.12.0"],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
