from setuptools import setup, find_packages
from pathlib import Path

module_dir = Path(__file__).resolve().parent

with open(module_dir / "README.md") as f:
    long_description = f.read()

if __name__ == "__main__":
    setup(
        name="htase",
        description="Enhance ASE for high-throughput DFT",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Andrew S. Rosen",
        author_email="asrosen93@gmail.com",
        url="https://github.com/arosen93/HT-ASE",
        python_requires=">=3.7.0",
        version="0.0.1",
        packages=find_packages(),
        license="modified BSD",
        keywords="high-throughput automated workflow dft vasp",
        data_files=["LICENSE.md"],
        zip_safe=False,
        install_requires=[
            "ase>=3.22.1",
            "pymatgen>2022.0.17",
            "custodian>=2021.12.2",
            "jobflow>=0.1.6",
            "atomate2>=0.0.4",
        ],
        extras_require={"fireworks": ["fireworks>=1.9.8", "monty>2021.12.1"],},
        tests_require=["pytest"],
        include_package_data=True,
    )
