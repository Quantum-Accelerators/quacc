from setuptools import setup, find_packages
from pathlib import Path

module_dir = Path(__file__).resolve().parent
# with open("requirements.txt") as f:
#     required = f.read().splitlines()

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
        version="0.0.1",
        packages=find_packages(),
        license="modified BSD",
        keywords="high-throughput automated workflow dft vasp",
        data_files=["LICENSE.md"],
        zip_safe=False,
        # install_requires=required,
        extras_require={
            "tests": ["pytest==7.0.0"],
            "docs": [
                "sphinx==4.4.0",
                "numpydoc==1.2",
                "m2r2==0.3.2",
                "mistune==2.0.2",
                "pydata-sphinx-theme==0.8.0",
                "autodoc_pydantic==1.6.1",
                "sphinx_panels==0.6.0",
            ],
            "dev": ["pytest==7.0.0", "black==22.1.0"],
        },
        tests_require=["pytest"],
        include_package_data=True,
    )
