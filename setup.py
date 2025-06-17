import setuptools
import glob

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
    "atomicwrites==1.4.1",
    "attrs==25.1.0",
    "colorama==0.4.6",
    "iniconfig==2.0.0",
    "nmrglue==0.11",
    "numpy==2.2.3",
    "packaging==24.2",
    "pandas==2.2.3",
    "pluggy==1.5.0",
    "py==1.11.0",
    "pyparsing==3.2.1",
    "pytest==8.3.4",
    "python-dateutil==2.9.0.post0",
    "pytz==2025.1",
    "scipy==1.15.2",
    "six==1.17.0",
    "toml==0.10.2",
    "rdkit==2024.9.5",
    "requests==2.32.3",
    "patool==3.1.3",
    "Unidecode>=1.3.7",
]

setuptools.setup(
    name="apvalidation",
    version="0.5.76",
    author="liningtonlabs",
    author_email="liningtonlabstest@gmail.com",
    description="Testing installation of Package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/liningtonlab/apvalidation",
    license="MIT",
    package_data={
        "apvalidation": [
            "npmrd_data_exchange/standardization_files/*.json",
            "npmrd_data_exchange/validation/validator.py",
            "npmrd_data_exchange/standardization/standardizer.py",
            "extract/*.py",
        ],
    },
    include_package_data=True,
    packages=["apvalidation"],
    install_requires=requirements,
)
