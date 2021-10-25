import os
import pathlib
from setuptools import setup, find_packages
import subprocess
from os import path

this_directory = path.abspath(path.dirname(__file__))

with open("requirements.txt") as f:
    requirements = f.read().splitlines()


setup(
    name="anndata2pagoda",
    description="Small converter to generate pagoda2 web object from anndata",
    version='0.1',
    setup_requires=["setuptools"],
    package_dir={"anndata2pagoda": "anndata2pagoda"},
    packages=find_packages(),
    package_data={'': ['*.r', '*.R']},
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': ['anndata2pagoda=anndata2pagoda.anndata2pagoda:main'],
    },
)
