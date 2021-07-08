# -*- coding: utf-8 -*-

# Always prefer setuptools over distutils
import setuptools

# To use a consistent encoding
from codecs import open
from os import path

# The directory containing this file
HERE = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="freeda",
    version="0.0.1",
    author="Damian Dudka",
    author_email="damiandudka0@gmail.com",
    description="Finder of Rapidly Evolving Exons in De novo Assemblies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DDudka9/freeda",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS (10.14.6)",
    ],
    include_package_data=True,
    package_dir={"": "freeda"},
    packages=setuptools.find_packages(where="freeda"),
    python_requires="==3.7",
)