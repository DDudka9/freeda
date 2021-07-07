# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="freeda-ddudka_.-9",
    version="0.0.1",
    author="Damian Dudka",
    author_email="damiandudka0@gmail.com",
    description="Finder of Rapidly Evolving Exons in De novo Assemblies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DDudka9/freeda_package",
    project_urls={
        "FREEDA": "https://github.com/DDudka9/freeda_package",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS (10.14.6)",
    ],
    package_dir={"": "freeda"},
    packages=setuptools.find_packages(where="freeda"),
    python_requires=">=3.7",
)