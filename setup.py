#!/usr/bin/env python3
"""
Setup script for mdCATH
"""

from setuptools import setup, find_packages

setup(
    name="mdcath",
    version="0.1.0",
    description="Process mdCATH dataset for ML applications",
    author="Biochemistry Team",
    author_email="info@example.com",
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=">=3.7",
    install_requires=[
        "h5py",
        "numpy",
        "pandas",
        "biopython",
        "pyyaml",
        "matplotlib",
        "seaborn",
        "tqdm",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
