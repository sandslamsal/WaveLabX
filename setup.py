from setuptools import setup, find_packages

setup(
    name="wavelabx",
    version="0.1.0",
    description="WaveLabX: A Python toolkit for wave-probe analysis",
    packages=find_packages(exclude=("tests",)),
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "pandas>=1.3",
        "matplotlib>=3.4",
    ],
)
