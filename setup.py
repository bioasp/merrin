from setuptools import setup, find_packages

setup(name="merrin",
    version = "0.1",
    description = "MEtabolic Regulation Rule INference",
    url = "https://github.com/bioas/merrin",
    author = "Kerian Thuillier",
    author_email = "kerian.thuillier@ens-rennes.fr",
    install_requires = [
        "bonesis",
        "clingo",
        "ginsim",
        "pandas",
        "pulp",
        "python-libsbml",
        "scipy",
    ],
    entry_points = {
        "console_scripts": [
            "merrin=merrin.cli:main",
            "bnet2flexflux=merrin.bnet2flexflux:main",
            "merrin-timeseries_generator=merrin.generators.timeseries:main",
        ],
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="computational systems biology",
    packages = find_packages()
)
