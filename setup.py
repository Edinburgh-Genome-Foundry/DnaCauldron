import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open("dnacauldron/version.py").read())  # loads __version__

setup(
    name="dnacauldron",
    version=__version__,
    author="Zulko",
    url="https://github.com/Edinburgh-Genome-Foundry/DnaCauldron",
    description="Cloning simulation for DNA assembly (Golden Gate, Gibson...)",
    long_description=open("pypi-readme.rst").read(),
    license="see LICENSE.txt",
    keywords="DNA assembly cloning simulator synthetic biology",
    scripts=["scripts/dnacauldron"],
    packages=find_packages(exclude="docs"),
    install_requires=[
        "Biopython",
        "numpy",
        "matplotlib",
        "fuzzywuzzy",
        "pandas",
        "scipy",
        "networkx",
        "dna_features_viewer",
        "flametree",
        "snapgene_reader",
        "proglog",
        "xlwt",
        "openpyxl",
        "python-Levenshtein",
        "xlrd",
    ],
)
