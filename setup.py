import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

version = {}
with open("dnacauldron/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="dnacauldron",
    version=version["__version__"],
    author="Zulko",
    url="https://github.com/Edinburgh-Genome-Foundry/DnaCauldron",
    description="Cloning simulation for DNA assembly (Golden Gate, Gibson...)",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    keywords="DNA assembly cloning simulator synthetic biology",
    scripts=["scripts/dnacauldron"],
    packages=find_packages(exclude="docs"),
    include_package_data=True,
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
    extras_require={"reports": ["pdf_reports"]},
)
