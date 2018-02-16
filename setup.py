import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('dnacauldron/version.py').read()) # loads __version__

setup(name='dnacauldron',
      version=__version__,
      author='Zulko',
    description='',
    long_description=open('pypi-readme.rst').read(),
    license='see LICENSE.txt',
    keywords="",
    packages= find_packages(exclude='docs'),
    install_requires=['Biopython', 'numpy', 'matplotlib', 'pandas', 'scipy',
                      'networkx', 'dna_features_viewer', 'flametree'])
