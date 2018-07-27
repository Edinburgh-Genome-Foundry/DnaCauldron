.. raw:: html

    <p align="center">
    <img alt="DNA Cauldron Logo" title="DNA Cauldron Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaCauldron/master/docs/_static/images/title.png" width="500">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaCauldron.svg?branch=master
  :target: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaCauldron
  :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/DnaCauldron/badge.svg?branch=master
  :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/DnaCauldron?branch=master


DNA Cauldron (complete documentation `here <https://edinburgh-genome-foundry.github.io/DnaCauldron/>`_)
is a Python library to simulate restriction-based assembly operations.
Provided the sequences of genetic parts and receptor vectors, DNA Cauldron will compute the assembli(es) that could result from the mix.

DNA Cauldron was written for Synthetic Biology applications - typically, to predict and validate batches of parts-based assemblies. It is simple to use, plays well with BioPython, can import and export Genbank (it conserves all features), and provides advanced methods such as connector part auto-selection, backbone selection for linear parts, methods to select constructs subsets when dealing with large combinatorial assemblies.

**Try it online !** Use `this web service <http://cuba.genomefoundry.org/simulate_gg_assemblies>`_ to predict the outcome of a batch of (possibly combinatorial) Type 2S assemblies.

Installation
-------------

You can install DnaCauldron through PIP


.. code:: shell

    sudo pip install dnacauldron

Alternatively, you can unzip the sources in a folder and type


.. code:: shell

    sudo python setup.py install

It works better with the Networkx development version, that you install with

.. code:: shell

    sudo pip3 install git+https://github.com/networkx/networkx.git

Usage
------

Single assembly
~~~~~~~~~~~~~~~

To assemble several parts and a receptor plasmid into a single construct,
use `single_assembly`. The parts can be provided either as paths to genbank
files or as Biopython records. Dna Cauldron returns a Biopython record of the
final assembly, and (optionally) writes it to a Genbank file.

.. code:: python

    from dnacauldron.utils import single_assembly
    final_construct = single_assembly(
        parts=["partA.gb", "partB.gb", "partC.gb", "partD.gb", "receptor.gb"]
        outfile="final_construct.gb", # Name of the output
        enzyme="BsmBI" # enzyme used for the assembly
    )

Combinatorial assembly
~~~~~~~~~~~~~~~~~~~~~~

The following example imports parts from Genbank files and outputs all
possible outcomes of BsmBI-based Golden-Gate assembly as new genbank files
`001.gb`, `002.gb`, etc. We ignore the final assemblies containing a BsmBI site
as these are unstable.

.. code:: python

    from Bio import SeqIO # for exporting to Genbank
    from dnacauldron import (RestrictionLigationMix, NoRestrictionSiteFilter,
                             load_record)

    # Load all the parts (including the receptor)
    parts_files = ["partA.gb", "partA2.gb", "partB.gb", "partB2.gb",
                   "partC.gb", "receptor.gb"]
    parts = [load_record(filename, linear=False, id=filename)
             for filename in parts_files]

    # Create the "reaction mix"
    mix = RestrictionLigationMix(parts, enzyme='BsmBI')

    # Find all final assemblies (containing no sites from the restriction enzyme)
    assemblies = mix.compute_circular_assemblies()

    # Iter through all possible constructs and write them on disk as Genbanks.
    for i, assembly in enumerate(assemblies):
        SeqIO.write(assembly, os.path.join("..", "%03d.gb" % i), "genbank")


Full Assembly report
~~~~~~~~~~~~~~~~~~~~

DNA Cauldron also implements routine to generate reports on the assemblies,
featuring the resulting constructs (in genbank and PDF format) as well as
figures for verifying that the parts assembled as expected and help troubleshoot
if necessary.

The following code produces a structured directory with various reports:

.. code:: python

    from dnacauldron import load_genbank, full_assembly_report
    parts = [
        load_genbank("partA.gb", linear=False, name="PartA"),
        load_genbank("partB.gb", linear=False, name="PartB"),
        load_genbank("partC.gb", linear=False, name="PartC"),
        load_genbank("receptor.gb", linear=False, name="Receptor")
    ]
    dc.full_assembly_report(parts, target="./my_report", enzyme="BsmBI",
                            max_assemblies=40, fragments_filters='auto',
                            assemblies_prefix='asm')

Result:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaCauldron/master/docs/_static/images/report_screenshot.jpg
   :alt: [logo]
   :align: center
   :width: 600px


How it works
------------

Dna Cauldron simulates enzyme digestions and computes sticky ends, then generates
a graph of the fragments that bind together, and explores circular paths in this graph
(which correspond to circular constructs), an idea also used in
`PyDNA <https://github.com/BjornFJohansson/pydna>`_ and first
described in `Pereira et al. Bioinf. 2015 <http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0544-x>`_ .
DNA Cauldron adds methods to deal with combinatorial assemblies,
selecting constructs based on a marker, routines for report generation, etc.


Licence
--------

Dna Cauldron is an open-source software originally written at the `Edinburgh Genome Foundry
<http://www.genomefoundry.io>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaCauldron>`_ under the MIT licence (¢ Edinburgh Genome Foundry).
Everyone is welcome to contribute !

More biology software
----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Cauldron is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
