Dna Cauldron
=============

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaCauldron/master/docs/logo.png
   :alt: [logo]
   :align: center
   :width: 600px


Dna Cauldron is a Python library to simulate restriction-based assembly operations.
You provide a set of parts and receptor vectors and Dna Cauldron will compute the
assembli(es) that could result from the mix.

Dna Cauldron can import and export from Genbank, plays well with BioPython, and provides
ways to select or randomize certain assemblies when dealing with large combinatorial
assemblies.

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
        parts_filenames=["partA.gb", "partB.gb", "partC.gb", "partD.gb"],
        receptor="receptor.gb", # Receptor plasmid for the final assembly
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
                             load_genbank)
    enzyme = "BsmBI"
    filters = [NoRestrictionSiteFilter(enzyme)]
    parts_files = ["partA.gb", "partA2.gb", "partB.gb", "partB2.gb", "partC.gb",
                "receptor.gb"]
    parts = [load_genbank(filename, linear=False) for filename in parts_files]
    mix = RestrictionLigationMix(parts, enzyme)
    assemblies = mix.compute_circular_assemblies(seqrecord_filters=filters)
    for i, assembly in enumerate(assemblies):
        SeqIO.write(assembly, os.path.join("..", "%03d.gb" % i), "genbank")

How it works
------------

Dna Cauldron simulates enzyme digestions and computes sticky ends, then generates
a graph of the fragments that bind together, and explores circular paths in this graph
(which correspond to circular constructs), an idea also used in
`PyDNA <https://github.com/BjornFJohansson/pydna>`_ and first
described in `Pereira et al. Bioinf. 2015 <http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0544-x>`_ .
DNA Cauldron adds methods to deal with combinatorial assemblies, selecting constructs based on a marker, etc.


Licence
--------

Dna Cauldron is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer>`_ under the MIT licence (¢ Edinburgh Genome Foundry).
Everyone is welcome to contribute !
