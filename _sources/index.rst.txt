DNA Cauldron Documentation
==========================

.. image:: _static/images/title.png
   :width: 500px
   :align: center


DNA Cauldron
===================

DNA Cauldron is a Python library to simulate restriction-based assembly operations. Provided the sequences of genetic parts and receptor vectors, DNA Cauldron will compute the assembli(es) that could result from the mix.

.. image:: _static/images/principle.png
   :width: 600px
   :align: center

DNA Cauldron was written for Synthetic Biology applications - typically, to predict and validate batches of parts-based assemblies. It is simple to use, plays well with BioPython, can import and export Genbank (it conserves all features), and provides advanced methods such as connector part auto-selection, backbone selection for linear parts, methods to select constructs subsets when dealing with large combinatorial assemblies.

**Try it online !** Use `this web service <http://cuba.genomefoundry.org/simulate_gg_assemblies>`_ to predict the outcome of a batch of (possibly combinatorial) Type 2S assemblies.

Installation
-------------

You can install DNA Cauldron through PIP


.. code:: shell

    sudo pip install dnacauldron

Alternatively, you can unzip the sources in a folder and type


.. code:: shell

    sudo python setup.py install


Usage
------

Single assembly
~~~~~~~~~~~~~~~

To assemble several parts and a receptor plasmid into a single construct,
use `single_assembly`. The parts can be provided either as paths to genbank
files or as Biopython records. DNA Cauldron returns a Biopython record of the
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

Full Assembly report
~~~~~~~~~~~~~~~~~~~~

DNA Cauldron also implements routine to generate reports on the assemblies,
featuring the resulting constructs (in genbank and PDF format) as well as
figures for verifying that the parts assembled as expected and help troubleshoot
if necessary.

The following code produces a structured directory with various reports:

.. code:: python

    import dnacauldron as dc
    parts = [
        dc.load_genbank("partA.gb", linear=False, name="PartA"),
        dc.load_genbank("partB.gb", linear=False, name="PartB"),
        dc.load_genbank("receptor.gb", linear=False, name="Receptor"),
    ]
    dc.full_assembly_report(parts, target="./my_report", enzyme="BsmBI",
                            max_assemblies=40, fragments_filters='auto',
                            assemblies_prefix='asm')

Result:

.. image:: _static/images/report_screenshot.jpg
   :width: 600px
   :align: center

How it works
------------

DNA Cauldron simulates enzyme digestions and computes sticky ends, then generates
a graph of the fragments that bind together, and explores circular paths in this graph
(which correspond to circular constructs), an idea also used in
`PyDNA <https://github.com/BjornFJohansson/pydna>`_ and first
described in `Pereira et al. Bioinf. 2015 <http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0544-x>`_ .
DNA Cauldron adds methods to deal with combinatorial assemblies, selecting constructs based on a marker, etc.

Contribute
----------

DNA Cauldron is an open-source library originally written at the
Edinburgh Genome Foundry by Zulko_ and is released on Github_ under the MIT
licence (¢ Edinburgh Genome Foundry), everyone is welcome to contribute.

.. raw:: html

       <a href="https://twitter.com/share" class="twitter-share-button"
       data-text="DNA Cauldron - A Python module for printing with living matter" data-size="large" data-hashtags="Bioprinting">Tweet
       </a>
       <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';
       if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';
       fjs.parentNode.insertBefore(js,fjs);}}(document, 'script', 'twitter-wjs');
       </script>
       <iframe src="http://ghbtns.com/github-btn.html?user=Edinburgh-Genome-Foundry&repo=DNA Cauldron&type=watch&count=true&size=large"
       allowtransparency="true" frameborder="0" scrolling="0" width="152px" height="30px" margin-bottom="30px"></iframe>


.. toctree::
    :hidden:
    :maxdepth: 3

    self

.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 3

    ref

.. toctree::
  :caption: Examples

  examples/single_assembly
  examples/full_report
  examples/backbone_autoselection
  examples/combinatorial_assembly
  examples/connector_selection



More biology software
----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Cauldron is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.

.. _Zulko: https://github.com/Zulko/
.. _Github: https://github.com/EdinburghGenomeFoundry/dnacauldron
.. _PYPI: https://pypi.python.org/pypi/dnacauldron
