Dna Cauldron
=============

Dna Cauldron is a Python library to simulate restriction-based assembly operations.
You provide a set of parts and receptor vectors and Dna Cauldron will compute the
assembli(es) that could result from the mix.

Dna Cauldron can import and export from Genbank, plays well with BioPython, and provides
ways to select or randomize certain assemblies when dealing with large combinatorial
assemblies.

Installation
-------------

You can install DnaCauldron through PIP
::
  sudo pip install dnacauldron

Alternatively, you can unzip the sources in a folder and type
::
  sudo python setup.py install


Usage
------

Single assembly
~~~~~~~~~~~~~~~

To assemble several parts and a receptor plasmid into a single construct,
use `single_assembly`. The parts can be provided either as paths to genbank
files or as Biopython records. Dna Cauldron returns a Biopython record of the
final assembly, and (optionally) writes it to a Genbank file.
::
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
::
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


Licence
--------

Dna Cauldron is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer>`_ under the MIT licence (¢ Edinburgh Genome Foundry).
Everyone is welcome to contribute !
