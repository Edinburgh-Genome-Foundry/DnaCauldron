Dna Cauldron
=============


Usage
-----

The simplest way to assemble Genbank files is with ``:

from dnachisel.utils import single_assembly
genbank_files_assembly(
    parts_filenames=["partA.gb", "partB.gb", "partC.gb", "partD.gb"],
    receptor="receptor.gb", # Receptor plasmid for the final assembly
    outfile="final_construct.gb", # Name of the output
    enzyme="BsmBI" # enzyme used for the assembly
)


Installation
--------------

dnacauldron can be installed by unzipping the source code in one directory and using this command: ::

    sudo python setup.py install

You can also install it directly from the Python Package Index with this command: ::

    sudo pip dnacauldron install


Licence
--------

Dna Cauldron is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer>`_ under the MIT licence (¢ Edinburgh Genome Foundry).
Everyone is welcome to contribute !
