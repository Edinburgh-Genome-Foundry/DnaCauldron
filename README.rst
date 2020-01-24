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


DNA Cauldron
============

Don't you hate it when you are planning a large DNA assembly batch and you're not sure how the many parts will assemble together, so when some assemblies fail at the bench you spend hours checking your files for design flaws?

DNA Cauldron might be for you! The library provides a generic cloning simulation framework to predict constructs sequences and detect assembly design flaws for various cloning methods:

- Advanced design flaw detection (missing parts, presence of problematic
restriction sites in the parts, wrong overhang designs, etc.)
- Awesome support for Type 2S assembly (automatic enzyme and connector part selection)
- Support for other restriction-based assemblies (BioBrick, BASIC) or homology based assemblies (Gibson Assembly, Ligase Cycling Reaction Assembly).

The library also enables the management of large, complex assembly batches:

- Import assembly plans from spreadsheets
- Simulate and validate hierarchical assembly plans where some assembly products are used as parts in the next level of assembly.
- Generate complete, self contained reports, in folders, zip files, in memory or in JSON format (for use on servers).

Usage
-----

Providing part sequences
~~~~~~~~~~~~~~~~~~~~~~~~

To simulate an assembly, you need first to provide parts sequences. In DNA Cauldron, these
are managed via a ``SequenceRepository``:

.. code:: python

    import dnacauldron as dc
    
    # Create a repository from BioPython records
    repository = dc.SequenceRepository(parts={"part_A": record_1, "part_B":...})
    
    # Or use "import_sequences" to import files from folders, zip files, etc.
    repository = dc.SequenceRepository()
    repository.import_sequences(folder="my_sequences/", use_file_names_as_ids=True)

One important thing is to make sure the parts topology (as set in Biopython records at ``record.annotations['topology']``) is accurate. This can be done at import time by setting ``topology='linear'``, ``topology='circular'``, or ``topology='auto'`` to use the topology specified by each Genbank file.  

Parts assembly
~~~~~~~~~~~~~~

You define and simulate an assembly by explaining which parts of the repository should
be assembled together:

.. code:: python

    assembly = dc.Type2sRestrictionAssembly(parts=["part_A", "part_B", "receptor"])
    simulation = assembly.simulate(sequence_repository=repository)

Note that we could have provided ``enzyme='BsmBI'`` in ``Type2sRestrictionAssembly``
but the enzyme will be auto-selected by the framework.

If you want to simulate other restriction-based assembly reactions such as IGEM Biobricks
or BASIC assembly instead of Type2s restriction, use the corresponding built-in Assembly subclass:

.. code:: python

    assembly = dc.BioBrickStandardAssembly(parts=['part_A', 'part_B'])
    assembly = dc.BASICAssembly(parts=['part_A', 'part_B'])
    

Now you can 

.. code:: python

    # Print the ID and length of the generated construct(s)
    for record in simulation.construct_records:
        print (record.id, len(record))
    
    # Get a list of dictionnaries with data on each construct
    constructs_data = simulation.get_all_constructs_data()
    
    # Write a full report with sequences and figures
    simulation.write_report("report/")

DNA Cauldron shines in the informative reports it produces to help you pinpoint problems
These reports can be configured and customized to show more or less information.
[Report screenshot]

Assembly Plans
~~~~~~~~~~~~~~

An assembly plan is simply defined by a list of assemblies:

.. code:: python

   # Define an assembly plan as a list of Assembly objects
   assembly_plan = dc.AssemblyPlan(assemblies=[...])
   
   # Or import an assembly plan from spreadsheets:
   assembly_plan = dc.AssemblyPlan.from_spreadsheet(
       spreadsheet="batch_1.csv", # could also be an xls(x) file
       assembly_class=dc.Type2sRestrictionAssembly
   )
See this example for a spreadsheet defining assemblies (or this example for a spreadsheet with more parameters).
   
A nice thing is that assembly plans can be hierarchical (i.e. have an assembly's construct serve as a part
of another assembly). DNA Cauldron will automatically figure out the dependencies
between assemblies and sort the order in which they should be simulated.

The simulation and reporting on an assembly plan is very similar to that of a single assembly:

.. code:: python

   plan_simulation = assembly_plan.simulate(sequence_repository=sequence_respository)
   
   # Get a list of dictionnaries with data on each construct
   plan_simulation.get_all_constructs_data()
   
   plan_simulation.write_report("my_assembly_simulation.zip")
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

    from dnacauldron import (Type2sRestrictionMix, NoRestrictionSiteFilter,
                             load_record, write_record)

    # Load all the parts (including the receptor)
    parts_files = ["partA.gb", "partA2.gb", "partB.gb", "partB2.gb",
                   "partC.gb", "receptor.gb"]
    parts = [load_record(filename, topology='circular', id=filename)
             for filename in parts_files]

    # Create the "reaction mix"
    mix = Type2sRestrictionMix(parts, enzyme='BsmBI')

    # Find all final assemblies (containing no sites from the restriction enzyme)
    assemblies = mix.compute_circular_assemblies()

    # Iter through all possible constructs and write them on disk as Genbanks.
    for i, assembly in enumerate(assemblies):
        out_path = os.path.join("..", "%03d.gb" % i)
        write_record(assembly, out_path, "genbank")


Full Assembly report
~~~~~~~~~~~~~~~~~~~~

DNA Cauldron also implements routine to generate reports on the assemblies,
featuring the resulting constructs (in genbank and PDF format) as well as
figures for verifying that the parts assembled as expected and help troubleshoot
if necessary.

The following code produces a structured directory with various reports:

.. code:: python

    from dnacauldron import load_record, full_assembly_report
    parts = [
        load_record("partA.gb", topology='circular', name="PartA"),
        load_record("partB.gb", topology='circular', name="PartB"),
        load_record("partC.gb", topology='circular', name="PartC"),
        load_record("receptor.gb", topology='circular', name="Receptor")
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
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaCauldron>`_ under the MIT licence (Copyright Edinburgh Genome Foundry).
Everyone is welcome to contribute !

More biology software
----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Cauldron is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
