.. raw:: html

    <p align="center">
    <img alt="DNA Cauldron Logo" title="DNA Cauldron Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaCauldron/master/docs/_static/images/title.png" width="500">
    <br /><br />
    </p>

.. image:: https://github.com/Edinburgh-Genome-Foundry/DnaCauldron/actions/workflows/build.yml/badge.svg
    :target: https://github.com/Edinburgh-Genome-Foundry/DnaCauldron/actions/workflows/build.yml
    :alt: GitHub CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/DnaCauldron/badge.svg?branch=master
    :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/DnaCauldron?branch=master


DNA Cauldron
============

DNA Cauldron (full documentation `here <https://edinburgh-genome-foundry.github.io/DnaCauldron/>`_)
is a generic cloning simulation framework to predict
final construct sequences and detect assembly flaws. It aims in particular at
automating the simulation and verification of large (and possibly multi-step)
DNA assembly batches, with the idea that, if your assembly plan works with
Cauldron, you'll get the same results at the bench.

- Great support for Golden Gate assembly (incl. MoClo, EMMA, Phytobrick, etc.), the simulation is as simple as
  `dropping sequence files in a web app <http://cuba.genomefoundry.org/simulate_gg_assemblies>`_.
- Good support for Gibson, Biobrick, BASIC, LCR Assembly.
- Provide genetic parts in any order, with or without annotations, in reverse or direct
  sense, linear or circular, single or combinatorial assembly... Cauldron will get it!
- Design flaw detection (missing parts, unwanted restriction sites, wrong overhang designs, etc.).
- If your assembly needs connector parts, Cauldron can select them for you!
- Import batch assembly plans from spreadsheets, including hierarchical (=multi-step) plans.
- Comprehensive reports with constructs sequences, fragment sequences, summary spreadsheets.
- Export multi-file reports in folders, zip files, or in-memory zip files (for use on servers).

.. raw:: html

    <p align="center">
    <img alt="DNA Cauldron Logo" title="DNA Cauldron Logo" src="https://github.com/Edinburgh-Genome-Foundry/DnaCauldron/raw/master/docs/_static/images/reports_elements.png" width="800">
    <br /><br />
    </p>


Usage tutorial
--------------

Providing part sequences
~~~~~~~~~~~~~~~~~~~~~~~~

To simulate an assembly, you need first to provide part sequences. In DNA Cauldron, sequences
are managed by a ``SequenceRepository``, created as follows:

.. code:: python

    import dnacauldron as dc
    
    # Create a repository from BioPython records
    repository = dc.SequenceRepository(collections={"parts": {"part_A": record_A, "part_B": record_B}})
    
    # Or use "import_sequences" to import files from folders, zip files, etc.
    repository = dc.SequenceRepository()
    repository.import_sequences(folder="my_sequences/", use_file_names_as_ids=True)

**Important:** ensure that the part's topology (as set in Biopython records at
``record.annotations['topology']``) is accurate. This can be done at import
time by setting ``topology='linear'``, ``topology='circular'``, or for instance
``topology='default_to_linear'`` to use the topology specified by each Genbank and
default to linear if none is specified.


Parts assembly
~~~~~~~~~~~~~~

An assembly simply defines which parts of the sequence repository should
be assembled together. Here we use a "Type-2s Restriction" class of assembly,
which will work for all Golden-Gate-style assemblies:

.. code:: python

    assembly = dc.Type2sRestrictionAssembly(parts=["part_A", "part_B", "receptor"])
    simulation = assembly.simulate(sequence_repository=repository)

Note that we could have provided the enzyme in ``Type2sRestrictionAssembly`` with
``enzyme='BsmBI'``, but it will be auto-selected by Cauldron.

If you want to simulate other restriction-based assembly reactions such as IGEM Biobricks
or Gibson Assembly instead of Type2S restriction, use the corresponding built-in Assembly subclass:

.. code:: python

    assembly = dc.BioBrickStandardAssembly(parts=['part_A', 'part_B'])
    assembly = dc.GibsonAssembly(parts=['part_A', 'part_B'])
    

Now you can explore the results of the simulation:

.. code:: python

    # Print the ID and length of the generated construct(s)
    for record in simulation.construct_records:
        print (record.id, len(record))
    
    # Get a list of dictionnaries with data on each construct
    constructs_data = simulation.compute_all_construct_data_dicts()
    
    # Write a full report with sequences and figures in a zip.
    simulation.write_report("report.zip")

DNA Cauldron aims at generating reports as useful as possible to help you
pinpoint any problem when you don't get the expected number of assemblies.


Assembly Plans
~~~~~~~~~~~~~~

An assembly plan is simply defined by a list of assemblies:

.. code:: python

   # Define an assembly plan as a list of Assembly objects
   assembly_plan = dc.AssemblyPlan(assemblies=[assembly_1, ...])
   
   # Or import an assembly plan from spreadsheets:
   assembly_plan = dc.AssemblyPlan.from_spreadsheet(
       spreadsheet="batch_1.csv", # could also be an xls(x) file
       assembly_class=dc.Type2sRestrictionAssembly
   )

See these different examples for a spreadsheet defining assemblies.
   
Assembly plans can be hierarchical (i.e. have an assembly's construct serve as a
part in another assembly). DNA Cauldron will automatically figure out the dependencies
between assemblies and sort the order in which they should be simulated.

The simulation and reporting on an assembly plan is very similar to that of a single assembly:

.. code:: python

   plan_simulation = assembly_plan.simulate(sequence_respository)
   
   # Get a list of dictionnaries with data on each construct
   plan_simulation.compute_all_construct_data_dicts()
   
   # Write a detailed report on each assembly and on the plan as a whole
   plan_simulation.write_report("my_assembly_simulation.zip")


Installation
------------

You can install DNA Cauldron through PIP:

.. code:: shell

    pip install dnacauldron

The full installation using ``dnacauldron[reports]`` is required for report generation.
Alternatively, you can unzip the sources in a folder and type:

.. code:: shell

    python setup.py install


How it works
------------

DNA Cauldron predicts circular constructs by finding circular paths in part
homology graphs, an idea first described in
`Pereira et al. (2015) <http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0544-x>`_,
and used in the Python library `PyDNA <https://github.com/BjornFJohansson/pydna>`_.


License = MIT
-------------

DNA Cauldron is an open-source software originally written at the `Edinburgh Genome Foundry
<http://www.genomefoundry.io>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaCauldron>`_ under the MIT license (Copyright 2020 Edinburgh Genome Foundry).
Everyone is welcome to contribute!


More biology software
---------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Cauldron is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
