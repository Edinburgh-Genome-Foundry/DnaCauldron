.. _reference:

DNA Cauldron Reference manual
=============================

**Classes dependencies:**

.. mermaid::
   :align: center

    graph TD;
        StickyEnd["StickyEnd<br/>(DNA fragment protusion)"]
        StickyEndSeq["StickyEndSeq<br/>(DNA fragment with sticky ends)"]
        StickyEndRecord["StickyEndRecord<br/>(BioPython record with StickyEndSeq sequence)"]
        AssemblyMix
        RestrictionLigationMix
        StickyEnd-->|can be found at the end of...| StickyEndSeq;
        StickyEndSeq-->StickyEndRecord;
        StickyEndRecord -->|several can be grouped into...| FragmentsChain
        AssemblyMix  --> |is base class of...| RestrictionLigationMix
        RestrictionLigationMix -->|computes digestions resulting in...| StickyEndRecord
        FragmentsChain --> |fragments are assembled together into a...| Biopython-SeqRecord
        style Biopython-SeqRecord fill:#fff;


High-level functions
--------------------

.. autofunction:: dnacauldron.utils.insert_parts_on_backbones
.. autofunction:: dnacauldron.utils.single_assembly
.. autofunction:: dnacauldron.utils.record_contains_backbone
.. autofunction:: dnacauldron.utils.autoselect_enzyme

Reporting
---------

.. automodule:: dnacauldron.reports.reports
  :members:


.. automodule:: dnacauldron.reports.plots
 :members:

AssemblyMix
-----------

.. automodule:: dnacauldron.AssemblyMix.AssemblyMix
  :members:

.. automodule:: dnacauldron.AssemblyMix.BASICLigationMix
  :members:


Lower-level classes
-------------------

Filters
~~~~~~~

.. automodule:: dnacauldron.AssemblyMix.Filter
   :members:

StickyEndsSeq
~~~~~~~~~~~~~

.. autoclass:: dnacauldron.StickyEndsSeq.StickyEnd
   :members:
.. autoclass:: dnacauldron.StickyEndsSeq.StickyEndsSeq
   :members:
.. autoclass:: dnacauldron.StickyEndsSeq.StickyEndsSeqRecord
   :members:

BackBoneChoice
~~~~~~~~~~~~~~

.. autoclass:: dnacauldron.utils.insert_parts_on_backbones.BackboneChoice
  :members:
