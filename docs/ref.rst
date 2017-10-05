.. reference ::

DnaCauldron Reference manual
==========================

Classes dependencies:
----------------------

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

.. automodule:: dnacauldron.utils
   :members:

Reporting
---------

.. automodule:: dnacauldron.reports.reports
  :members:


.. automodule:: dnacauldron.reports.plots
 :members:


AssemblyMix
-----------

.. automodule:: dnacauldron.AssemblyMix
   :members:


Filters
---------

.. automodule:: dnacauldron.Filter
   :members:


StickyEndsSeq
--------------

.. automodule:: dnacauldron.StickyEndsSeq
   :members:
