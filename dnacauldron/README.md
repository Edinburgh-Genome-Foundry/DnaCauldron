# Code organization

## Relationship between classes

The short story is that an ***AssemblyPlan*** contains ***Assemblies*** simulated using **Mixes** of DNA ***Fragments***. Now for more details:

- An ***AssemblyPlan*** is mostly a list of ***Assembly*** instances
- A ***SequenceRepository*** provides parts sequences to simulate an ***Assembly***.
- Simulating an ***Assembly*** results in an ***AssemblySimulation***, which contains:
  - Biopython *SeqRecords* of the predicted assembly constructs.
  - The list of ***AssemblyMix*** instances created during the simulation
  - The list of ***AssemblySimulationError*** detected during the simulation
- Simulating an ***AssemblyPlan*** results in an ***AssemblyPlanSimulation*** which is mostly a list of ***AssemblySimulation*** instances.
- An ***AssemblySimulation*** or ***AssemblyPlanSimulation*** can be written to a file using an ***AssemblyReportWriter***

The biology of the different cloning techniques is implemented in the subclasses of ***Assembly*** and ***AssemblyMix***.

### AssemblyMix subclasses

An ***AssemblyMix*** contains a list of ***Fragment*** instances (which subclass BioPython's *SeqRecord*), and can be:
- A ***StickyEndFragmentMix***:
  - Such mixes contain a list of ***StickyEndFragment*** instances.
  - A ***StickyEndFragment*** has a ***StickyEndSeq*** (subclass of BioPython's *Seq* with additional ***StickyEnd*** sequences on the left and right).
  - Usable Subclasses include ***RestrictionLigationMix***, ***Type2sRestrictionMix***, ***BASICAssemblyMix***.
- An ***HomologousFragmentMix***:
  - Such mixes contain a list of ***HomologousFragment*** instances.
  - Mixes also require a ***HomologyChecker*** to detect homologies.
  - Subclasses include ***LCRAssemblyMix***

### Assembly subclasses
