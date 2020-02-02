import os
import dnacauldron as dc

this_directory = os.path.dirname(os.path.realpath(__file__))
parts_fasta= os.path.join(this_directory, "RFP_GFP_plasmid_parts.fa")
oligos_fasta = os.path.join(this_directory, "./RFP_GFP_plasmid_BOs.fa")
assembly_plan_path = os.path.join(this_directory, "assembly_plan.csv")

def test_lcr_assembly():
    repo = dc.SequenceRepository()
    repo.import_records(files=[oligos_fasta, parts_fasta])
    plan = dc.AssemblyPlan.from_spreadsheet(path=assembly_plan_path)
    simulation = plan.simulate(repo)
    stats = simulation.compute_stats()
    assert stats["valid_assemblies"] == 1
    # The second assembly is flawed on purpose
    assert stats["errored_assemblies"] == 1

    # Coverage!
    simulation.write_report("@memory")