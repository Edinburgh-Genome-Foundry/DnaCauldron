import dnacauldron as dc
import os

parts_dir = os.path.join("data", "assemblies")
records_dict = {
    filename.split('.')[0]: dc.load_record(
        os.path.join(parts_dir, filename),
        topology="circular",
        id=filename.split('.')[0],
    )
    for filename in os.listdir(parts_dir)
}

constructs = {
    "C1": ["receptor", "partA", "partB", "partC"],
    "C2": ["receptor", "partA2", "partB2", "partC"],
    "C3": ["receptor", "partA", "partA2", "partB", "partC"],
}
part_names = set([p for parts in constructs.values() for p in parts])
parts = {name: records_dict[name] for name in part_names}
plan = [
    (construct, [parts[p] for p in construct_parts])
    for construct, construct_parts in constructs.items()
]
dc.full_assembly_plan_report(
    plan,
    target=os.path.join("output_data", "report.zip"),
    enzyme="autoselect",
    assert_single_assemblies=False,
    fail_silently=True,
)
