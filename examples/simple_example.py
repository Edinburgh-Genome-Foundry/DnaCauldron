import os

from Bio import SeqIO
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicRecord

from dnacauldron import (RestrictionLigationMix,
                         NoRestrictionSiteFilter,
                         TextSearchFilter,
                         load_genbank)


parts = [load_genbank(os.path.join("data", "%s.gb" % name), linear=False)
         for name in ["partA", "partB", "partC", "receptor"]]
enzyme = "BsmBI"
mix = RestrictionLigationMix(parts, enzyme)
assemblies = list(mix.compute_circular_assemblies(
    seqrecord_filters=[
        NoRestrictionSiteFilter(enzyme),
    ])
)

SeqIO.write(assemblies[0], "final_sequence.gb", "genbank")


def plot_constructs(constructs, title, filename):
        fig, ax = plt.subplots(1, len(constructs),
                               figsize=(max(7,2*len(constructs)),2))
        if len(constructs) == 1:
            ax = [ax]
        ax[0].set_title(title)
        for i, construct in enumerate(constructs):
            GraphicRecord.from_biopython_record(
            construct,
            fun_label=lambda f: f.qualifiers["label"][0],
            fun_color=lambda f: {"cds": "red",
                                 "rep._origin": "green"
                                }.get(f.type, "blue"),
        ).plot(ax=ax[i], with_ruler=False)
        fig.tight_layout()
        fig.savefig(filename, bbox_inches="tight")



plot_constructs(mix.constructs, "ORIGINAL CONSTRUCTS", "original.png")
plot_constructs(mix.fragments, "FRAGMENTS", "fragments.png")
plot_constructs(assemblies, "ASSEMBLIES", "assemblies.png")
