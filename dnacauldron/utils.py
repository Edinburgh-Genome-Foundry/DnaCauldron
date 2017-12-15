"""Useful functions to simplify the most common operations."""

from Bio import SeqIO, Restriction
import pandas

from .Filter import NoRestrictionSiteFilter
from .AssemblyMix import RestrictionLigationMix, AssemblyError
from .tools import reverse_complement

def autoselect_enzyme(parts, enzymes):
    """Finds the enzyme that the parts were probably meant to be assembled with

    Parameters
    ----------

    parts
      A list of SeqRecord files. They should have a "linear" attribute set to
      True or False, otherwise

    Returns
    --------
    The enzyme that has as near as possible as exactly 2 sites in the different
    constructs.
    """
    def enzyme_fit_score(enzyme_name):
        enz = Restriction.__dict__[enzyme_name]
        return sum([
            abs(2-len(enz.search(part.seq,
                                 linear=part.__dict__.get('linear', False))))
            for part in parts
        ])
    return min(enzymes, key=enzyme_fit_score)

def single_assembly(parts, receptor, outfile=None,
                    enzyme="BsmBI", annotate_homologies=True,
                    mix_class="restriction"):
    """Return the single assembly obtained by assembling together different
    parts on a receptor vector.

    Parameters
    ----------

    parts
      A list of either filenames or Biopython records or parts. They are
      assumed circular (i.e. parts on backbones) but linear parts should
      work too as long as the receptor is circular. To make sure that
      a part will be treated as linear DNA, provide a Biopython record
      for that part, with a ``.linear`` attribute set to true.

    outfile
      Name of a genbank file where to output the result

    enzyme
      Name of the enzyme used for the assembly.
    """

    if mix_class == "restriction":
        mix_class = RestrictionLigationMix

    def load_genbank(filename):
        """Specific loag_genbank flavor"""
        record = SeqIO.read(filename, "genbank")
        record.linear = False
        record.name = filename.split("/")[-1].split(".")[0].lower()
        if isinstance(receptor, str) and (filename == receptor):
            record.name += " (RECEPTOR)"
        return record
    part_records = [
        load_genbank(part) if isinstance(part, str) else part
        for part in parts + [receptor]
    ]
    biopython_enzyme = Restriction.__dict__[enzyme]
    sites_in_receptor = \
        len(biopython_enzyme.search(part_records[-1].seq, linear=False))

    def exactly_one_receptor_vector(fragments):
        receptor_fragments = [
            fragment for fragment in fragments
            if "(RECEPTOR)" in fragment.original_construct.name
        ]
        return len(receptor_fragments) == sites_in_receptor - 1

    mix = mix_class(part_records, enzyme)
    assemblies = mix.compute_circular_assemblies(
        fragments_sets_filters=[exactly_one_receptor_vector],
        annotate_homologies=annotate_homologies
    )
    assemblies = list(assemblies)
    N = len(assemblies)
    if N != 1:
        raise AssemblyError('Found %d assemblies instead of 1 expected' % N)
    assembly = assemblies[0]
    if outfile is not None:
        SeqIO.write(assembly, outfile, "genbank")
    return assembly

def swap_donor_vector_part(donor_vector, insert, enzyme):
    """Return the records obtained by cloning inserts into a donor vector.

    Meant for Type-2S assembly standards only (Golden Gate, etc.)

    This method is meant to quickly go from a linearized sequence of a part
    to a circular vector (the part in its donor vector) by starting from
    an existing donor vector (with same overhangs) and swapping this vector's
    part for the insert of interest.

    Parameters
    ----------
    donor_vector
      Biopython record of a donor vector. must have an insert producing a
      restriction-free fragment

    insert
      Biopython record of a plasmid or a linear DNA sequence containing an
      insert (i.e. a fragment that is cut out)

    enzyme
      The name of the enzyme to use e.g. 'BsmBI', 'BsaI', ...

    """

    mix = RestrictionLigationMix([donor_vector], enzyme=enzyme)
    donor_fragments = [
        f for f in mix.fragments
        if len(mix.enzyme.search('A' + f.seq.to_standard_sequence())) > 0
    ]
    for fr in donor_fragments:
        fr.features = [f for f in fr.features
                       if 'source' not in f.qualifiers]
    assert(len(donor_fragments) == len(mix.fragments) - 1)

    mix = RestrictionLigationMix([insert], enzyme=enzyme)
    insert_fragments = [
        f for f in mix.fragments
        if len(mix.enzyme.search('A' + f.seq.to_standard_sequence())) == 0
    ]
    assert(len(insert_fragments) == 1)
    mix = RestrictionLigationMix(
        fragments=[insert_fragments[0]] + donor_fragments,
        enzyme=enzyme,
        fragments_filters=()
    )
    assemblies = list(mix.compute_circular_assemblies())
    assert (len(assemblies) == 1)
    return assemblies[0]


class BackboneChoice:
    """Class to represent the result of a backbone autoselection"""

    def __init__(self, record, already_on_backbone,
                 backbone_record=None, final_record=None):
        self.record = record
        self.already_on_backbone = already_on_backbone
        self.backbone_record = backbone_record
        self.final_record = final_record

    def __repr__(self):
        if self.already_on_backbone:
            return "%s (already on backbone)" % self.record.id
        else:
            return "%s inserted on %s" % (self.record.id,
                                          self.backbone_record.id)
    def to_dict(self):
        return dict(original_record=self.record.id,
                    already_on_backbone=self.already_on_backbone,
                    detected_backbone=self.backbone_record.id if
                                      hasattr(self.backbone_record, 'id')
                                      else '',
                    final_record_length=len(self.final_record) if
                                        self.final_record else
                                        len(self.record))

    @staticmethod
    def list_to_infos_spreadsheet(choices):
        return pandas.DataFrame.from_records(
            [
                choice.to_dict()
                for choice in choices
            ],
            columns=['original_record', 'already_on_backbone',
                     'detected_backbone', 'final_record_length']
        )

    @staticmethod
    def write_final_records(choices, directory):
        for choice in choices:
            pass


def insert_parts_on_backbones(part_records, backbone_records,
                              enzyme='autodetect',
                              min_backbone_length=500,
                              process_parts_with_backbone=False):
    """Autodetect the right backbone for each Golden Gate part.

    This method is meant to process a batch of genbank files, some of
    which might represent a part on a backbone, and some of which
    represent simply a part (and enzyme-flanked overhangs) which needs
    to be complemented with the right backbone.

    It will return, for each part, whether it has already a backbone, and
    if not, which backbone was selected and what the final sequence is.

    Parameters
    ----------

    part_records
      List of genbanks of the parts to put on vectors.

    backbone_vectors
      Vectors to insert parts in, typically donor vectors for different
      positions of an assembly standard.

    enzyme
      Enzyme to use. Use autodetect for autodetection.

    min_backbone_length
      Minimal length of a backbone. Used to determine if a part is
      represented alone or with a backbone.

    process_parts_with_backbone
      If true, parts will be inserted in an autoselected backbone even
      when they already have a backbone (it will be replaced).

    """

    # HELPER FUNCTIONS

    def record_contains_backbone(record, enzyme='BsmBI',
                                 min_backbone_length=500):
        mix = RestrictionLigationMix([record], enzyme='BsmBI')
        insert = [
            frag for frag in mix.filtered_fragments
            if not frag.is_reverse
        ][0]
        return (len(record) - len(insert)) > min_backbone_length

    def get_insert_from_record(record, enzyme='BsmBI'):
        mix = RestrictionLigationMix([record], enzyme=enzyme)
        return [
            frag for frag in mix.filtered_fragments
            if not frag.is_reverse
        ][0]

    def standardize_overhangs(overhangs):
        o1, o2 = overhangs
        ro1, ro2 = [reverse_complement(o) for o in (o1, o2)]
        return min((o1, o2), (ro2, ro1))

    def get_overhangs_from_record(record, enzyme='BsmBI'):
        insert = get_insert_from_record(record, enzyme=enzyme)
        overhangs = str(insert.seq.left_end), str(insert.seq.right_end)
        return standardize_overhangs(overhangs)

    def records_to_overhangs_dict(records, allow_multiple_choices=False):
        result = {}
        for record in records:
            overhangs = get_overhangs_from_record(record)
            if overhangs in result:
                if allow_multiple_choices:
                    result[overhangs].append(record)
                else:
                    raise ValueError("Vector %s has same overhangs as %s"
                                     % (record.id, result[overhangs].id))
            else:
                if allow_multiple_choices:
                    result[overhangs] = []
                else:
                    result[overhangs] = record
        return result

    # MAIN SCRIPT

    if enzyme == 'autodetect':
        enzyme = autoselect_enzyme(part_records, ['BsmBI', 'BsaI', 'BbsI'])

    overhangs_dict = records_to_overhangs_dict(backbone_records)
    backbone_choices = []
    for record in part_records:
        if (not process_parts_with_backbone) and record_contains_backbone(
            record, enzyme=enzyme, min_backbone_length=min_backbone_length):
            choice = BackboneChoice(record, already_on_backbone=True)
        else:
            overhangs = get_overhangs_from_record(record, enzyme=enzyme)
            if overhangs in overhangs_dict:
                backbone_record = overhangs_dict[overhangs]
                final_record = swap_donor_vector_part(
                    donor_vector=backbone_record, insert=record, enzyme=enzyme)
                final_record.id = record.id
                choice = BackboneChoice(record=record,
                                        already_on_backbone=False,
                                        backbone_record=backbone_record,
                                        final_record=final_record)
            else:
                choice = BackboneChoice(record=record,
                                        already_on_backbone=False,
                                        backbone_record='none found',
                                        final_record=None)
        backbone_choices.append(choice)

    return backbone_choices
