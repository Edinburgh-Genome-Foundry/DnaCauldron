"""Useful functions built on top of the DnaCauldron classes to simplify the
most common operations."""

import pandas
from ..AssemblyMix import RestrictionLigationMix, AssemblyError
from ..tools import reverse_complement, write_record
from .utils import autoselect_enzyme
import flametree

class BackboneChoice:
    """Class to represent the result of a backbone autoselection.

    Parameters
    ----------

    record
      Record that was analyzed, containing an insert, to be backbonized,
      and potentially an original backbone too.

    already_on_backbone
      Was a backbone detected in that record ?

    error

    backbone_record
      Record of the backone that was selected for this record among all
      provided choices.

    final_record
      Record featuring the original insert in the given ``record`` cloned into
      the auto-selected backbone.
    """

    def __init__(self, record, already_on_backbone=None, error=None,
                 backbone_record=None, final_record=None):
        """Initialize."""
        self.record = record
        self.already_on_backbone = already_on_backbone
        self.backbone_record = backbone_record
        self.final_record = final_record
        self.error = error

    def __repr__(self):
        """Write [recordname] (already/inserted on backbone)."""
        if self.already_on_backbone:
            return "%s (already on backbone)" % self.record.id
        elif (isinstance(self.backbone_record, str)):
            return "%s: no backbone found." % self.record.id
        else:
            return "%s inserted on %s" % (self.record.id,
                                          self.backbone_record.id)
    def to_dict(self):
        """Return the object as a dict. Used for dataframe/spreadsheet."""
        backbone, final_record = self.backbone_record, self.final_record
        detected = backbone.id if hasattr(backbone, 'id') else ''
        final_length = len(final_record) if final_record else  len(self.record)
        return dict(
            original_record=self.record.id,
            already_on_backbone='yes' if self.already_on_backbone else 'no',
            detected_backbone=detected,
            final_record_length=final_length,
            error=str(self.error.args[0]) if self.error else ''
        )

    @staticmethod
    def list_to_infos_spreadsheet(choices):
        """Return a pandas dataframe summarizing a list of BackboneChoices.

        The dataframe's columns are 'original_record', 'already_on_backbone',
        'detected_backbone', 'final_record_length', 'error'.
        """
        return pandas.DataFrame.from_records(
            [
                choice.to_dict()
                for choice in choices
            ],
            columns=['original_record', 'already_on_backbone',
                     'detected_backbone', 'final_record_length', 'error']
        )

    @staticmethod
    def write_final_records(choices, directory):
        """Write a list of BackboneChoices final records as genbanks."""
        target_dir = flametree.file_tree(directory)
        for choice in choices:
            if choice.already_on_backbone:
                record = choice.record
            else:
                record = choice.final_record
            if record is not None:
                record.name = choice.record.name
                record.id = choice.record.id
                target_file = target_dir._file('%s.gb' % choice.record.id)
                write_record(record, target_file, 'genbank')


def _get_insert_from_record(record, enzyme='BsmBI'):
    """Return the record of the one digested fragment without enzyme site."""
    mix = RestrictionLigationMix([record], enzyme=enzyme)
    inserts = [
        frag for frag in mix.filtered_fragments
        if not frag.is_reverse
    ]
    if len(inserts) != 1:
        raise ValueError("")
    return inserts[0]

def _standardize_overhangs_pair(overhangs):
    """Standardize a pair of overhangs (o1, o2).

    Returns either ``(o1, o2)`` or its reverse complement ``(rev_o2, rev_o1)``,
    whichever is smaller in alphabetical order.
    """
    o1, o2 = overhangs
    ro1, ro2 = [reverse_complement(o) for o in (o1, o2)]
    return min((o1, o2), (ro2, ro1))

def get_overhangs_from_record(record, enzyme='BsmBI', standardize=True):
    """Return a pair (o1, o2) of overhangs found by record digestion.

    If ``standardized=True`` the pair returned is either ``(o1, o2)`` or
    its reverse complement ``(rev_o2, rev_o1)``, which ever is smaller
    in alphabetical order.
    """
    insert = _get_insert_from_record(record, enzyme=enzyme)
    overhangs = str(insert.seq.left_end), str(insert.seq.right_end)
    return _standardize_overhangs_pair(overhangs) if standardize else overhangs

def _records_to_overhangs_dict(records, allow_multiple_choices=False):
    """Return ``{(o1, o2): rec}`` where o1, o2 are standardized overhangs.

    If ``allow_multiple_choices=True`` it will return ``{(o1, o2): [r1, r2]}``
    when several records share the same overhangs pair.

    If ``allow_multiple_choices=False`` and several records share the same
    overhangs pair, an error is raised.

    """
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
                result[overhangs] = [record]
            else:
                result[overhangs] = record
    return result

def record_contains_backbone(record, enzyme='BsmBI',
                             min_backbone_length=500):
    """Return True iff it believes the given record contains a backbone.

    A backbone is detected if, when cutting the circularized record with the
    given enzyme, there is one fragment with no site (the insert), and the rest
    has a total size above the given ``min_backbone_length``.
    """
    mix = RestrictionLigationMix([record], enzyme='BsmBI')
    fragments = [
        frag for frag in mix.filtered_fragments
        if not frag.is_reverse
    ]
    if fragments == []:
        raise AssemblyError('No site-less fragment found digesting record '
                            + record.id, record.id)
    insert = fragments[0]
    return (len(record) - len(insert)) > min_backbone_length

def swap_donor_vector_part(donor_vector, insert, enzyme,
                           allow_overhangs_edits=False):
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

def insert_parts_on_backbones(part_records, backbone_records,
                              enzyme='autodetect',
                              min_backbone_length=500,
                              process_parts_with_backbone=False,
                              default_backbone_choice=None):
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

    if enzyme == 'autodetect':
        enzyme = autoselect_enzyme(part_records, ['BsmBI', 'BsaI', 'BbsI'])

    overhangs_dict = _records_to_overhangs_dict(backbone_records)
    backbone_choices = []
    for record in part_records:
        try:
            has_backbone = record_contains_backbone(
                record, enzyme=enzyme, min_backbone_length=min_backbone_length)
            if (not process_parts_with_backbone) and has_backbone:
                choice = BackboneChoice(record, already_on_backbone=True)
            else:
                overhangs = get_overhangs_from_record(record, enzyme=enzyme)
                if overhangs in overhangs_dict:
                    backbone_record = overhangs_dict[overhangs]
                    final_record = swap_donor_vector_part(
                        donor_vector=backbone_record, insert=record,
                        enzyme=enzyme)
                    choice = BackboneChoice(record=record,
                                            already_on_backbone=has_backbone,
                                            backbone_record=backbone_record,
                                            final_record=final_record)
                else:
                    if default_backbone_choice is not None:
                        choice = default_backbone_choice(record)
                    else:
                        choice = BackboneChoice(
                            record=record,
                            already_on_backbone=has_backbone,
                            backbone_record='none found',
                            final_record=None
                        )
        except AssemblyError as e:
            choice = BackboneChoice(record, error=e)
        backbone_choices.append(choice)

    return backbone_choices
