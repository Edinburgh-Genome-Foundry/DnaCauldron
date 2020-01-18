from Bio.SeqRecord import SeqRecord
from dna_features_viewer import BiopythonTranslator

class Fragment(SeqRecord):
    
    def plot(self, ax=None):
        graphic_record = BiopythonTranslator().translate_record(self)
        ax, _ = graphic_record.plot(ax=ax, strand_in_label_threshold=7)
        return ax
