from altools.dna_sequence import DNASequence

class GlobalAlignment:
    def __init__(self, seq_a: DNASequence, seq_b: DNASequence):
        self.seq_a = seq_a
        self.seq_b = seq_b

    def __str__(self):
        return str(self.seq_a) + "\n" + str(self.seq_b) + "\n"