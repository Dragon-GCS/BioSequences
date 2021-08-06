from config import TABLE, START_CODON

class DNA():
    def __init__(self, sequence):
        self.sequence = sequence.upper() or ''
        self.length = len(sequence)
        self.orf = []
        self.peptide = []

    @property
    def complemented(self):
        """
        以5'-3'方向返回序列的互补链
        """
        bp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        return ''.join([bp[b] for b in self.sequence[::-1]])
    
    @property
    def reversed(self):
        """
        返回序列的3'-5'方向
        """
        return reversed(self.sequence)

    def get_orf(self):
        """
        Get open reading frame +1~+3 of the sequence
        Args:
            reverse: whether get reversed orf
        Returns:
            orf: the orf of the sequence
        """
        orf_len = int(self.length/3) * 3
        if not self.orf:
            for orf in range(3):

                for i in range(orf, orf_len, 3):

                    if self.sequence[i:i+3] == START_CODON:
                        seq = START_CODON

                        for j in range(i+3, orf_len, 3):
                            codon = self.sequence[j:j+3]
                            if len(codon) != 3:
                                break
                            seq += codon
                            if TABLE[codon] == "Stop":
                                self.orf.append(seq)
                                break
                            
            if not self.orf:
                self.orf = f'Sequence "{self.sequence:.8}..." has no ORF'
        return self.orf

    def translate(self, filter=True):
        """
        Translate the sequence to peptide(ORF +1 +2 +3)
        Args:
            filter: bool, return all or longest peptide 
        Returns:
            peptide: the translate product
        """
        if not self.orf:
            self.get_orf()

        if not self.peptide:
            if isinstance(self.orf, str):
                self.peptide = f'Sequence "{self.sequence:.8}..." can not translate because of no ORF '
            else:
                for orf in self.orf:
                    peptide = ''
                    for i in range(0, len(orf), 3):
                            aa = TABLE[orf[i:i+3]]
                            if aa == "Stop":
                                self.peptide.append(peptide)
                                break
                            peptide += TABLE[orf[i:i+3]]

                if filter:
                    self.peptide = sorted(self.peptide, key=len, reverse=True)[0]

        return self.peptide
    

if __name__ == "__main__":
    from config import test_dna
    d = DNA(test_dna)
    d.get_orf()
    d.translate()
    print(d.orf)
    print(d.peptide)
    c = DNA("ATCGTAGCAGS")
    c.translate()
    print(c.orf)
    print(c.peptide)
    n = DNA(d.complemented)
    n.translate(filter=False)
    print(n.orf)
    print(n.peptide)
