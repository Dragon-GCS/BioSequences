import re

from collections import Counter
from abc import ABC, abstractmethod

from bioseq import config
from bioseq.align import alignment
from bioseq.config import HYDROPATHY, MW, PK, NC_INFO, SYMBOL


class Sequence(ABC):
    def __init__(self, seq = ""):
        self._seq = seq.upper()

    @property
    def composition(self):
        """
        Analysis the composition of sequence
        Returns:
            composition: Dict of each element's appearance's times
        """
        if not hasattr(self, "_composition"):
            self._composition = {}
            counter = dict(Counter(self._seq))
            for key in sorted(counter):
                self._composition[key] = counter[key]
        return self._composition

    @property
    def length(self):
        return len(self._seq)

    @property
    def seq(self):
        """
        Seq can only be modified by mutation
        """
        return self._seq

    @property
    def weight(self):
        """
        Calculate the molecular Weight
        Returns:
            weight: Sequence's Mole Weight in Daltons
        """
        if not hasattr(self, "_weight"):
            weight_table = MW[self.__class__.__name__ + "_MW"]
            self._weight = round(sum([weight_table[e] for e in self._seq], 18), 2)
        return self._weight

    def align(self, subject, mode = 1, return_score = False):
        """
        Align two sequence
        Args:
            subject: Sequence to align
            mode: 1 - Use Needleman-Wunsch to global alignment
                  2 - Use Smith-Waterman to partial alignment
            return_score: whether return align score
        Returns:
            query: Self sequence after alignmnt
            subject: Subject sequence afer alignment
            score: align score if choose return_score
        """
        if isinstance(subject, self.__class__):
            subject = subject._seq
        elif not isinstance(subject, str):
            raise TypeError(f"Only str or {self.__class__.__name__} can be aligned to {self.__class__.__name__}")

        return alignment(self._seq, subject, mode, return_score)

    def find(self, target):
        """
        Find the target sequence in sequence and return the position
        Returns:
            : All position of target appearance in self.sequence
        """
        if isinstance(target, Sequence):
            target = target._seq

        elif not isinstance(target, (list, tuple)):
            raise TypeError("Target should be str or Sequence")

        return [i.start() for i in re.finditer(target, self._seq)]

    def mutation(self, position, target):
        """
        Modify the sequence
        Args:
            target: the sequence after mutation
            position: allow a position or a list of positions or the string to be modified
        Returns:
            seq: Sequence after modified
        """
        if isinstance(position, str):
            self._seq = self._seq.replace(position, target)
            return self._seq

        elif isinstance(position, int):
            position = [position]

        elif not isinstance(position, (list, tuple)):
            raise TypeError("Position should be str, int or list of int")

        if isinstance(target, Sequence):
            target = target._seq

        elif not isinstance(target, (list, tuple)):
            raise TypeError("Target should be str or Sequence")

        seq_list = list(self._seq)
        for pos in position:
            seq_list[pos: pos + len(target)] = target
        self._seq = "".join(seq_list)

        return self._seq

    @abstractmethod
    def _print(self):
        """
        Output sequence info
        """
        return self._seq

    def __add__(self, s):
        if isinstance(s, str):
            return self.__class__(self._seq + s)
        elif isinstance(s, self.__class__):
            return self.__class__(self._seq + s._seq)
        else:
            raise TypeError(f"Only str or {self.__class__.__name__} can be added to {self.__class__.__name__}")

    def __radd__(self, s):
        if isinstance(s, str):
            return self.__class__(s + self._seq)
        elif isinstance(s, self.__class__):
            return self.__class__(s._seq + self._seq)
        else:
            raise TypeError(f"Only str or {self.__class__.__name__} can add {self.__class__.__name__}")

    def __eq__(self, o):
        if isinstance(o, str):
            return self._seq == o
        elif isinstance(o, self.__class__):
            return self._seq == o._seq
        else:
            raise TypeError(f"Only str or {self.__class__.__name__} can compare to {self.__class__.__name__}")

    def __getitem__(self, index):
        return self._seq[index]

    def __len__(self):
        return len(self._seq)

    def __str__(self):
        return self._print()

    def __repr__(self):
        return self._print()


class RNA(Sequence):
    @property
    def reversed(self):
        """
        Return reversed sequence
        """
        return self.__class__(self._seq[::-1])

    @property
    def complemented(self):
        """
        Return complemented sequence
        """
        complement_seq = self.__class__()
        for bp in self._seq:
            complement_seq = NC_INFO[self.__class__.__name__ + "_COMPLEMENT"][bp] + complement_seq
        return complement_seq

    @property
    def GC(self):
        """
        Calculate the GC percentage
        """
        if not hasattr(self, "_GC"):
            self._GC = round((self.composition["C"] + self.composition["G"]) / self.length, 4)
        return self._GC

    def complement(self):
        """
        Change the sequence to complemented
        """
        self._seq = self.complemented._seq

    def reverse(self):
        """
        Change the sequence to reversed
        """
        self._seq = self.reversed._seq

    def getOrf(self, multi = False, replace = False):
        """
        Find the Open Reading Frame in sequence and save in self.orf
        Args:
            multi:      Return all Orf if true else only the longest
            replace:    Replace origin sequence with the longest Orf
        Returns:
            Orf:        Sequence's Orf got by search
        """
        self.orf = []
        # Traverse all Orf Frame
        for i in range(self.length - 5):
            if self._seq[i:i + 3] in config.START_CODON:
                # when start coden exits, from end start finding end coden
                end = i + (self.length - i) // 3 * 3
                for j in range(end, i + 5, -3):
                    if config.TABLE.get(self._seq[j - 3: j]) == "*":
                        orf = self._seq[i:j]
                        self.orf.append(orf)

                        if not multi:
                            if replace:
                                self._seq = orf
                            return self.orf  # Stop search if only request longest orf

        return self.orf

    def transcript(self, filter=True):
        """
        Transcript the sequence to peptide, the result will save in self.peptide
        Args:
            filter:  Return all product or the longest
        Returns:
            peptide: List of transcript product
        """
        orf = self.getOrf(multi = not filter) if not hasattr(self, "orf") else self.orf
        self.peptide = []
        for frame in orf:
            peptide = ""
            for i in range(0, len(frame), 3):
                peptide += config.TABLE.get(frame[i:i + 3], SYMBOL["printAlign"][1])
            self.peptide.append(Peptide(peptide[:-1]))

            if filter:
                return self.peptide

        return self.peptide

    def _print(self):
        if self._seq:
            return f"5'-{super()._print()}-3'"
        return ""


class DNA(RNA):
    def translate(self):
        if not hasattr(self, "_translate") or config.START_CODON != ["AUG"]:
            self._translate = RNA(self._seq.replace("T", "U"))
        return self._translate

    def transcript(self, filter = True):
        return self.translate().transcript(filter)



class Peptide(Sequence):
    @property
    def pI(self):
        """
        Calculate the peptide's pI which is a pH make peptide's charge equal zero
        """
        if not hasattr(self, "_pl"):
            min = 0.
            max = 14.
            pH = (min + max) / 2
            charge = self.chargeInpH(pH)
            while abs(charge) > 1e-8:
                if charge > 0:
                    min = pH
                else:
                    max = pH
                pH = (min + max) / 2
                charge = self.chargeInpH(pH)
            self._pl = round(pH, 3)
        return self._pl

    def chargeInpH(self, pH):
        """
        Calculate the charge amount of peptide at pH
        Henderson Hasselbalch equation: pH = pKa + log([A-]/[HA])
        Rearranging: [HA]/[A-] = 10 ** (pKa - pH)
        partial_charge =
            [A-]/[A]total = [A-]/([A-] + [HA]) = 1 / { ([A-] + [HA])/[A-] } =
            1 / (1 + [HA]/[A-]) = 1 / (1 + 10 ** (pKa - pH)) for acidic residues;
                                  1 / (1 + 10 ** (pH - pKa)) for basic residues
        """
        pos_charge = 1.0 / (10 ** (pH - PK["Nterm"]) + 1.0)
        for pos_aa, aa_pK in PK["pos_pK"].items():
            pos_charge += self.composition.get(pos_aa, 0) / (10 ** (pH - aa_pK) + 1.0)

        neg_charge = 1.0 / (10 ** (PK["Cterm"] - pH) + 1.0)
        for neg_aa, aa_pK in PK["neg_pK"].items():
            neg_charge += self.composition.get(neg_aa, 0) / (10 ** (aa_pK - pH) + 1.0)

        return pos_charge - neg_charge

    def getHphob(self, window_size = 9, show_img = False):
        """
        Calculate the Hydropathy Score.The lager the score, the higher the hydrophobicity
        Each aa's score is the average score of all aa in window_size.
        So part of Amino Acid at begin and end don't have score
        Args:
            window_size: the number for calculate averge hydropathy value
            show_img: whether to draw the result
        Returns:
            Hphob_list: the result of peptide's Hydropathy Score
        """
        if not hasattr(self, "_Hphob_list"):
            _Hphob = [HYDROPATHY[aa] for aa in self._seq]
            half_part = window_size // 2
            self._Hphob_list = []
            for i in range(half_part, len(_Hphob) - half_part):
                self._Hphob_list.append(round(sum(_Hphob[i - half_part: i + window_size - half_part]) / window_size, 3))

        if show_img:
            from matplotlib import pyplot as plt
            plt.title(f"Hydropathy Score for {self._seq[:4]}...{self._seq[-4:]}")
            plt.plot(range(window_size // 2 + 1, self.length - window_size // 2 + 1), self._Hphob_list, linewidth = 0.8)
            plt.xlim((1, self.length + 1))
            plt.grid(linestyle = "--")
            plt.xlabel("Position")
            plt.ylabel("Score")
            plt.show()

        return self._Hphob_list

    def _print(self):
        return f"N-{super()._print()}-C"
