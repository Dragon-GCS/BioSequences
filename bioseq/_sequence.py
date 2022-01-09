import re

from collections import Counter
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, TypeVar, Union, Iterable

from bioseq import config, algorithm
from bioseq.config import AlignmentConfig


class Sequence:
    info: str
    _seq: str
    _composition: Dict[str, int]
    _weight: float

    def __init__(self, seq: str = "", info:str = ""):
        self._seq = seq.upper()
        self.info = info
        self.reset_cache()

    def reset_cache(self):
        """
        Reset some cache property
        """
        self._composition = {}
        self._weight = 0.

    @property
    def composition(self) -> Dict[str, int]:
        """
        Analysis the composition of sequence
        Returns:
            composition: Dict of each element's appearance's times
        """
        if not self._composition:
            counter = dict(Counter(self._seq))
            for key in sorted(counter):
                self._composition[key] = counter[key]
        return self._composition

    @property
    def length(self) -> int:
        return len(self._seq)

    @property
    def seq(self) -> str:
        """
        Sequence.seq is a read-only property, can only be modified by mutation
        """
        return self._seq

    @property
    def weight(self) -> float:
        """
        Calculate the molecular Weight
        Returns:
            weight: Sequence's Mole Weight in Daltons
        """
        if not self._weight:
            weight_table = config.MW[self.__class__.__name__ + "_MW"]
            self._weight = round(sum([weight_table[e]
                                 for e in self._seq], 18), 2)
        return self._weight

    def align(self,
              subject: Union[str, "Sequence"],
              mode: int = 1) \
            -> Tuple[str, str, float]:
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
            raise TypeError(
                f"Only str or {self.__class__.__name__} can be aligned to {self.__class__.__name__}")

        args = (self._seq, subject,
                AlignmentConfig.MATCH, AlignmentConfig.MISMATCH,
                AlignmentConfig.GAP_OPEN,  AlignmentConfig.GAP_EXTEND)

        if mode == 1:
            return algorithm.NeedlemanWunsch(*args)
        elif mode == 2:
            return algorithm.SmithWaterman(*args)
        else:
            raise TypeError("\
                Please choose alignment mode:\n\
                1-Global alignment by Needleman-Wunsch\n\
                2-Local alignment by Smith-Waterman")

    def find(self, target: Union[str, "Sequence"]) -> List[int]:
        """
        Find the target sequence in sequence and return the position
        Returns:
            All position of target appearance in self.sequence
        """
        if isinstance(target, Sequence):
            target = target._seq

        elif not isinstance(target, str):
            raise TypeError("Target should be str or Sequence")

        return [i.start() for i in re.finditer(target, self._seq)]

    def mutation(self,
                 position: Union[str, int, List[int]],
                 target: Union[str, "Sequence"]) -> str:
        """
        Modify the sequence, reset the cached property.
        Args:
            position: char, index(s) to be mutation
            target: the target char of mutation
        Returns:
            seq: Sequence after modified
        """
        if isinstance(target, Sequence):
            target = target._seq

        elif not(target, str):
            raise TypeError("Target should be str or Sequence")

        if isinstance(position, str):
            self._seq = self._seq.replace(position, target)
            self.reset_cache()  # reset cached property
            return self._seq

        elif isinstance(position, int):
            position = [position]

        elif not isinstance(position, List):
            raise TypeError("Position should be str, int or list of int")

        seq_list, length = list(self._seq), len(target)
        prev_end = 0    # previous end of mutation

        while position:
            pos = position.pop(0)
            if pos < prev_end:
                print(
                    f"WARNING: muation <{pos}~{pos + length}> overlaped previous mutation <{prev_end - length}~{prev_end}>")

            if not isinstance(pos, int):
                raise TypeError("Position must be int")

            seq_list[pos: pos + len(target)] = target
            prev_end = pos + length
        self._seq = "".join(seq_list)

        self.reset_cache()  # reset cached property
        return self._seq

    def toDNA(self) -> "DNA":
        """
        Convert Sequence to a DNA sequence
        """
        return DNA(self.seq, self.info)

    def toRNA(self) -> "RNA":
        """
        Convert Sequence to a RNA sequence
        """
        return DNA(self.seq, self.info)

    def toPeptide(self) -> "Peptide":
        """
        Convert Sequence to a Peptide
        """
        return Peptide(self.seq, self.info)

    @abstractmethod
    def _print(self) -> str:
        """
        Output sequence info
        """
        return self._seq

    def __add__(self, s: Union[str, "Sequence"]) -> "Sequence":
        if isinstance(s, str):
            return self.__class__(self._seq + s)
        elif isinstance(s, self.__class__):
            return self.__class__(self._seq + s._seq)
        else:
            raise TypeError(
                f"Only str or {self.__class__.__name__} can be added to {self.__class__.__name__}")

    def __radd__(self, s: Union[str, "Sequence"]) -> "Sequence":
        if isinstance(s, str):
            return self.__class__(s + self._seq)
        elif isinstance(s, self.__class__):
            return self.__class__(s._seq + self._seq)
        else:
            raise TypeError(
                f"Only str or {self.__class__.__name__} can add {self.__class__.__name__}")

    def __eq__(self, o: object) -> bool:
        if isinstance(o, str):
            return self._seq == o
        elif isinstance(o, self.__class__):
            return self._seq == o._seq
        else:
            raise TypeError(
                f"Only str or {self.__class__.__name__} can compare to {self.__class__.__name__}")

    def __getitem__(self, index: int) -> str:
        return self._seq[index]

    def __len__(self) -> int:
        return len(self._seq)

    def __str__(self) -> str:
        return self._print()

    def __repr__(self) -> str:
        return self._print()


class Peptide(Sequence):
    _pI: float
    _Hphob_list: List[float]

    def reset_cache(self):
        super().reset_cache()
        self._pI = 0.
        self._Hphob_list = []

    @property
    def pI(self) -> float:
        """
        Calculate the peptide's pI which is a pH make peptide's charge equal zero
        """
        if not self._pI:
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
            self._pI = round(pH, 3)
        return self._pI

    def chargeInpH(self, pH: float) -> float:
        """
        Calculate the charge amount of peptide at pH
        Henderson Hasselbalch equation: pH = pKa + log([A-]/[HA])
        Rearranging: [HA]/[A-] = 10 ** (pKa - pH)
        partial_charge =
            [A-]/[A]total = [A-]/([A-] + [HA]) = 1 / { ([A-] + [HA])/[A-] } =
            1 / (1 + [HA]/[A-]) = 1 / (1 + 10 ** (pKa - pH)) for acidic residues;
                                  1 / (1 + 10 ** (pH - pKa)) for basic residues
        """
        pos_charge = 1.0 / (10 ** (pH - config.PK["Nterm"]) + 1.0)
        for pos_aa, aa_pK in config.PK["pos_pK"].items():
            pos_charge += self.composition.get(pos_aa, 0) / \
                (10 ** (pH - aa_pK) + 1.0)

        neg_charge = 1.0 / (10 ** (config.PK["Cterm"] - pH) + 1.0)
        for neg_aa, aa_pK in config.PK["neg_pK"].items():
            neg_charge += self.composition.get(neg_aa, 0) / \
                (10 ** (aa_pK - pH) + 1.0)

        return pos_charge - neg_charge

    def getHphob(self,
                 window_size: int = 9,
                 show_img: bool = False) -> List[float]:
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
        if not self._Hphob_list:
            _Hphob = [config.HYDROPATHY[aa] for aa in self._seq]
            half_part = window_size // 2

            if len(_Hphob) < half_part:
                print("WARNING: Windows size is to big to calculate")

            for i in range(half_part, len(_Hphob) - half_part):
                self._Hphob_list.append(
                    round(sum(_Hphob[i - half_part: i + window_size - half_part]) / window_size, 3))

        if show_img:
            try:
                import matplotlib.pyplot as plt
            except ImportError:
                print(
                    "Can't import matplotlib.pyplot, please use 'pip install matplotlib' to install.")
            else:
                plt.title(
                    f"Hydropathy Score for {self._seq[:4]}...{self._seq[-4:]}")
                plt.plot(range(window_size // 2 + 1, self.length -
                               window_size // 2 + 1), self._Hphob_list, linewidth=0.8)
                plt.xlim((1, self.length + 1))
                plt.grid(linestyle="--")
                plt.xlabel("Position")
                plt.ylabel("Score")
                plt.show()

        return self._Hphob_list

    def _print(self) -> str:
        if self._seq:
            return f"N-{super()._print()}-C"
        return self.__class__.__name__ + "()"


T = TypeVar("T", "RNA", "DNA")

class RNA(Sequence):
    _GC: float
    orf: List[str]
    peptide: List[Peptide]

    def reset_cache(self):
        self._GC, self.orf = 0., []
        super().reset_cache()

    def complement(self):
        """
        Change the sequence to complemented
        """
        self._seq = self.complemented._seq

    @property
    def complemented(self: T) -> T:
        """
        Return complemented sequence
        """
        complement_dict = config.NC_INFO[self.__class__.__name__ + "_COMPLEMENT"]
        complement_seq = "".join([complement_dict[bp]
                                 for bp in self._seq[::-1]])
        return self.__class__(complement_seq)

    @property
    def GC(self) -> float:
        """
        Calculate the GC percentage
        """
        if not self._GC:
            self._GC = round(
                (self.composition["C"] + self.composition["G"]) / self.length, 4)
        return self._GC

    def reverse(self):
        """
        Change the sequence to reversed
        """
        self._seq = self.reversed._seq

    @property
    def reversed(self: T) -> T:
        """
        Return reversed sequence
        """
        return self.__class__(self._seq[::-1])

    def getOrf(self,
               multi: bool = False,
               replace: bool = False) -> List[str]:
        """
        Find the Open Reading Frame in sequence and save in self.orf
        Args:
            multi: Return all Orf if true else only the longest
            replace: Replace origin sequence with the longest Orf
        Returns:
            Orf: Sequence's Orf got by search
        """
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

    def transcript(self, filtered: bool) -> List[Peptide]:
        """
        Transcript the sequence to peptide, the result will save in self.peptide
        Args:
            filtered:  Return all product or the longest
        Returns:
            peptide: List of transcript product
        """
        orf = self.orf if self.orf else self.getOrf(multi=not filter)
        self.peptide = []
        for frame in orf:
            peptide = ""
            for i in range(0, len(frame), 3):
                peptide += config.TABLE.get(frame[i:i + 3],
                                            config.SYMBOL["printAlign"][1])
            self.peptide.append(Peptide(peptide[:-1]))

            if filtered:
                return self.peptide

        return self.peptide

    def _print(self):
        if self._seq:
            return f"5'-{super()._print()}-3'"
        return self.__class__.__name__ + "()"


class DNA(RNA):
    _translate: Optional[RNA]

    def reset_cache(self):
        self._translate = None
        super().reset_cache()

    def translate(self) -> RNA:
        if not self._translate:
            self._translate = RNA(self._seq.replace("T", "U"))
        return self._translate

    def transcript(self, filtered: bool = True) -> List[Peptide]:
        return self.translate().transcript(filtered)
