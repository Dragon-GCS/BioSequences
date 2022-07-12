import re

from collections import Counter
from typing import Dict, List, Optional, Tuple, TypeVar, Union

from bioseq import config, algorithm
from bioseq.config import AlignmentConfig


class Sequence:
    info: str
    _seq: str
    _composition: Dict[str, Union[int, float]]
    _weight: float

    def __init__(self, seq: str, info: str = ""):
        """Base class of all sequence type

        Args:
            seq(str): Sequence
            info(str): Some information string about the sequence
        """
        self._seq = seq.upper()
        self.info = info
        self.reset_cache()

    def reset_cache(self):
        """
        Reset cached property related with sequence, called when sequence changed
        """
        self._composition = {}
        self._weight = 0.

    @property
    def composition(self) -> Dict[str, Union[int, float]]:
        """Analysis the composition of sequence

        Returns:
            Dict: Each element's appearance times or percentage in sequence
        """
        if not self._composition:
            counter = dict(Counter(self._seq))
            self._composition = {key: counter[key] for key in sorted(counter)}
        return self._composition

    @property
    def length(self) -> int:
        """
        The length of sequence
        """
        return len(self._seq)

    @property
    def seq(self) -> str:
        """
        Read-only property, sequence can only be modified by ``mutation()``
        """
        return self._seq

    @property
    def weight(self) -> float:
        """Calculate the sequence's Molar mass by below function

        .. math::
            weight_{seq} = \\sum_i^{length}weight_i - 18 * (length - 1)

        Returns:
            weight(float): Molar mass with unit of Dalton
        """
        if not self._weight:
            weight_table = config.MW[self.__class__.__name__]
            self._weight = sum(
                [weight_table[e] for e in self._seq]) - 18 * (self.length - 1)

        return self._weight

    def align(self,
              subject: Union[str, "Sequence"],
              mode: int = 1) \
            -> Tuple[str, str, float]:
        """Align two sequence. Use ``bioseq.config.AlignmentConfig`` to set the alignment score,
        including match(2), mismatch(-3), gap_open(-3), gap_extend(-3). number in brackets is default value

        Args:
            subject(str|Sequence): Sequence to align
            mode(int): 1: Use Needleman-Wunsch to global alignment\n
                  2: Use Smith-Waterman to partial alignment
        Returns:
            tuple:
            query(str): Self sequence after alignment\n
            subject(str): Subject sequence after alignment\n
            score（int): align score if choose return_score
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
        """Find the target sequence in this sequence and return the positions

        Returns:
            List[int]: All position of target in this sequence
        """
        if isinstance(target, Sequence):
            target = target._seq

        elif not isinstance(target, str):
            raise TypeError("Target should be str or Sequence")

        return [i.start() for i in re.finditer(target, self._seq)]

    def mutation(self,
                 position: Union[str, int, List[int]],
                 target: Union[str, "Sequence"]) -> str:
        """Change this sequence, ``self.seq`` can only be modified by this function

        Args:
            position(Union[str, int, List[int]): char, index(s) to be mutation
            target(Union[str, Sequence]): the target char of mutation
        Returns:
            str: Sequence after modified
        """
        if isinstance(target, Sequence):
            target = target._seq

        elif not isinstance(target, str):
            raise TypeError("Target should be str or Sequence")

        if isinstance(position, str):
            self._seq = self._seq.replace(position, target)
            # Reset cached property related with sequence(weight, composition...)
            self.reset_cache()
            return self._seq

        elif isinstance(position, int):
            position = [position]

        elif not isinstance(position, List):
            raise TypeError("Position should be str, int or list of int")

        seq_list, length = list(self._seq), len(target)
        prev_end = 0    # previous end of mutation

        while position:
            pos = position.pop(0)
            if not isinstance(pos, int):
                raise TypeError("Position must be int")

            if pos < prev_end:
                print(
                    f"WARNING: mutation <{pos}~{pos + length}> overlapped previous mutation <{prev_end - length}~{prev_end}>")

            prev_end = pos + length
            seq_list[pos: prev_end] = target
        self._seq = "".join(seq_list)

        # Reset cached property related with sequence(weight, composition...)
        self.reset_cache()
        return self._seq

    def toDNA(self) -> "DNA":
        """
        Convert to a DNA instance
        """
        return DNA(self.seq, self.info)

    def toRNA(self) -> "RNA":
        """
        Convert to an RNA instance
        """
        return RNA(self.seq, self.info)

    def toPeptide(self) -> "Peptide":
        """
        Convert to a Peptide instance
        """
        return Peptide(self.seq, self.info)

    def _print(self) -> str:
        """
        Sequence's output format, if length more than 30, only show the first and end 30 chars
        """
        if self._seq and len(self._seq) > 30:
            return f"{self._seq[:10]}...{self._seq[-10:]}"
        return self._seq

    def __add__(self, s: Union[str, "Sequence"]) -> "Sequence":
        if isinstance(s, str):
            return self.__class__(self._seq + s)
        elif isinstance(s, self.__class__):
            return self.__class__(self._seq + s._seq)
        else:
            raise TypeError(
                f"Expect str or {self.__class__.__name__} add with {self.__class__.__name__}, not {s.__class__.__name__}")

    def __radd__(self, s: Union[str, "Sequence"]) -> "Sequence":
        if isinstance(s, str):
            return self.__class__(s + self._seq)
        elif isinstance(s, self.__class__):
            return self.__class__(s._seq + self._seq)
        else:
            raise TypeError(
                f"Expect str or {self.__class__.__name__} add with {self.__class__.__name__}, not {s.__class__.__name__}")

    def __eq__(self, o: object) -> bool:
        if isinstance(o, str):
            return self._seq == o
        elif isinstance(o, self.__class__):
            return self._seq == o._seq
        else:
            raise TypeError(
                f"Expect str or {self.__class__.__name__} compare to {self.__class__.__name__}, not {o.__class__.__name__}")

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

        .. math::
            pH = pK_a + \\lg{\\frac{[A^-]}{[HA]}}\\implies\\frac{[HA]}{[A^-]} = 10^{pK_a - pH}

            charge = \\frac{[A^-]}{[A]_{total}} = \\frac{[A^-]}{[A^-] + [HA]}
            = \\frac{1}{\\frac{[A^-] + [HA]}{[A^-]}}
            = \\frac{1}{1 + \\frac{[HA]}{[A-]}}

        for acidic residues:  :math:`charge = 1 / (1 + 10 ^{pK_a - pH})`

        for basic residues:  :math:`charge = 1 / (1 + 10^{pH - pK_a})`

        Args:
            pH(float): pH value
        Returns:
            float: charge in specific pH
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
        Calculate the Hydropathy Score.The lager the score, the higher the hydrophobicity.
        Each aa's score is the average score of all aa in window_size.
        So part of Amino Acid at begin and end don't have score

        Args:
            window_size(int): the number for calculate average hydropathy value
            show_img(book): whether to draw the result, require ``matplotlib``
        Returns:
            List[float]: the result of peptide's Hydropathy Score
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
        """
        Peptide print starts with "N-" and then ends with "-C", means sequence is from N-terminal to C-terminal
        """
        if self._seq:
            return f"N-{super()._print()}-C"
        return self.__class__.__name__ + "()"


T = TypeVar("T", "RNA", "DNA")


class RNA(Sequence):
    _GC: float
    orf: List[str]          #: Can only visit after called `get_orf()`
    peptide: List[Peptide]  #: Can only visit after called `transcript()`

    def reset_cache(self):
        self._GC, self.orf, self.peptide = 0., [], []
        super().reset_cache()

    @property
    def complement(self: T) -> T:
        """
        Return self complementary sequence
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

    @property
    def reversed(self: T) -> T:
        """
        Return self reversed sequence
        """
        return self.__class__(self._seq[::-1])

    def getOrf(self,
               topn: int = 1,
               replace: bool = False) -> List[str]:
        """
        Find the Open Reading Frame in sequence and save in ``self.orf``

        Args:
            topn(int): the num of orfs, sorted by length of each orf, default is 1
            replace(bool): Replace origin sequence with the longest Orf
        Returns:
            List[str]: Orf found on self.
        """
        # Traverse all Orf Frame
        i, step = 0, 1
        starts: List[int] = []
        end_points: List[Tuple[int, int]] = []

        while i < self.length:
            if self._seq[i: i+3] in config.START_CODON:
                starts.append(i)
                step = 3

            i += step

            if config.CODON_TABLE.get(self._seq[i: i+3]) == "*":
                end_points.extend([(s, i + 3) for s in starts])
                starts = []
                step = 1

        end_points = sorted(
            end_points, key=lambda se: se[1] - se[0], reverse=True)
        self.orf = [self._seq[se[0]: se[1]] for se in end_points[:topn]]

        if replace:
            self._seq = self.orf[0]

        return self.orf

    def transcript(self, topn: int = 1) -> List[Peptide]:
        """
        Transcript the sequence to peptide, the result will save in ``self.peptide``

        Args:
            topn(int):  filter num of transcripts, sorted by length of each transcript, default is 1
        Returns:
            List[Peptide]: List of transcript product
        """
        orf = self.orf if self.orf else self.getOrf(topn=topn)
        self.peptide = []
        for frame in orf:
            peptide = "".join([
                config.CODON_TABLE.get(
                    frame[i:i + 3], config.SYMBOL["printAlign"][1])
                for i in range(0, len(frame), 3)
            ])

            self.peptide.append(Peptide(peptide[:-1]))

        return self.peptide

    def _print(self):
        """
        Peptide print starts with " 5'- " and then ends with " -3' ", means sequence is from 5' to 3'
        """
        if self._seq:
            return f"5'-{super()._print()}-3'"
        return self.__class__.__name__ + "()"


class DNA(RNA):
    _translate: Optional[RNA]

    def reset_cache(self):
        self._translate = None
        super().reset_cache()

    def translate(self) -> RNA:
        """
        Translate the sequence to RNA, which replace the T with U
        """
        if not self._translate:
            self._translate = RNA(self._seq.replace("T", "U"))
        return self._translate

    def getOrf(self, topn: int = 1, replace: bool = False) -> List[str]:
        """Return the open reading frame of mRNA which is translated from this sequence.

        Args:
            topn(int): the num of orfs, sorted by length of each orf, default is 1
            replace(bool): Replace origin sequence with the longest Orf
        Returns:
            List[str]: Orf found on mRNA
        """
        self.orf = self.translate().getOrf(topn, replace)
        return self.orf

    def transcript(self, topn: int = 1) -> List[Peptide]:
        """Return the transcript product of mRNA which is translated from this sequence.

        Args:
            topn(int):  filter num of transcripts, sorted by length of each transcript, default is 1
        Returns:
            List[Peptide]: List of transcript product
        """
        self.peptide = self.translate().transcript(topn)
        return self.peptide
