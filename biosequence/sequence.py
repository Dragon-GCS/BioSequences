# Base class Sequence for DNA,RNA,peptide...
import re

from collections import Counter
from typing import Tuple, Union, List
from abc import ABC, abstractmethod
from matplotlib import pyplot as plt

from biosequence.align import alignment
from biosequence import config

class Sequence(ABC):
    def __init__(self, seq: str = "") -> None:
        self._seq = seq.upper()

    @property
    def seq(self) -> str:
        """
        Seq can only be modified by mutation
        """
        return self._seq

    @property
    def length(self) -> int:
        return len(self._seq)

    @property
    def weight(self) -> float:
        """
        Calculate the molecular Weight
        Returns:
            weight: Sequence's Mole Weight in Daltons
        """
        if not hasattr(self, "_weight"):
            weight_table = config.MW[self.__class__.__name__ + "_MW"]
            self._weight = round(sum([weight_table[e] for e in self._seq], 18),2)
        return self._weight

    def align(self, subject: Union[str, "Sequence"], mode: int = 1, return_score: bool = False) -> Tuple[str, str, None] :
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

    def analysis(self) -> dict:
        """
        Analysis the composition of sequence
        Returns:
            composition: Dict of each element's appearance's times
        """
        if not hasattr(self, "_composition"):
            for key in sorted(counter:=dict(Counter(self._seq))):
                self._composition[key] = counter[key]
        return self._composition

    def find(self, target:str) -> List[int]:
        """
        Find the target sequence in sequence and return the position
        Returns:
            : All position of target appearence in self.seuquence
        """
        return [i.start() for i in re.finditer(target, self._seq)]

    def mutation(self, position: Union[int, List[int], str], target: str) -> str:
        """
        Modify the seqence
        Args:
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

        seq_list = list(self._seq)
        for pos in position:
            seq_list[pos: pos+len(target)] = target
        self._seq = "".join(seq_list)

        return self._seq

    @abstractmethod
    def _print(self) -> str:
        """
        Output sequence info
        """
        ...

    def __str__(self) -> str:
        return self._print()

    def __repr__(self) -> str:
        return self._print()

    def __add__(self, s: Union[str, "Sequence"]) -> "Sequence":
        """
        Add a new sequence to the end and return a new Sequence
        """
        if isinstance(s, str):
            return self.__class__(self._seq + s)
        elif isinstance(s, self.__class__):
            return self.__class__(self._seq + s._seq)
        else:
            raise TypeError(f"Only str or {self.__class__.__name__} can be added to {self.__class__.__name__}")

    def __radd__(self, s: Union[str, "Sequence"]) -> "Sequence":
        """
        Add a new sequence to the start and return a new Sequence
        """
        if isinstance(s, str):
            return self.__class__(s + self._seq)
        elif isinstance(s, self.__class__):
            return self.__class__(s._seq + self._seq)
        else:
            raise TypeError(f"Only str or {self.__class__.__name__} can be added to {self.__class__.__name__}")


class Peptide(Sequence):
    @property
    def pI(self):
        if not hasattr(self, "_pl"):
            pass
    

    def getHphob(self, window_size: int = 9, show_img: bool = True):
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
        if not hasattr(self, "_Hphob_lsit"):
            _Hphob = [config.HYDROPATHY[aa] for aa in self._seq]
            half_part = window_size // 2
            self._Hphob_lsit = []
            for i in range(half_part, len(_Hphob) - half_part):
                self._Hphob_lsit.append(round(sum(_Hphob[i - half_part: i + window_size - half_part]) / window_size, 3))
            
        if show_img:
            plt.title(f"Hydropathy Score for {self._seq[:4]}...{self._seq[-4:]}")
            plt.plot(range(window_size // 2 + 1, self.length - window_size // 2 + 1), self._Hphob_lsit, linewidth=0.8)
            plt.xlim((1, self.length + 1))
            plt.grid(linestyle="--")
            plt.xlabel("Position")
            plt.ylabel("Score")
            plt.show()

        return self._Hphob_lsit
    
    def _print(self) -> str:
        return f"N'-{self._seq}-C'"


class RNA(Sequence):
    @property
    def reverse(self):
        """
        Reverse the sequence
        """
        self._seq = self._seq.reversed()

    @property 
    def complement(self):
        "将该序列变为其互补序列并返回序列"
        pass

    @property 
    def GC(self):
        " 返回序列中GC含量,计算后保存在_GC中"
        pass
    def reversed(self):
        "计算返回序列的反向序列"
        pass

    def complemented(self):
        "计算返回序列的互补序列"
        pass

    def getOrf(self, multi=False, replace=False):
        "multi是否查找所有frame +1~+3的orf，默认值为仅查找最长的orf。 replace 当multi=False是生效，是否使用最长的orf替换原序列"
        if not hasattr(self, "orf"):
            pass
        
    def transcript(self, filter):
        "filter是否仅返回最长的翻译产物。返回值为一个或多个Peptide对象。"
        if not hasattr(self, "peptide"):
            pass
    
    def _print(self):
        return f"5'-{self._seq}-3'"


class DNA(RNA):
    def translate(self):
        if not hasattr(self, "_translate"):
            self._translate = RNA(self._seq.replace("T", "U"))
            
        return self._translate

    def transcript(self, filter=True):
        return self.translate().transcript(filter)