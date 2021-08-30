# Base class Sequence for DNA,RNA,peptide...
from typing import Union, Any

class Sequence():
    def __init__(self, seq:str ="")-> None: 
        self._seq = seq
        self._length = len(seq)
    
    @property
    def seq(self):
        """
        Seq can only be modified by mutation
        """
        return self._seq
    
    def mutation():
        pass
    
    def _print(self) -> str:
        """
        Output sequence info
        """
        return f"Seq({self._seq})"

    def __str__(self) -> str:
        return self._print()
    
    def __repr__(self) -> str:
        return self._print()
    
    def __add__(self, s:Union[str, Any]):
        """
        Add a new sequence to the end and return a new Sequence
        """
        if isinstance(s, str):
            return self.__class__(self._seq + s)
        elif isinstance(s, self.__class__):
            return self.__class__(self._seq + s._seq)
        else:
            raise TypeError(f"Only str or {self.__class__.__name__} can be add to {self.__class__.__name__}")
    
    def __radd__(self, s:Union[str, Any]):
        """
        Add a new sequence to the start and return a new Sequence
        """
        if isinstance(s, str):
            return self.__class__(s + self._seq)
        elif isinstance(s, self.__class__):
            return self.__class__(s._seq+self._seq)
        else:
            raise TypeError(f"Only str or {self.__class__.__name__} can be add to {self.__class__.__name__}")

