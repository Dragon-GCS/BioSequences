from typing import List, Tuple, Union
from bioseq.sequence import Sequence

def read_fasta(filename:str) -> Tuple[List[str], List[str]]:    ...
def printAlign(sequence1: Union[str, Sequence], 
               sequence2:Union[str, Sequence],
               spacing:int = 10, 
               line_width:int = 30, 
               show_seq:bool = True) -> None:    ...
