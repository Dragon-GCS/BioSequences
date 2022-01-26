from typing import Tuple


def NeedlemanWunsch(query: str,
                    subject: str,
                    match: float,
                    mismatch: float,
                    gap_open: float,
                    gap_extend: float) -> Tuple[str, str, float]:
    ...


def SmithWaterman(query: str,
                  subject: str,
                  match: float,
                  mismatch: float,
                  gap_open: float,
                  gap_extend: float) -> Tuple[str, str, float]:
    ...
