from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Iterator, Iterable, List, Literal, Tuple, Union, overload
from urllib.parse import urlencode
from urllib.request import urlopen
from urllib.error import HTTPError

from bioseq.config import SYMBOL
from bioseq import DNA, RNA, Peptide, Sequence


# TODO: Merge fetch function
@overload
def fetchENS(uid: str) -> DNA:
    ...


@overload
def fetchENS(uid: Iterable[str]) -> List[DNA]:
    ...


def fetchENS(uid):
    """
    Fetch sequence corresponding to UID from Ensemble REST api.

    Args:
        uid: One or list of ENS's unique id
    Returns:
        One of list of DNA
    """
    ENS_DATABASE_URL = "https://rest.ensembl.org/sequence/id/{}?content-type=text/plain"

    def fetch(uid):
        raw_info = ""
        for _ in range(3):
            try:
                raw_info = urlopen(ENS_DATABASE_URL.format(uid)).read().decode()
                break
            except HTTPError as e:
                if e.code == 400:
                    print(f"{uid} not found")
                    break
                else:
                    print(f"{uid} retry: {e}.")
            except Exception as e:
                print(e)
                break
        return DNA(raw_info, uid)

    if isinstance(uid, str):
        return fetch(uid)
    elif isinstance(uid, Iterable):
        with ThreadPoolExecutor(5) as executor:
            return list(executor.map(fetch, uid))
    else:
        raise ValueError(f"{uid} is not a str or list of str")


@overload
def fetchNCBI(uid: str) -> Union[DNA, RNA, Peptide]:
    ...


@overload
def fetchNCBI(uid: Iterable[str]) -> List[Sequence]:
    ...


def fetchNCBI(uid):
    """
    Fetch sequence corresponding to UID from NCBI E-utilities. Only support RNA, mRNA(DNA), Protein.

    NM_(mRNA):      Protein-coding transcripts (usually curated)
    NR_(RNA ):      Non-protein-coding transcripts
    XM_(mRNA):      Predicted model protein-coding transcript
    XR_(RNA ):      Predicted model non-protein-coding transcript
    AP_(Protein):   Annotated on AC_ alternate assembly
    NP_(Protein):   Associated with an NM_ or NC_ accession
    YP_(Protein):   Annotated on genomic molecules without an instantiated transcript record
    XP_(Protein):   Predicted model, associated with an XM_ accession
    WP_(Protein):   Non-redundant across multiple strains and species

    Args:
        uid: One or list of NCBI's unique id
    Returns:
        If uid is a list, the return is a list of Sequence(excluded the uid not found data on NCBI)
        without ensure sequence's type,else the return is a Sequence corresponding to UID.
    """
    EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    EUTILS_POST = {
        "db": "",               # database
        "rettype": "fasta",     # text type
        "retmode": "text",      # data type
        "id": "",               # uid
    }

    def checkUID(uid: str) -> Tuple[str, str]:
        """
        Check whether the uid is valid
        Args:
            uid: uid
        Returns:
            returns: ncbi_db_name, Sequence's type corresponding to uid
        """
        if uid[:3] in ["AP_", "NP_", "YP_", "XP_", "WP_"]:
            return "protein", "Peptide"
        elif uid[:3] in ["NM_", "XM_"]:
            return "nuccore", "DNA"
        elif uid[:3] in ["NR_", "XR_"]:
            return "nuccore", "RNA"
        else:
            raise ValueError(f"{uid} is not a support uid")

    def fetch(data: Dict) -> List[Sequence]:
        try:
            raw_info = urlopen(
                EUTILS_URL + urlencode(data)
            ).read().decode()
        except HTTPError as e:
            print(e)
            return [Sequence()]
        else:
            return parseFasta(raw_info)

    if isinstance(uid, str):
        EUTILS_POST["id"] = uid
        EUTILS_POST["db"], seq_type = checkUID(uid)
        sequence = fetch(EUTILS_POST)[0]
        return getattr(sequence, f"to{seq_type}")()
    elif isinstance(uid, Iterable):
        uids = defaultdict(list)
        for id in uid:
            db, seq_type = checkUID(id)
            uids[db].append((id, seq_type))
        result = []
        for db, seqs in uids.items():
            EUTILS_POST["db"] = db
            for i in range(0, len(seqs), 200):
                EUTILS_POST["id"] = ",".join([info[0] for info in seqs[i:i + 200]])
                result.extend(fetch(EUTILS_POST))
        return result
    else:
        raise ValueError(f"{uid} is not a str or list of str")


def parseFasta(fasta_text: str) -> List[Sequence]:
    """
    Parse a string in FASTA format

    Args:
        fasta_text: string to be parsed
    Returns:
        A list[Sequence] of parsing result
    """
    return [Sequence("".join((text := fasta.split("\n"))[1:]), text[0])
            for fasta in fasta_text.split(">")
            if fasta.strip()]


@overload
def loadFasta(filename: str) -> List[Sequence]:
    ...


@overload
def loadFasta(filename: str, iterator: Literal[True]) -> Iterator[Sequence]:
    ...


def loadFasta(filename, iterator=False):
    """
    Read fasta file

    Args:
        filename: the fasta file's name.
        iterator: Set to True as reading a large file, it will return a iterator.
    Returns:
        A Iterator(when iterator=True) or a List of Sequence.
    """
    def iter_parse(filename: str) -> Iterator[Sequence]:
        seq = ""
        info = ""
        with open(filename, encoding="utf8") as f:
            while line := f.readline():
                # process the blank line
                if not line.strip() and seq:
                    yield Sequence(seq, info)
                    seq = ""    # reset seq
                    continue

                if line.startswith(">"):
                    # yield previous sequence
                    if seq:
                        yield Sequence(seq, info)
                        seq = ""    # reset seq
                    info = line[1:].strip()
                else:
                    seq += line.strip()

            if seq:
                yield Sequence(seq, info)

    if iterator:
        return iter_parse(filename)
    else:
        with open(filename, encoding="utf8") as f:
            return parseFasta(f.read())


def printAlign(
        sequence1: str,
        sequence2: str,
        spacing: int = 10,
        line_width: int = 30,
        show_seq: bool = True):
    """
    Print two sequence by a pretty format

    Args:
        sequence1:  Sequence1
        sequence2:  Sequence2
        spacing:    A space each $spacing char
        line_width: the width of each line
        show_seq:   if False, only print the alignment result
    """
    symbol_line = ""
    format_seq1 = ""
    format_seq2 = ""
    count = 0
    length = len(sequence1) if len(sequence1) < len(
        sequence2) else len(sequence2)
    match_symbol, mismatch_symbol, gap_symbol = SYMBOL["printAlign"]

    for i in range(length):
        base1 = sequence1[i]
        base2 = sequence2[i]
        if base1 == base2:
            symbol_line += match_symbol
        elif base1 == "-" or base2 == "-":
            symbol_line += gap_symbol
        else:
            symbol_line += mismatch_symbol

        format_seq1 += base1
        format_seq2 += base2
        count += 1

        if count % spacing == 0:
            format_seq1 += " "
            format_seq2 += " "
            symbol_line += " "
        if count % line_width == 0:
            count = 0

    space_num = 0
    space_each_line = line_width // spacing

    for i in range(0, length, line_width):
        start = i + space_num
        end = start + line_width + space_each_line
        space_num += space_each_line
        if show_seq:
            print(f"{i + 1:>5}", end=" ")
            print(format_seq1[start: end])

            print(f" " * 5, end=" ")
            print(symbol_line[start: end])

            print(f"{i + 1:>5}", end=" ")
            print(format_seq2[start: end])
        else:
            print(f"{i + 1:>5}", end=" ")
            print(symbol_line[start: end])

        print()
