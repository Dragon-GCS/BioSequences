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
def fetchENS(uid: List[str]) -> List[DNA]:
    ...


def fetchENS(uid):
    """Fetch sequence corresponding to UID from Ensemble REST api.

    Args:
        uid(str | List[str]): One or list of ENS's unique id
    Returns:
        DNA | List[DNA]: One or list of DNA sequence corresponding to UID
    """
    ens_db_url = "https://rest.ensembl.org/sequence/id/{}?content-type=text/plain"

    def fetch(uid):
        raw_info = ""
        for i in range(3):
            try:
                print(
                    f"[Try {i + 1}/3]Fetching {uid} from Ensemble REST API...")
                raw_info = urlopen(ens_db_url.format(uid)).read().decode()
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
        raise ValueError(f"{uid}'s type <{type(uid)}> is not a str or list of str")


@overload
def fetchNCBI(uid: str) -> Union[DNA, RNA, Peptide]:
    ...


@overload
def fetchNCBI(uid: List[str]) -> List[Sequence]:
    ...


def fetchNCBI(uid):
    """Fetch sequence corresponding to UID from NCBI E-utilities. Only support RNA, mRNA(DNA), Protein.

    =============   ==============================================================================
    Prefix          Explanation
    =============   ==============================================================================
    NM_(mRNA)       Protein-coding transcripts (usually curated)
    NR_(RNA )       Non-protein-coding transcripts
    XM_(mRNA)       Predicted model protein-coding transcript
    XR_(RNA )       Predicted model non-protein-coding transcript
    AP_(Protein)    Annotated on AC alternate assembly
    NP_(Protein)    Associated with an NM or NC accession
    YP_(Protein)    Annotated on genomic molecules without an instantiated transcript record
    XP_(Protein)    Predicted model, associated with an XM accession
    WP_(Protein)    Non-redundant across multiple strains and species
    =============   ==============================================================================

    | NCBI RefSeq's document: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole
    | some NCBI E-utilities's api: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/

    Args:
        uid(str|List[str]): One or list of NCBI's unique id
    Returns:
        DNA | RNA | Peptide | List[Sequence]:
        If uid is a list, the return is a list of Sequence(excluded the uid not found data on NCBI)
        without ensure sequence's type, else the return is a Sequence corresponding to UID.
    """
    eutils_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    eutils_post = {
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
            print(f"Fetching {data['id']} from NCBI E-utilities...")
            raw_info = urlopen(
                eutils_url + urlencode(data)
            ).read().decode()
        except HTTPError as e:
            print(e)
            return [Sequence("")]
        else:
            return parseFasta(raw_info)

    if isinstance(uid, str):
        eutils_post["id"] = uid
        eutils_post["db"], seq_type = checkUID(uid)
        sequence = fetch(eutils_post)[0]
        return getattr(sequence, f"to{seq_type}")()
    elif isinstance(uid, Iterable):
        uids = {}
        for id in uid:
            db, seq_type = checkUID(id)
            uids.setdefault(db, []).append((id, seq_type))
        result = []
        for db, seqs in uids.items():
            eutils_post["db"] = db
            for i in range(0, len(seqs), 200):
                eutils_post["id"] = ",".join(
                    [info[0] for info in seqs[i:i + 200]])
                result.extend(fetch(eutils_post))
        return result
    else:
        raise ValueError(f"{uid}'s type <{type(uid)}> is not a str or list of str")


def parseFasta(fasta_text: str) -> List[Sequence]:
    """
    Parse a FASTA formatted string.

    Args:
        fasta_text(str): string to be parsed
    Returns:
        List[Sequence]: Parsing result
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
    """Load fasta file

    Args:
        filename: the fasta file's name.
        iterator: Set to True as reading a large file, it will return a iterator.
    Returns:
        List[Sequence] | Iterator
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
        sequence1(Iterator)
        sequence2(Iterator)
        spacing(int):    A space each $spacing char
        line_width(int): the width of each line
        show_seq(bool):   if False, only print the alignment result
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
