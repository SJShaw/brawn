""" Contains various alphabet handling code """

from enum import auto, Enum

UNKNOWN_CHAR = "?"

RESIDUE_GROUP = {  # via MAFFT with some MUSCLE modifications
    "A": 0,
    "C": 5,
    "D": 2,
    "E": 2,
    "F": 4,
    "G": 0,
    "H": 3,
    "I": 1,
    "K": 3,
    "L": 1,
    "M": 1,
    "N": 2,
    "P": 0,
    "Q": 2,
    "R": 3,
    "S": 0,
    "T": 0,
    "V": 1,
    "W": 4,
    "Y": 4,
    "B": 2,
    "Z": 2,
    "X": 0,
    None: 0,  # gaps
}

AMINO_INDICES = {k: i for i, k in enumerate(RESIDUE_GROUP)}
AMINO_INDICES[chr(0)] = len(RESIDUE_GROUP) + 1
AMINO_INDICES[UNKNOWN_CHAR] = len(RESIDUE_GROUP) + 2

DNA_INDICES = {k: i for i, k in enumerate("GCAT")}


class Alphabet(Enum):
    """ Possible sequence alphabets """
    AMINO = auto()
    DNA = auto()
    RNA = auto()

    def __str__(self) -> str:
        return super().__str__().split(".")[1]

    @classmethod
    def from_string(cls, string: str) -> "Alphabet":
        """ Selects an alphabet from various string options """
        string = string.upper()
        try:
            return Alphabet[string]
        except KeyError:
            if string in ["AA", "PROT", "PROTEIN"]:
                return Alphabet.AMINO
        raise ValueError(f"unknown alphabet: '{string}'")

    def get_size(self) -> int:
        """ Returns the number of residues in the given alphabet """
        if self == Alphabet.AMINO:
            return 20
        if self in [Alphabet.RNA, Alphabet.DNA]:
            return 4
        raise NotImplementedError(f"unhandled alphabet type: '{self}'")


def is_gap_char(char: str) -> bool:
    """ Returns True if the given character is a gap character """
    return char in "-."


def is_wildcard(char: str, alphabet: Alphabet) -> bool:
    """ Returns True if the given character is a wildcard/ambiguity character
        in the given alphabet
    """
    if char is None:
        return False
    if alphabet == Alphabet.AMINO:
        return char in "BXZ"
    if alphabet == Alphabet.DNA:
        return char not in "GATC"
    raise NotImplementedError(f"unhandled alphabet type: '{alphabet}'")


def is_residue(char: str, alphabet: Alphabet) -> bool:
    """ Returns True if the given character is a residue in the given alphabet
    """
    if not char:
        return False
    if alphabet == Alphabet.AMINO:
        return char in "ACDEFGHIKLMNPQRSTVWY"
    if alphabet == Alphabet.DNA:
        return char in "GATC"
    raise NotImplementedError(f"unhandled alphabet type: '{alphabet}'")
