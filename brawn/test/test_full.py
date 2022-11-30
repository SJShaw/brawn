# pylint: disable=missing-function-docstring,missing-module-docstring,missing-class-docstring
# pylint: disable=protected-access

from io import StringIO

from brawn.alphabet import (
    Alphabet,
    is_gap_char,
    is_residue,
    is_wildcard,
)
from brawn.alignment import (
    Alignment,
    combine_alignments,
    get_aligned_pair,
    insert_into_alignment,
)
from brawn.path import Modification

from .shared import (
    MEDIUM_DB_PATH,
    MEDIUM_OUTPUT_PATH,
    MEDIUM_QUERY_PATH,
    SMALL_DB,
    SMALL_QUERY,
)


def test_full_run_medium():
    with open(MEDIUM_QUERY_PATH, encoding="utf-8") as handle:
        query = Alignment.from_file(handle)
    with open(MEDIUM_DB_PATH, encoding="utf-8") as handle:
        ref = Alignment.from_file(handle)
    result = combine_alignments(query, ref)
    fasta = StringIO()
    result.to_file(handle=fasta, columns=60)
    with open(MEDIUM_OUTPUT_PATH, encoding="utf-8") as handle:
        expected = handle.read()
    fasta.seek(0)
    assert expected
    assert fasta.read() == expected


def test_full_run_small():
    result = combine_alignments(Alignment(SMALL_QUERY), Alignment(SMALL_DB))
    assert result.to_dict() == {
        "query": "GT--IV",
        "A": "GT-DVG",
        "B": "GTK-VG",
    }


def test_modification_to_string():
    assert str(Modification.MATCH) == "M"
    assert str(Modification.DELETION) == "D"


def test_dna_as_amino():
    query = "ATGCAGGCATGGCAGGTGCACGAGAACGGCGAGCCGGGCGAGGTGATGCGCCT"
    ref = {
        "a": "---ATGCGC---GCCGTACAGTTCGACCGTTTTGGCCCGCCCGACGTTCTGCGC",
        "b": "GTGATGCGT---GCCGTGCAGTTCGACCGGTACGGGGACCCCGACGTGCTG---",
        "c": "---ATGCCCAAGGCCGTAGCCATCCACCAGTTCGGCGGGCCGGACGTACTG---",
    }
    query_expected = "---ATGCAG---GCATGGCAGGTGCACGAGAACGGCGAGCCGGGCGAGGTGATGCGCCT"
    ref_expected = {
        "a": "---ATGCGC---GCCGTACAGTTCGACCGTTTTGGCCCGCCCGAC---GTTCTGCGC--",
        "b": "GTGATGCGT---GCCGTGCAGTTCGACCGGTACGGGGACCCCGAC---GTGCTG-----",
        "c": "---ATGCCCAAGGCCGTAGCCATCCACCAGTTCGGCGGGCCGGAC---GTACTG-----",
    }
    reference = Alignment(ref)
    aligned_query, aligned_reference = insert_into_alignment(query, reference)
    assert aligned_query == query_expected
    assert aligned_reference == ref_expected


def test_alphabet_sizes():
    assert Alphabet.AMINO.get_size() == 20
    assert Alphabet.DNA.get_size() == 4
    assert Alphabet.RNA.get_size() == 4


def test_residues():
    assert not is_residue(None, Alphabet.AMINO)
    for char in "ACDEFGHIKLMNPQRSTVWY":
        assert is_residue(char, Alphabet.AMINO)
    for char in "GCAT":
        assert is_residue(char, Alphabet.DNA)
    for char in "EFH":
        assert not is_residue(char, Alphabet.DNA)
    for char in "XO":
        assert not is_residue(char, Alphabet.AMINO)


def test_wildcard():
    assert not is_wildcard(None, Alphabet.AMINO)
    for char in "BXZ":
        assert is_wildcard(char, Alphabet.AMINO)
        assert is_wildcard(char, Alphabet.DNA)  # everything is a wildcard except ACGT
    for char in "ACGT":
        assert not is_wildcard(char, Alphabet.AMINO)
        assert not is_wildcard(char, Alphabet.DNA)


def test_gaps():
    assert is_gap_char("-")
    assert is_gap_char(".")
    assert not is_gap_char("A")
    assert not is_gap_char("N")


def test_insertion():
    with open(MEDIUM_QUERY_PATH, encoding="utf-8") as handle:
        query = Alignment.from_file(handle).to_dict()["query"]
    with open(MEDIUM_DB_PATH, encoding="utf-8") as handle:
        alignment = Alignment.from_file(handle)
    with open(MEDIUM_OUTPUT_PATH, encoding="utf-8") as handle:
        result_alignment = Alignment.from_file(handle)
    expected = result_alignment.to_dict()
    aligned_query = expected.pop("query")

    query_result, ref_align = insert_into_alignment(query, alignment)
    assert query_result == aligned_query
    assert ref_align == expected

    target = list(expected)[1]
    query_result, ref = get_aligned_pair(query, target, alignment)
    assert query_result == aligned_query
    assert ref == expected[target]

    # and just check that asking for a non-existent alignment fails as expected
    try:
        get_aligned_pair(query, "missing", alignment)
        assert False  # should have failed above
    except ValueError as err:
        assert "not in reference alignment" in str(err)
