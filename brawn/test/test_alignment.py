# pylint: disable=missing-function-docstring,missing-module-docstring,missing-class-docstring
# pylint: disable=protected-access

from io import StringIO
import json
import unittest
from unittest.mock import patch

from brawn.alignment import (
    Alignment,
    Alphabet,
    MismatchedCacheVersionError,
    combine_alignments,
)

from .shared import SMALL_DB, SMALL_QUERY, dict_to_fasta, near_equal


class TestAlignment(unittest.TestCase):
    def test_equality(self):
        first = Alignment({"A": "AA", "AA": "CC"})
        # names match, sequences don't
        second = Alignment({"A": "BB", "AA": "DD"})
        assert first != second
        # names don't match, sequences do
        second = Alignment({"B": "AA", "BB": "CC"})
        assert first != second
        # both sequences and names match
        second = Alignment({"A": "AA", "AA": "CC"})
        assert first == second

    def test_from_file(self):
        from_f = Alignment.from_file(StringIO(dict_to_fasta(SMALL_DB)))
        assert isinstance(from_f, Alignment)
        assert from_f.to_dict() == SMALL_DB
        assert from_f.column_count == 6

    def test_from_dict(self):
        from_d = Alignment(SMALL_DB)
        assert isinstance(from_d, Alignment)
        assert from_d.to_dict() == SMALL_DB
        assert from_d.column_count == 6

    def test_empty(self):
        with self.assertRaisesRegex(ValueError, "at least one"):
            Alignment({})

    def test_inconsistent_lengths(self):
        with self.assertRaisesRegex(ValueError, "not of consistent length"):
            Alignment({"A": "AA", "B": "BBB"})

    def test_bad_weights(self):
        weights = Alignment(SMALL_DB).get_weights()
        with self.assertRaisesRegex(ValueError, "does not match number of sequences"):
            Alignment(SMALL_DB, weights=weights[:-1])  # trimmed to different length

    def test_bad_positions(self):
        good = Alignment(SMALL_DB)
        weights = good.get_weights()
        positions = good.positions[:-1]  # trimmed to different length
        with self.assertRaisesRegex(ValueError, "does not match sequence length"):
            Alignment(SMALL_DB, weights=weights, positions=positions[:-1])

    def test_weight_building(self):
        alignment = Alignment(SMALL_DB)
        with patch.object(Alignment, "_build_weights", return_value=[0, 1]) as patched:
            assert alignment._weights is None  # shouldn't be built
            assert not patched.called  # sanity check, shouldn't be touched from the above
            assert alignment.get_weights() == alignment._weights == [0, 1]
            assert patched.called  # definitely should have been called now

        # and now to test caching
        with patch.object(Alignment, "_build_weights", side_effect=ValueError) as patched:
            assert alignment._weights  # sanity check that it's still built
            assert alignment.get_weights() == alignment._weights == [0, 1]

    def test_sequence_fetching(self):
        alignment = Alignment({"A": "AA", "C": "CC"})
        assert alignment.get_sequence_by_name("A") == ["A", "A"]
        assert alignment.get_sequence_by_name("C") == ["C", "C"]

    def test_pid_zeroes(self):
        alignment = Alignment({"A": "A-", "B": "-B"})
        assert alignment.get_percentage_identity_pair(0, 1) == 0

    def test_to_file(self):
        handle = StringIO()
        alignment = Alignment(SMALL_DB)
        alignment.to_file(handle=handle)
        handle.seek(0)
        assert handle.read() == ">A\nGT-DVG\n>B\nGTK-VG\n"

    def test_to_file_with_columns(self):
        handle = StringIO()
        alignment = Alignment(SMALL_DB)
        alignment.to_file(handle=handle, columns=4)
        expected = ">A\nGT-D\nVG\n>B\nGTK-\nVG\n"
        handle.seek(0)
        assert handle.read() == expected

    def test_combining_mismatching_alphabets(self):
        first = Alignment({"A": "A"})
        second = Alignment({"C": "C"})
        second._alphabet = Alphabet.DNA
        with self.assertRaisesRegex(ValueError, "alphabets must be the same"):
            combine_alignments(first, second)

    def test_weighted_counts(self):
        alignment = Alignment({"A": "BA-", "B": "AZX"})
        alignment._weights = [0.2, 0.8]
        result = alignment.get_fractional_weighted_counts(0)
        expected = [0.] * 20
        expected[0] = 0.8  # A in second sequence
        expected[2] = 0.1  # from the ambiguous B in first sequence
        expected[11] = 0.1  # from the ambiguous B in first sequence
        assert result == expected

        result = alignment.get_fractional_weighted_counts(1)
        expected = [0.] * 20
        expected[0] = 0.2  # A in first sequence
        expected[3] = 0.4  # from the ambiguous Z in second sequence
        expected[13] = 0.4  # from the ambiguous Z in second sequence
        assert result == expected

        result = alignment.get_fractional_weighted_counts(2)
        # the X in the first sequence, spread evenly and zero from the gap in the first
        expected = [0.05] * 20
        assert near_equal(result, expected)

    def test_weighted_counts_dna(self):
        alignment = Alignment({"A": "RA-", "B": "AYN"})
        alignment._alphabet = Alphabet.DNA
        alignment._weights = [0.2, 0.8]
        result = alignment.get_fractional_weighted_counts(0)
        expected = [0.1, 0.0, 0.9, 0.0]  # R from the first spreads
        assert result == expected

        result = alignment.get_fractional_weighted_counts(1)
        expected = [0.0, 0.4, 0.2, 0.4]  # Y from the first spreads
        assert result == expected

        result = alignment.get_fractional_weighted_counts(2)
        expected = [0.05] * 4  # the N spreads over 20 rather than 4, matching MUSCLE
        assert near_equal(result, expected)

    def test_caching(self):
        reference = Alignment(SMALL_DB)
        handle = StringIO()
        reference.to_cache_file(handle)

        normal_result = combine_alignments(Alignment(SMALL_QUERY), reference).to_dict()

        handle.seek(0)
        cached_reference = Alignment.from_cache_file(handle)
        cached_result = combine_alignments(Alignment(SMALL_QUERY), cached_reference).to_dict()

        assert normal_result == cached_result

        handle.seek(0)
        as_json = json.load(handle)
        as_json["version"] = "different"
        with patch.object(json, "loads", return_value=as_json):
            with self.assertRaisesRegex(MismatchedCacheVersionError, "outdated cache version"):
                Alignment.from_cache_file(handle)

    def test_from_bad_files(self):
        handle = StringIO(">A\nAAA\n>B\n")
        with self.assertRaisesRegex(ValueError, "missing sequence for B"):
            Alignment.from_file(handle)

        handle = StringIO(">A\n>B\nBBB")
        with self.assertRaisesRegex(ValueError, "missing sequence for A"):
            Alignment.from_file(handle)

        handle = StringIO("AAA\n>B\nBBB")
        with self.assertRaisesRegex(ValueError, "sequence without name"):
            Alignment.from_file(handle)

        handle = StringIO(">A\nAAA\n>B\n>C\nCCC\n")
        with self.assertRaisesRegex(ValueError, "missing sequence for B"):
            Alignment.from_file(handle)

        with self.assertRaisesRegex(TypeError, "handle must be readable"):
            Alignment.from_file("")

    def test_columns(self):
        for seq in ["AAA", "AAAAA", "A"]:
            alignment = Alignment({"A": seq})
            assert alignment.column_count == len(seq)


def test_gap_open_weights():
    alignment = Alignment({"a": "-ERF", "b": "M-RF", "c": "-E--"}, weights=[0.7, 0.2, 0.1])
    weights = [alignment.get_gap_open_weight_total(i) for i in range(alignment.column_count)]
    assert near_equal(weights, [
        .8,  # opened in both first and last sequence
        .2,  # opened only in second
        .1,  # opened only in third
        0.,  # no new gap is opened, the third sequence gap extends
    ])


def test_gap_close_weights():
    alignment = Alignment({"a": "-ERF", "b": "M-RF", "c": "-E--"}, weights=[0.7, 0.2, 0.1])
    weights = [alignment.get_gap_close_weight_total(i) for i in range(alignment.column_count)]
    assert near_equal(weights, [.8, .2, 0, .1])
    assert near_equal(weights, [
        .8,  # a one-length gap closes in both first and last sequence
        .2,  # a one-length gap closes in second
        0.,  # a new gap opens in third, but extends beyond this column
        .1,  # the gap in third closes
    ])
