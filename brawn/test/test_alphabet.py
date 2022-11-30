# pylint: disable=missing-function-docstring,missing-module-docstring,missing-class-docstring
# pylint: disable=protected-access

from brawn.alphabet import Alphabet


def test_alphabet_from_string():
    for value in Alphabet:
        assert Alphabet.from_string(str(value)) == value
        assert Alphabet.from_string(str(value).lower()) == value
    assert Alphabet.from_string("PROT") == Alphabet.AMINO


def test_alphabet_conversion():
    assert str(Alphabet.AMINO) == "AMINO"
    for value in Alphabet:
        assert Alphabet.from_string(str(value)) == value
