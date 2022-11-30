# pylint: disable=missing-function-docstring,missing-module-docstring,missing-class-docstring
# pylint: disable=protected-access

import argparse
from io import StringIO
from unittest.mock import patch

from brawn.alignment import Alignment
from brawn import _cmd as cmd
from brawn._cmd import (
    _main,
    _main_build_cache,
    _swap_muscle_args as swap_muscle_args,
)

from .shared import SMALL_DB, SMALL_QUERY, dict_to_fasta


def test_muscle_args():
    # checks that mimicing MUSCLE args works
    args = swap_muscle_args(["executable", "-profile", "-quiet", "-in1",
                             "some.fasta", "-in2", "some_ref.fa"])
    assert args == ["executable", "some.fasta", "--reference-alignment", "some_ref.fa"]


@patch.object(cmd, "_main", side_effect=ValueError("stop here"))
def test_muscle_args_via_parse(_patched_main):
    muscle_args = ["executable", "-profile", "-quiet", "-in1", "some.fasta", "-in2", "some_ref.fa"]
    brawn_args = ["executable", "some.fasta", "--reference-alignment", "some_ref.fa"]

    # for some reason "called_once_with" always passes on this check, so do it the long way
    namespace = argparse.Namespace()
    namespace.build_cache = False
    namespace.query = None
    namespace.reference_alignment = None
    namespace.output_columns = 1

    def patch_parse(args=None):
        assert args == brawn_args  # check the args are correct
        return namespace

    with patch.object(argparse.ArgumentParser, "parse_args", side_effect=patch_parse):
        try:
            cmd.entrypoint(muscle_args)
            assert False  # the mock didn't function
        except ValueError as err:
            assert "stop here" in str(err)  # only the mocked value is ok


def test_main_build_cache():
    in_handle = StringIO(dict_to_fasta(SMALL_DB))
    out_handle = StringIO()
    status = _main_build_cache(in_handle, out_handle)
    assert status == 0
    out_handle.seek(0)
    assert Alignment.from_cache_file(out_handle) == Alignment(SMALL_DB)


def test_main():
    ref_handle = StringIO(dict_to_fasta(SMALL_DB))
    query_handle = StringIO(dict_to_fasta(SMALL_QUERY))
    out_handle = StringIO()
    status = _main(query_handle, ref_handle, -1, out_handle)
    assert status == 0
    out_handle.seek(0)
    result = out_handle.read()
    assert result == ">query\nGT--IV\n>A\nGT-DVG\n>B\nGTK-VG\n"

    # and if a cache file is given instead, it should match
    ref_handle = StringIO()
    Alignment(SMALL_DB).to_cache_file(ref_handle)
    ref_handle.seek(0)
    query_handle.seek(0)
    out_handle = StringIO()

    status = _main(query_handle, ref_handle, -1, out_handle)
    assert status == 0
    out_handle.seek(0)
    cache_result = out_handle.read()
    assert result == cache_result


def test_main_invalid_ref():
    ref_handle = StringIO("{>the worst of both: worlds}")
    query_handle = StringIO(dict_to_fasta(SMALL_QUERY))
    status = _main(query_handle, ref_handle, -1, StringIO())
    assert status != 0


def test_main_invalid_query():
    query_handle = StringIO("{>the worst of both: worlds}")
    ref_handle = StringIO(dict_to_fasta(SMALL_DB))
    status = _main(query_handle, ref_handle, -1, StringIO())
    assert status != 0
