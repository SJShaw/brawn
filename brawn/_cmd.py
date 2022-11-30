""" Command line helpers for brawn """

import argparse
import sys
from typing import IO

from .alignment import (
    Alignment,
    InvalidCacheFormatError,
    combine_alignments,
)


def _main(query_file: IO, reference_file: IO, output_columns: int,
          output_file: IO = sys.stdout) -> int:
    try:
        try:
            reference = Alignment.from_cache_file(reference_file)
        except InvalidCacheFormatError:
            reference_file.seek(0)
            reference = Alignment.from_file(reference_file)
    except ValueError as err:
        print(f"Unknown reference alignment format: {err}", file=sys.stderr)
        return 1

    try:
        query_alignment = Alignment.from_file(query_file)
    except ValueError as err:
        print(f"Invalid query format: {err}", file=sys.stderr)
        return 1
    result = combine_alignments(query_alignment, reference)
    try:
        result.to_file(handle=output_file, columns=output_columns)
    except BrokenPipeError:  # if the user doesn't want more, don't raise an error
        pass
    return 0


def _main_build_cache(input_file: IO, output_file: IO) -> int:
    try:
        alignment = Alignment.from_file(input_file)
        alignment.to_cache_file(output_file)
    except ValueError as err:
        print(f"Could not build cache file: {err}", file=sys.stderr)
        return 1
    return 0


def _swap_muscle_args(args: list[str]) -> list[str]:
    args = list(args)  # since some modification is about to take place
    # convert MUSCLE args
    if set(args).intersection({"-profile", "-in1", "-in2"}):
        # the optional quiet mode is irrelevant
        if "-quiet" in args:
            args.pop(args.index("-quiet"))
        # profile-profile mode is the only part replicated, so drop it
        args.pop(args.index("-profile"))
        # convert the query sequence/alignment into the default arg for argparse
        args.pop(args.index("-in1"))
        # and convert the reference/second alignment into the matching arg
        if "-in2" in args:
            args[args.index("-in2")] = "--reference-alignment"
    return args


def _parse_args(args: list[str]) -> argparse.Namespace:
    args = _swap_muscle_args(args)

    parser = argparse.ArgumentParser()
    parser.add_argument("query", metavar="FASTA",
                        type=argparse.FileType("r"))
    parser.add_argument("--reference-alignment", metavar="DB_FASTA",
                        type=argparse.FileType("r"))
    parser.add_argument("--build-cache", metavar="CACHE_PATH", default=None,
                        type=argparse.FileType("w"))
    parser.add_argument("--output-columns", default=60, type=int)
    return parser.parse_args(args)


def entrypoint(arguments: list[str] = None) -> int:
    """ The main brawn entry point, if no arguments are given, sys.argv will be used """
    if arguments is None:
        arguments = list(sys.argv)
    args = _parse_args(arguments)
    if args.build_cache:
        return _main_build_cache(args.query, args.build_cache)

    return _main(args.query, args.reference_alignment, args.output_columns)
