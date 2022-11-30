""" A collection of classes and functions relating to multiple sequence alignments
    and operations between them
"""
import json
import os
import sys
from typing import Any, IO, Iterator, Optional, Union

from .align import global_align
from .alphabet import (
    AMINO_INDICES,
    DNA_INDICES,
    Alphabet,
    is_gap_char,
    is_residue,
    is_wildcard,
)
from .constants import (
    AMINO_SCORE_MATRIX,
    OTHER_SCORE_MATRIX,
    ID_GUARD,
    LENGTH_GUARD,
    VERSION,
)
from .distance import (
    calculate_distance,
)
from .path import Modification, Path
from .positions import AlignmentPosition, AlignmentPositions
from .tree import Tree


class MismatchedCacheVersionError(ValueError):
    """ Used when a cache file load is attempted from a mismatching version """


class InvalidCacheFormatError(ValueError):
    """ Used when a cache file load is attempted from a file with an invalid format """


class Sequence(list):
    """ A wrapper for sequences to handle conversion to/from gap characters
    """
    @classmethod
    def from_string(cls, line: str, alphabet: Alphabet) -> "Sequence":
        """ Constructs an instance from a string and an alphabet """
        wildcard = "X" if alphabet == Alphabet.AMINO else "N"
        seq: Sequence = cls()
        for char in line:
            if is_gap_char(char):
                seq.append(None)
            elif char and not is_residue(char, alphabet) and not is_wildcard(char, alphabet):
                seq.append(wildcard)
            else:
                seq.append(char)
        return seq

    def __str__(self) -> str:
        return "".join([s or "-" for s in self])


class Alignment:  # pylint: disable=too-many-instance-attributes
    """ Represents a multiple sequence alignment.

        Arguments:
            sequence_by_name: a mapping of entry name to entry sequence
            alphabet: the alphabet of the sequences
    """
    def __init__(self, sequence_by_name: dict[str, str], *,
                 tree: Tree = None, weights: list[float] = None,
                 positions: AlignmentPositions = None) -> None:
        if not sequence_by_name:
            raise ValueError("at least one sequence must be provided")
        self._col_count = len(list(sequence_by_name.values())[0])
        if not all(len(seq) == self._col_count for seq in sequence_by_name.values()):
            raise ValueError("alignment sequences not of consistent length")
        if weights and len(weights) != len(sequence_by_name):
            raise ValueError("number of weights does not match number of sequences")
        if positions and len(positions) != self._col_count:
            raise ValueError("number of positions does not match sequence length")
        self._alphabet = Alphabet.AMINO
        self._sequence_by_name = {}
        for name, sequence in sequence_by_name.items():
            self._sequence_by_name[name] = Sequence.from_string(sequence, self.alphabet)
        self._names = tuple(self._sequence_by_name)
        self._sequences = tuple(self._sequence_by_name.values())

        # cache sections
        self._tree: Optional[Tree] = tree
        self._weights: Optional[list[float]] = weights
        self._positions: Optional[AlignmentPositions] = positions

    @property
    def alphabet(self) -> Alphabet:
        """ The alphabet used in the alignment """
        return self._alphabet

    @property
    def column_count(self) -> int:
        """ The number of columns (i.e. the length of a sequence) in
            the alignment
        """
        return self._col_count

    @property
    def names(self) -> tuple[str, ...]:
        """ The names of the sequences in the alignment """
        return self._names

    @property
    def sequences(self) -> tuple[Sequence, ...]:
        """ The sequences of the alignment """
        return self._sequences

    def get_sequence_by_name(self, name: str) -> Sequence:
        """ Returns the sequence in the alignment that matches the given name
        """
        return self._sequence_by_name[name]

    @property
    def positions(self) -> AlignmentPositions:
        """ Returns the positional data for each position in the alignment
        """
        if self._positions is None:
            self._positions = tuple(self._build_positions())
        return self._positions

    def _build_weights(self) -> list[float]:
        if self._weights is not None:
            return self._weights
        if self._tree is None:
            self._tree = tree_from_alignment(self)
        return self._tree.get_weights()

    def _build_positions(self) -> Iterator[AlignmentPosition]:
        if self.alphabet == Alphabet.AMINO:
            score_matrix = AMINO_SCORE_MATRIX
        else:
            score_matrix = OTHER_SCORE_MATRIX
        alphabet_size = self.alphabet.get_size()
        for column in range(self.column_count):
            counts = self.get_fractional_weighted_counts(column)
            sort_order = get_indices_by_decreasing_value(counts)
            scores = []
            for i in range(alphabet_size):
                score = 0.
                for count, matrix_value in zip(counts, score_matrix[i]):
                    score += count * matrix_value
                scores.append(score)
            pos = AlignmentPosition(
                ungapped_weight=1.0 - self.get_column_ungapped_weight(column),
                gap_opens=1.0 - self.get_gap_open_weight_total(column),
                gap_closes=1.0 - self.get_gap_close_weight_total(column),
                base_counts=counts, sort_order=sort_order, scores=scores,
            )
            yield pos

    def get_weights(self) -> list[float]:
        """ Gets a list of sequence weights matching the sequence order"""
        if self._weights is None:
            self._weights = self._build_weights()
        return list(self._weights)

    def get_sequence_weight(self, index: int) -> float:
        """ Returns the tree weight of the sequence with the given index
        """
        if self._weights is None:
            self._weights = self._build_weights()
        return self._weights[index]

    def get_percentage_identity_pair(self, i: int, j: int) -> float:
        """ Returns the percentage of similar characters in the two sequences,
            skipping any pair which include a gap
        """
        count = 0
        same = 0
        for first, second in zip(self._sequences[i], self._sequences[j]):
            if first and second:  # neither can be a gap
                count += 1
                if first == second:
                    same += 1
        if not count:
            result = 0.
        else:
            result = same / count
        return result

    @classmethod
    def from_file(cls, handle: IO) -> "Alignment":
        """ Constructs an alignment instance from a file handle with
            content in FASTA format
        """
        if not hasattr(handle, "read"):
            raise TypeError("file handle must be readable")
        names: list[str] = []
        seqs: list[list[str]] = []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if seqs and not seqs[-1]:
                    raise ValueError(f"alignment missing sequence for {names[-1]}")
                names.append(line[1:])
                seqs.append([])
            elif line:
                if not seqs:
                    raise ValueError("sequence without name in alignment")
                seqs[-1].append(line)
        if len(names) != len(seqs):  # the checks above should catch this, but just in case
            raise ValueError("mismatching sequence and name count in file")
        for name, seq in zip(names, seqs):
            if not seq:
                raise ValueError(f"alignment missing sequence for {name}")
        combined = ["".join(seq) for seq in seqs]
        return cls(dict(zip(names, combined)))

    @classmethod
    def from_cache_file(cls, handle: IO) -> "Alignment":
        """ Rebuilds an alignment instance from a cache file """
        try:
            data = json.load(handle)
        except json.decoder.JSONDecodeError as err:
            raise InvalidCacheFormatError() from err
        if data["version"] != VERSION:
            msg = f"loading outdated cache version ({data['version']}) in version {VERSION}"
            raise MismatchedCacheVersionError(msg)

        positions = tuple(AlignmentPosition.from_json(pos) for pos in data["positions"])
        alignment = cls(data["sequences"], weights=data["weights"], positions=positions)
        return alignment

    def to_cache_file(self, handle: IO) -> None:
        """ Saves a cached version of alignment to the given file """
        positions = [pos.to_json() for pos in self.positions]
        assert self._weights
        data = {
            "alphabet": str(self.alphabet),
            "sequences": dict(zip(self.names, map(str, self.sequences))),
            "weights": self._weights,
            "positions": positions,
            "version": VERSION,
        }
        json.dump(data, handle)

    def get_column_ungapped_weight(self, column: int) -> float:
        """ Returns the total sequence weight of all sequences which have a gap in the
            given column.
        """
        total = 0.
        for i, seq in enumerate(self._sequences):
            if not seq[column]:
                total += self.get_sequence_weight(i)
        return total

    def get_fractional_weighted_counts(self, col_index: int  # pylint: disable=too-many-branches
                                       ) -> list[float]:
        """ Returns the residue weight totals of the given column,
            summing for each residue the sequence weight of the sequence
            contributing each residue
        """
        alphabet_size = self.alphabet.get_size()
        counts = [0.] * alphabet_size
        total_weight = 0.
        for seq, weight in zip(self._sequences, self.get_weights()):
            char = seq[col_index]
            if not char:
                continue

            if is_wildcard(char, self.alphabet):
                if self.alphabet == Alphabet.AMINO:
                    if char == "B":
                        counts[AMINO_INDICES["D"]] += weight / 2
                        counts[AMINO_INDICES["N"]] += weight / 2
                    elif char == "Z":
                        counts[AMINO_INDICES["E"]] += weight / 2
                        counts[AMINO_INDICES["Q"]] += weight / 2
                    else:
                        # completely unhandled, apply it equally to everything in the alphabet
                        avg = weight / alphabet_size
                        counts = [c + avg for c in counts]
                elif self.alphabet in [Alphabet.DNA, Alphabet.RNA]:
                    if char == "R":
                        counts[DNA_INDICES["G"]] += weight / 2
                        counts[DNA_INDICES["A"]] += weight / 2
                    elif char == "Y":
                        counts[DNA_INDICES["C"]] += weight / 2
                        counts[DNA_INDICES["T"]] += weight / 2
                    else:
                        # completely unhandled, apply it equally to everything in the alphabet
                        avg = weight / 20  # should be 4, but this matches MUSCLE
                        counts = [c + avg for c in counts]
                else:
                    raise NotImplementedError(f"unhandled alphabet: {self.alphabet}")
            else:
                if self.alphabet == Alphabet.AMINO:
                    counts[AMINO_INDICES[char]] += weight
                else:
                    counts[DNA_INDICES[char]] += weight
            total_weight += weight

        if total_weight:
            assert total_weight <= 1.001, f"{total_weight} {self._weights}"
            counts = [c / total_weight for c in counts]

        return counts

    def get_gap_open_weight_total(self, column: int) -> float:
        """ Returns the total sequence weight of all sequences which have a gap opening in
            the given column.
        """
        total = 0.
        if column < 1:  # no previous column to check
            for i, seq in enumerate(self._sequences):
                if seq[column] is None:
                    total += self.get_sequence_weight(i)
            return total

        for i, seq in enumerate(self._sequences):
            if seq[column] is None and seq[column - 1] is not None:
                total += self.get_sequence_weight(i)
        return total

    def get_gap_close_weight_total(self, column: int) -> float:
        """ Returns the total sequence weight of all sequences which have a gap closing in
            the given column.
        """
        total = 0.
        if self.column_count - 1 == column:
            for i, seq in enumerate(self._sequences):
                if not seq[column]:
                    total += self.get_sequence_weight(i)
            return total

        for i, seq in enumerate(self._sequences):
            if not seq[column] and seq[column + 1]:
                total += self.get_sequence_weight(i)
        return total

    def to_file(self, *, columns: int = -1, handle: IO = sys.stdout) -> None:
        """ Writes the alignment to file in FASTA format

            Arguments:
                columns: the maximum number of columns per line
                handle: the file to write to if not stdout
        """
        assert self.column_count != 0
        seqs = list(map(str, self._sequences))
        if columns < 1:
            columns = self.column_count
        for name, seq in zip(self.names, seqs):
            handle.write(f">{name}")
            handle.write(os.linesep)
            for i in range(0, self.column_count, columns):
                handle.write(seq[i:i + columns])
                handle.write(os.linesep)

    def to_dict(self) -> dict[str, str]:
        """ Returns a dictionary mapping sequence name to sequence for
            each sequence in the alignment
        """
        seqs = list(map(str, self._sequences))
        return dict(zip(self.names, seqs))

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, Alignment) and self._sequences == other.sequences and \
               self.names == other.names and self.alphabet == other.alphabet


class ResultAlignment:
    """ The result of combining two existing alignments via a path.
        Members are generally caculated in a lazy fashion, on demand.
    """
    def __init__(self, path: Path, query_alignment: Alignment, reference_alignment: Alignment
                 ) -> None:
        self.path = path
        self._query = query_alignment
        self._reference = reference_alignment

        self._query_built: list[Optional[str]] = [None] * len(self._query.sequences)
        self._reference_built: list[Optional[str]] = [None] * len(self._reference.sequences)

    @property
    def column_count(self) -> int:
        """ The number of columns in the resulting alignment """
        return len(self.path.edges)

    @property
    def names(self) -> Iterator[str]:
        """ The names of the inputs, starting with those from the query """
        for name in self._query.names:
            yield name
        for name in self._reference.names:
            yield name

    @property
    def sequences(self) -> Iterator[str]:
        """ The newly aligned sequences, starting with those from the query """
        for i, seq in enumerate(self._query.sequences):
            if self._query_built[i] is None:
                self._query_built[i] = build_query_result(seq, self.path)
            result = self._query_built[i]
            assert result
            yield result
        for i, seq in enumerate(self._reference.sequences):
            if self._reference_built[i] is None:
                self._reference_built[i] = build_reference_result(seq, self.path)
            result = self._reference_built[i]
            assert result
            yield result

    def get_aligned_references(self) -> dict[str, str]:
        """ Returns a dictionary of the newly (re-)aligned reference sequences """
        sequences = []
        for i, seq in enumerate(self._reference.sequences):
            if self._reference_built[i] is None:
                self._reference_built[i] = build_reference_result(seq, self.path)
            result = self._reference_built[i]
            assert result
            sequences.append(result)
        return dict(zip(self._reference.names, sequences))

    def to_file(self, handle: IO = sys.stdout, columns: int = 60) -> None:
        """ Writes a FASTA file of the result

            Arguments:
                handle: the file handle to write to (defaults to stdout)
                columns: the maximum columns per line for sequences

            Returns:
                None
        """
        if columns < 1:
            columns = len(self.path.edges)
        for name, seq in zip(self.names, self.sequences):
            handle.write(f">{name}")
            handle.write(os.linesep)
            for i in range(0, len(seq), columns):
                handle.write(seq[i:i + columns])
                handle.write(os.linesep)

    def to_dict(self) -> dict[str, str]:
        """ Returns a dictionary mapping name to sequence for each of the inputs """
        names = list(self.names)
        result = dict(zip(names, self.sequences))
        if len(result) != len(names):
            raise ValueError("cannot build dictionary, sequence names are not unique")
        return result


def get_indices_by_decreasing_value(values: list[Union[int, float]]) -> list[int]:
    """ Returns a list of indices, sorted by decreasing value in the
        given list
    """
    return [i for i, v in sorted(enumerate(values),
                                 key=lambda x: (x[1], -x[0]),
                                 reverse=True)]


def tree_from_alignment(  # pylint: disable=too-many-statements,too-many-locals,too-many-branches
                        alignment: Alignment) -> Tree:
    """ Builds a tree of sequences from the given alignment

        Arguments:
            alignment: the alignment to build a tree from

        Returns:
            the resulting Tree
    """
    # all positions are leaves, which means N-1 internal nodes in the tree
    leaf_count = len(alignment.sequences)
    internal_node_count = leaf_count - 1

    # a 1D representation of the triangular matrix
    distances = [0.] * ((leaf_count * internal_node_count) // 2)

    def get_flat_index(i: int, j: int) -> int:
        """ Converts a 2D index for a triangular matrix into a 1D index
            into an array
        """
        if i >= j:
            return i * (i - 1) // 2 + j
        return j * (j - 1) // 2 + i

    node_indices = list(range(leaf_count))
    nearest_neighbours = [ID_GUARD] * leaf_count
    min_dists = [LENGTH_GUARD] * leaf_count

    lefts = [ID_GUARD] * internal_node_count
    rights = [ID_GUARD] * internal_node_count
    parents = [ID_GUARD] * internal_node_count
    heights = [LENGTH_GUARD] * internal_node_count
    left_lengths = [LENGTH_GUARD] * internal_node_count
    right_lengths = [LENGTH_GUARD] * internal_node_count

    for i in range(1, leaf_count):
        row_start = get_flat_index(i, 0)
        for j in range(i):
            pid = alignment.get_percentage_identity_pair(i, j)
            distances[row_start + j] = calculate_distance(pid)
        for j in range(i):
            distance = distances[row_start + j]
            if distance < min_dists[i]:
                min_dists[i] = distance
                nearest_neighbours[i] = j
            if distance < min_dists[j]:
                min_dists[j] = distance
                nearest_neighbours[j] = i

    for internal_node_index in range(internal_node_count):
        left_min = ID_GUARD
        right_min = ID_GUARD
        min_dist = LENGTH_GUARD
        for j in range(leaf_count):
            if node_indices[j] == ID_GUARD:
                continue
            distance = min_dists[j]
            if distance < min_dist:
                min_dist = distance
                left_min = j
                right_min = nearest_neighbours[j]

        assert left_min != ID_GUARD
        assert right_min != ID_GUARD

        new_min_dist = LENGTH_GUARD
        new_nearest = ID_GUARD
        for j in range(leaf_count):
            if j in (left_min, right_min):
                continue
            if node_indices[j] == ID_GUARD:
                continue
            left_index = get_flat_index(left_min, j)
            distance_left = distances[left_index]
            distance_right = distances[get_flat_index(right_min, j)]
            new_dist = (.1 * ((distance_left + distance_right) / 2)
                        + .9 * min(distance_left, distance_right))
            if nearest_neighbours[j] == right_min:
                nearest_neighbours[j] = left_min

            distances[left_index] = new_dist
            if new_dist < new_min_dist:
                new_min_dist = new_dist
                new_nearest = j

        assert internal_node_index < leaf_count - 1 or new_min_dist != LENGTH_GUARD
        assert internal_node_index < leaf_count - 1 or new_nearest != ID_GUARD

        new_height = distances[get_flat_index(left_min, right_min)] / 2
        left = node_indices[left_min]
        right = node_indices[right_min]
        height_left = 0 if left < leaf_count else heights[left - leaf_count]
        height_right = 0 if right < leaf_count else heights[right - leaf_count]

        lefts[internal_node_index] = left
        rights[internal_node_index] = right
        left_lengths[internal_node_index] = new_height - height_left
        right_lengths[internal_node_index] = new_height - height_right
        heights[internal_node_index] = new_height

        node_indices[left_min] = leaf_count + internal_node_index
        nearest_neighbours[left_min] = new_nearest
        min_dists[left_min] = new_min_dist

        node_indices[right_min] = ID_GUARD

    node_count = 2 * leaf_count - 1
    root = node_count - 1

    lefts = [ID_GUARD] * leaf_count + lefts
    rights = [ID_GUARD] * leaf_count + rights
    left_lengths = [LENGTH_GUARD] * leaf_count + left_lengths
    right_lengths = [LENGTH_GUARD] * leaf_count + right_lengths

    parents = [ID_GUARD] * node_count
    parent_lengths = [LENGTH_GUARD] * node_count
    for i in range(leaf_count, node_count):
        left = lefts[i]
        right = rights[i]
        parents[left] = i
        parents[right] = i
        parent_lengths[left] = left_lengths[i]
        parent_lengths[right] = right_lengths[i]
    return Tree(
        node_count=node_count, root_node_index=root,
        lefts=lefts, rights=rights,
        left_lengths=left_lengths, right_lengths=right_lengths,
        names=list(alignment.names),
        parent_lengths=parent_lengths, parents=parents,
    )


def build_query_result(sequence: Sequence, path: Path) -> str:
    """ Aligns a query sequence to match the alignment path

        Arguments:
            sequence: the sequence to align
            path: the path of the alignment

        Returns:
            the aligned sequence as a string
    """
    result: list[str] = []
    char = iter(sequence)
    for edge in path.edges:
        if edge.type == Modification.MATCH:
            result.append(next(char) or "-")
        elif edge.type == Modification.INSERTION:
            result.append("-")
        else:
            assert edge.type == Modification.DELETION
            result.append(next(char) or "-")
    return "".join(result)


def build_reference_result(sequence: Sequence, path: Path) -> str:
    """ Aligns a reference sequence to match the alignment path

        Arguments:
            sequence: the sequence to align
            path: the path of the alignment

        Returns:
            the aligned sequence as a string
    """
    result: list[str] = []
    char = iter(sequence)
    for edge in path.edges:
        if edge.type == Modification.MATCH:
            result.append(next(char) or "-")
        elif edge.type == Modification.INSERTION:
            result.append(next(char) or "-")
        else:
            assert edge.type == Modification.DELETION
            result.append("-")
    return "".join(result)


def combine_alignments(query: Alignment, reference: Alignment) -> ResultAlignment:
    """ Aligns the sequences in the query and reference alignments,
        keeping individual columns of each alignment unchanged.

        Arguments:
            query: the query alignment
            reference: the reference alignment

        Returns:
            the combined alignment
    """

    if query.alphabet != reference.alphabet:
        raise ValueError("alignment alphabets must be the same")

    path = global_align(
        query.positions,
        reference.positions,
        reference.alphabet,
    )
    return ResultAlignment(path, query, reference)


def insert_into_alignment(query_sequence: str, alignment: Alignment,
                          ) -> tuple[str, Union[str, dict[str, str]]]:
    """ Inserts the given query sequence into an existing alignment.

        Arguments:
            query_sequence: the sequence to insert
            alignment: the alignment to insert the query into

        Returns:
            a tuple of
                the aligned query sequence
                a dictionary mapping reference name to reference sequence OR single sequence
    """
    query = Alignment({"query": query_sequence})
    result = combine_alignments(query, alignment)
    references_aligned = result.get_aligned_references()
    query_aligned = build_query_result(query.sequences[0], result.path)
    return query_aligned, references_aligned


def get_aligned_pair(query_sequence: str, reference_name: str, alignment: Alignment,
                     ) -> tuple[str, Union[str, dict[str, str]]]:
    """ Inserts the given query sequence into an existing alignment and calculates
        the aligned sequences of the query and the named reference

        Arguments:
            query_sequence: the sequence to insert
            reference_name: the name of the reference to return the newly aligned sequence of
            alignment: the alignment to insert the query into

        Returns:
            a tuple of
                the aligned query sequence
                the aligned reference sequence of the named reference
    """
    try:
        ref_seq = alignment.get_sequence_by_name(reference_name)
    except KeyError as err:
        raise ValueError("reference of interest not in reference alignment:"
                         f"{reference_name}") from err
    query = Alignment({"query": query_sequence})
    result = combine_alignments(query, alignment)
    path = result.path
    query_seq = Sequence.from_string(query_sequence, alignment.alphabet)
    query_align = build_query_result(query_seq, path)
    ref_align = build_reference_result(ref_seq, path)
    return query_align, ref_align
