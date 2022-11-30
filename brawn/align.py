""" Functions for running a global alignment of two existing alignments """

from .alphabet import Alphabet
from .constants import (
    SCORE_GUARD,
    GAP_EXTEND,
)
from .positions import (
    AlignmentPositions,
    compare_profile_positions,
    set_terminal_gaps,
)
from .path import Edge, Modification, Path


def global_align(query_positions: AlignmentPositions,  # pylint: disable=too-many-statements,too-many-locals
                 reference_positions: AlignmentPositions,
                 alphabet: Alphabet) -> Path:
    """ Builds and returns a Path for the two given alignments being merged
        into a single alignment
    """
    set_terminal_gaps(query_positions)
    set_terminal_gaps(reference_positions)
    query_length = len(query_positions)
    reference_length = len(reference_positions)
    pref1 = query_length + 1
    pref2 = reference_length + 1
    current_match = [0.] * pref2
    next_match = [0.] * pref2
    prev_match = [SCORE_GUARD] * pref2
    delete_row = [SCORE_GUARD] * pref2
    traceback = [[0] * pref2 for i in range(pref1)]

    def recurse_d(row: list[int], i: int, j: int) -> None:
        dd = delete_row[j] + GAP_EXTEND
        md = prev_match[j] + query_positions[i - 1].score_gap_open
        if dd > md:
            delete_row[j] = dd
        else:
            delete_row[j] = md
            row[j] = (row[j] & ~BIT__D) | BIT_MD

    def recurse_i(iij: float, row: list[int], _i: int, j: int) -> float:
        iij += GAP_EXTEND
        mi = current_match[j - 1] + reference_positions[j - 1].score_gap_open
        if mi >= iij:
            iij = mi
            row[j] = (row[j] & ~BIT__I) | BIT_MI
        return iij

    def recurse_m(iij: float, i: int, j: int) -> None:
        dm = delete_row[j] + query_positions[i - 1].score_gap_close
        im = iij + reference_positions[j - 1].score_gap_close
        mm = current_match[j]
        if mm >= dm and mm >= im:
            next_match[j + 1] += mm
            bit = BIT_MM
        elif dm >= mm and dm >= im:
            next_match[j + 1] += dm
            bit = BIT_DM
        else:
            next_match[j + 1] += im
            bit = BIT_IM
        traceback[i + 1][j + 1] = (traceback[i + 1][j + 1] & ~BIT__M) | bit

    def set_bit_tbm(i: int, j: int, modification: Modification) -> None:
        if modification == Modification.MATCH:
            bit = BIT_MM
        elif modification == Modification.DELETION:
            bit = BIT_DM
        elif modification == Modification.INSERTION:
            bit = BIT_IM
        else:
            raise ValueError(f"unexpected modification type: {modification}")
        traceback[i][j] = ((traceback[i][j]) & ~BIT__M) | bit

    # start/init
    prev_match[0] = 0
    current_match[0] = SCORE_GUARD
    current_match[1] = compare_profile_positions(query_positions[0],
                                                 reference_positions[0], alphabet)
    set_bit_tbm(1, 1, Modification.MATCH)

    for j in range(2, pref2):
        current_match[j] = (
            compare_profile_positions(query_positions[0], reference_positions[j - 1],
                                      alphabet)
            + reference_positions[0].score_gap_open
            + (j - 2) * GAP_EXTEND
            + reference_positions[j-2].score_gap_close
        )
        set_bit_tbm(1, j, Modification.INSERTION)

    # mid
    for i in range(1, query_length):
        row = traceback[i]
        iij = SCORE_GUARD
        delete_row[0] = (
            query_positions[0].score_gap_open
            + (i - 1) * GAP_EXTEND
        )
        current_match[0] = SCORE_GUARD

        if i == 1:
            current_match[1] = compare_profile_positions(query_positions[0],
                                                         reference_positions[0],
                                                         alphabet)
            set_bit_tbm(i, 1, Modification.MATCH)
        else:
            current_match[1] = (
                compare_profile_positions(query_positions[i - 1],
                                          reference_positions[0], alphabet)
                + query_positions[0].score_gap_open
                + (i - 2) * GAP_EXTEND
                + query_positions[i - 2].score_gap_close
            )
            set_bit_tbm(i, 1, Modification.DELETION)

        for j in range(1, reference_length):
            next_match[j + 1] = compare_profile_positions(query_positions[i],
                                                          reference_positions[j],
                                                          alphabet)

        for j in range(1, reference_length):
            recurse_d(row, i, j)
            iij = recurse_i(iij, row, i, j)
            recurse_m(iij, i, j)

        recurse_d(row, i, reference_length)
        iij = recurse_i(iij, row, i, reference_length)
        prev_match, current_match, next_match = current_match, next_match, prev_match

    # final
    row = traceback[query_length]
    current_match[0] = SCORE_GUARD
    comparison = compare_profile_positions(query_positions[query_length - 1],
                                           reference_positions[0], alphabet)
    current_match[1] = comparison + query_positions[0].score_gap_open
    current_match[1] += ((query_length - 2) * GAP_EXTEND
                         + query_positions[query_length - 2].score_gap_close)
    set_bit_tbm(query_length, 1, Modification.DELETION)

    delete_row[0] = SCORE_GUARD
    for j in range(1, pref2):
        recurse_d(row, query_length, j)

    iij = SCORE_GUARD
    for j in range(1, pref2):
        iij = recurse_i(iij, row, query_length, j)

    dab = delete_row[reference_length]
    iab = iij

    score = current_match[reference_length]
    edge_type = Modification.MATCH

    if dab > score:
        score = dab
        edge_type = Modification.DELETION
    if iab > score:
        score = iab
        edge_type = Modification.INSERTION

    assert len(traceback) == query_length + 1, (query_length + 1, len(traceback))
    assert len(traceback[0]) == reference_length + 1
    return build_path(traceback, query_length, reference_length, edge_type)


# traceback matrix related code

BIT_MM = 0x00
BIT_DM = 0x01
BIT_IM = 0x02
BIT__M = 0x03
BIT_DD = 0x00
BIT_MD = 0x04
BIT__D = 0x04
BIT_II = 0x00
BIT_MI = 0x08
BIT__I = 0x08


def get_modification(bits: int, previous: Modification,  # pylint: disable=too-many-return-statements
                     ) -> Modification:
    """ Gets the matching modification instance from an  """
    if previous == Modification.MATCH:
        val = bits & BIT__M
        if val == BIT_MM:
            return Modification.MATCH
        if val == BIT_DM:
            return Modification.DELETION
        if val == BIT_IM:
            return Modification.INSERTION
        raise ValueError(f"incompatible matrix value for match: {val}")

    if previous == Modification.DELETION:
        val = bits & BIT__D
        if val == BIT_MD:
            return Modification.MATCH
        if val == BIT_DD:
            return Modification.DELETION
        raise ValueError(f"incompatible matrix value for deletion: {val}")

    if previous == Modification.INSERTION:
        val = bits & BIT__I
        if val == BIT_MI:
            return Modification.MATCH
        if val == BIT_II:
            return Modification.INSERTION
        raise ValueError(f"incompatible matrix value for insertion: {val}")

    raise ValueError(f"unknown modification type: '{previous}'")


def build_path(traceback: list[list[int]], query_length: int, reference_length: int,
               last_edge: Modification) -> Path:
    """ Builds a path through the matrix, beginning at the most distant point

        Arguments:
            traceback: the traceback matrix
            query_length: the number of columns in the query
            reference_length: the number of columns in the reference
            last_edge: the edge type to use for that last position

        Returns:
            a Path, with edges ordered from the beginning of the sequences
            to the end
    """
    edge = Edge(last_edge, query_length, reference_length)
    edges = [edge.copy()]
    while True:
        pla = edge.query_length
        plb = edge.reference_length
        bits = traceback[pla][plb]
        next_edge_type = get_modification(bits, edge.type)
        if edge.type == Modification.MATCH:
            edge.query_length -= 1
            edge.reference_length -= 1
        elif edge.type == Modification.DELETION:
            edge.query_length -= 1
        elif edge.type == Modification.INSERTION:
            edge.reference_length -= 1
        else:
            raise ValueError(f"unexpected modification type: {edge.type}")

        if not edge.query_length and not edge.reference_length:
            break

        edge.type = next_edge_type
        edges.append(edge.copy())
    return Path(edges[::-1])
