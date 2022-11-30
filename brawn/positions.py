""" Classes and functions for handling data and calculations of particular
    positions/columns within alignments
"""

from dataclasses import dataclass
from math import log
from typing import Any

from .alphabet import Alphabet
from .constants import (
    GAP_OPEN,
    SCORE_CENTER,
    SCORE_GUARD,
)


@dataclass
class AlignmentPosition:  # pylint: disable=too-many-instance-attributes
    """ Data for a single character position within an alignment """
    sort_order: list[int]
    base_counts: list[float]
    scores: list[float]
    ungapped_weight: float
    gap_opens: float
    gap_closes: float
    score_gap_open: float = 0.
    score_gap_close: float = 0.

    def __post_init__(self) -> None:
        self.score_gap_open = self.gap_opens * GAP_OPEN / 2
        self.score_gap_close = self.gap_closes * GAP_OPEN / 2

    def to_json(self) -> dict[str, Any]:
        """ Convert the position a JSON-friendly representation """
        return dict(vars(self))

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "AlignmentPosition":
        """ Reconstructs a position from a JSON-friendly representation """
        return cls(**data)


AlignmentPositions = tuple[AlignmentPosition, ...]


def compare_profile_positions(query_positions: AlignmentPosition,
                              reference_positions: AlignmentPosition,
                              alphabet: Alphabet) -> float:
    """ Scores two positions from two different alignments
    """
    score = 0.
    for index in query_positions.sort_order:
        count = query_positions.base_counts[index]
        # since it's sorted by count, will always be zero after the first zero is found
        if not count:
            break
        score += count * reference_positions.scores[index]
    if alphabet == Alphabet.AMINO:
        if not score:
            score = -2.5
        else:
            score = (
                (log(score) - SCORE_CENTER)
                * query_positions.ungapped_weight
                * reference_positions.ungapped_weight
            )
    else:
        score -= SCORE_CENTER
    return score


def set_terminal_gaps(positions: AlignmentPositions) -> None:
    """ (Re)Sets the gap scores for the start/end positions of the given list """
    if not positions:
        return

    first = positions[0]
    last = positions[-1]

    if first.score_gap_open != SCORE_GUARD:
        first.score_gap_open = 0
    if len(positions) > 1 and last.score_gap_open != SCORE_GUARD:
        last.score_gap_close = 0
