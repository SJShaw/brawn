""" A collection of useful data types for paths through a matrix """

from dataclasses import dataclass, field
from enum import auto, Enum


class Modification(Enum):
    """ A modification type for a sequence, based on the path direction """
    DELETION = auto()
    INSERTION = auto()
    MATCH = auto()

    def __str__(self) -> str:
        # don't use the full name, the first character is enough
        return super().__str__().split(".")[1][0]


@dataclass
class Edge:
    """ An collection of data of an edge in a path, with:
        type: the type of edge
        query_length: the number of steps to reach this in the query
        reference_length: the number of steps to reach this in the reference
    """
    type: Modification
    query_length: int
    reference_length: int

    def __post_init__(self) -> None:
        assert self.query_length >= 0
        assert self.reference_length >= 0

    def copy(self) -> "Edge":
        """ Creates a copy of the instance """
        return type(self)(self.type, self.query_length, self.reference_length)


@dataclass
class Path:
    """ A path of edges through a matrix """
    edges: list[Edge] = field(default_factory=list)
