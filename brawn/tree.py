""" Tree related code """

from dataclasses import dataclass

from .constants import ID_GUARD


@dataclass
class Tree:  # pylint: disable=too-many-instance-attributes
    """ A binary tree for sequences within an alignment """
    node_count: int
    parents: list[int]
    lefts: list[int]
    rights: list[int]
    left_lengths: list[float]
    right_lengths: list[float]
    parent_lengths: list[float]
    names: list[str]
    root_node_index: int
    leaf_count: int = -1

    def __post_init__(self) -> None:
        self.leaf_count = (self.node_count + 1) // 2
        assert self.node_count > 0, self.node_count
        assert len(self.parents) == self.node_count
        assert len(self.lefts) == self.node_count
        assert len(self.rights) == self.node_count
        assert len(self.left_lengths) == self.node_count
        assert len(self.right_lengths) == self.node_count
        assert len(self.names) == self.leaf_count
        assert self.parents[self.root_node_index] == ID_GUARD

    def get_parent(self, index: int) -> int:
        """ Returns the parent of the given node index """
        parent = self.parents[index]
        assert 0 <= parent < self.node_count
        return parent

    def get_left(self, index: int) -> int:
        """ Returns the left child of the given node index """
        return self.lefts[index]

    def get_right(self, index: int) -> int:
        """ Returns the right child of the given node index """
        return self.rights[index]

    def is_leaf(self, index: int) -> bool:
        """ Returns True if it's a leaf """
        assert index < self.node_count, f"{index} > {self.node_count}"
        return (self.node_count == 1
                or self.lefts[index] == ID_GUARD and self.rights[index] == ID_GUARD)

    def get_edge_length(self, first: int, second: int) -> float:
        """ Returns the length between two neighbouring nodes
        """
        if self.lefts[first] == second:
            return self.left_lengths[first]
        if self.rights[first] == second:
            return self.right_lengths[first]
        if self.parents[first] != second:
            raise ValueError(f"Nodes are not neighbours: {first}, {second}")
        return self.parent_lengths[first]

    @property
    def node_child_counts(self) -> list[int]:
        """ The counts of leaves under a nodes in the tree (including the node itself)

            Returns:
                a list of counts, each one the number of leaves under the relevant node
        """
        leaves_under_node = [0] * self.node_count

        def find_count(index: int) -> int:
            if len(leaves_under_node) == 1:
                return 1
            if self.is_leaf(index):
                leaves_under_node[index] = 1
                return 1

            left = self.get_left(index)
            left_count = find_count(left)

            right = self.get_right(index)
            right_count = find_count(right)
            count = left_count + right_count
            leaves_under_node[index] = count
            return count

        find_count(self.root_node_index)

        return leaves_under_node

    def get_weights(self) -> list[float]:
        """ Calculates sequence weights from a tree of the sequences

            Arguments:
                tree: the tree of sequences

            Returns:
                a normalised list of weights, summing to 1
        """
        weights: list[float] = []

        leaf_count = self.leaf_count
        if not leaf_count:
            return weights

        if leaf_count == 1:
            return [1.0]

        if leaf_count == 2:
            return [0.5, 0.5]

        strengths = []
        for i, leaves in enumerate(self.node_child_counts):
            if self.root_node_index == i:
                strengths.append(0.)
                continue
            parent = self.get_parent(i)
            length = self.get_edge_length(i, parent)
            strengths.append(length / leaves)

        for node in range(leaf_count):
            weight = 0.
            while node != self.root_node_index:
                weight += strengths[node]
                node = self.get_parent(node)
            if weight < 0.0001:
                weight = 1.
            weights.append(weight)
        return _normalise(weights)


def _normalise(values: list[float]) -> list[float]:
    """ Constructs a normalised list of values from the input,
        with the sum of the result being 1
    """
    total = sum(values)
    return [v / total for v in values]
