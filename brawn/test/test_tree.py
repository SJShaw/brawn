# pylint: disable=missing-function-docstring,missing-module-docstring,missing-class-docstring
# pylint: disable=protected-access

import unittest

from brawn.alignment import Alignment, tree_from_alignment
from brawn.constants import ID_GUARD, LENGTH_GUARD
from brawn.tree import Tree, _normalise as normalise

from .shared import MEDIUM_DB_PATH, near_equal


class TestTree(unittest.TestCase):
    def setUp(self):
        self.lefts = ([ID_GUARD] * 4) + [0, 3, 1]
        self.rights = ([ID_GUARD] * 4) + [2, 4, 5]
        self.parents = [4, 6, 4, 5, 5, 6, ID_GUARD]
        self.tree = Tree(
                len(self.lefts), lefts=self.lefts, rights=self.rights, parents=self.parents,
                root_node_index=6,
                # arbitrarily set lengths with unique values for directionality tests
                right_lengths=[1, 3, 5, 7, 9, 11, 13],
                left_lengths=[2, 4, 6, 8, 10, 12, 14],
                parent_lengths=[50, 60, 70, 80, 90, 100, 110],
                names=list("ABCD"),
        )

    def test_simple(self):
        tree = self.tree
        # the values of the tree better be the values given
        assert tree.leaf_count == 4
        assert tree.lefts == self.lefts
        assert tree.rights == self.rights
        assert tree.parents == self.parents
        assert tree.root_node_index == 6

        # test leaves register correctly as leaves and the node ids map correctly
        is_leaf = [tree.is_leaf(i) for i in range(tree.node_count)]
        assert not any(is_leaf[-(tree.leaf_count - 1):])
        assert all(is_leaf[:tree.leaf_count - 1])
        leaves = [i for i in range(tree.node_count) if tree.is_leaf(i)]
        assert leaves == [0, 1, 2, 3]
        assert self.tree.node_child_counts == [1, 1, 1, 1, 2, 3, 4]

    def test_edge_lengths(self):
        assert self.lefts[5] == 3
        assert self.parents[3] == 5
        assert self.tree.get_edge_length(5, 3) == 12
        # this ought to be the same, but it's forced to be different above
        assert self.tree.get_edge_length(3, 5) == 80

        assert self.rights[4] == 2
        assert self.parents[2] == 4
        assert self.tree.get_edge_length(4, 2) == 9
        # same again, forced to be different
        assert self.tree.get_edge_length(2, 4) == 70

        # and then test non-connected chunks
        with self.assertRaisesRegex(ValueError, "not neighbours"):
            self.tree.get_edge_length(0, 6)


def test_tree_full():
    with open(MEDIUM_DB_PATH, encoding="utf-8") as handle:
        alignment = Alignment.from_file(handle)
    tree = tree_from_alignment(alignment)
    # check that the lengths have been calculated correctly
    assert near_equal(tree.parent_lengths,
                      [0.8389, 1.33617, 0.8389, 1.1715, 0.3326, 0.16467, LENGTH_GUARD])
    assert near_equal(tree.left_lengths,
                      ([LENGTH_GUARD] * 4) + [0.8389, 1.1715, 1.33617])
    assert near_equal(tree.right_lengths,
                      ([LENGTH_GUARD] * 4) + [0.8389, 0.3326, 0.16467])


def test_normalise():
    assert normalise([3, 7, 0]) == [0.3, 0.7, 0.]
    assert normalise([2, 3]) == [0.4, 0.6]
