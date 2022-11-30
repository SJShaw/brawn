# pylint: disable=missing-function-docstring,missing-module-docstring,missing-class-docstring
# pylint: disable=protected-access

import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
SMALL_QUERY = {
    "query": "GTIV"
}
SMALL_DB = {
    "A": "GT-DVG",
    "B": "GTK-VG",
}
MEDIUM_QUERY_PATH = os.path.join(DATA_DIR, "query.fasta")
MEDIUM_DB_PATH = os.path.join(DATA_DIR, "db.fasta")
MEDIUM_OUTPUT_PATH = os.path.join(DATA_DIR, "out.fasta")


def near_equal(list1, list2):
    if len(list1) != len(list2):
        return False
    for i, j in zip(list1, list2):
        if not -1e-5 < i - j < 1e-5:
            print(i, "!=", j)
            return False
    return True


# just in case, since it's no good if it doesn't work as intended
assert not near_equal([1], [1, 1])
assert not near_equal([1, 2], [1 + 1e-4, 2])
assert near_equal([1 - 1e-6, 2 + 1e-6], [1, 2])


def dict_to_fasta(data):
    lines = []
    for name, seq in data.items():
        lines.append(f">{name}")
        lines.append(f"{seq}")
    return "\n".join(lines)


# again, check it just in case
assert dict_to_fasta({"A": "AB", "C": "CD"}) == ">A\nAB\n>C\nCD"
