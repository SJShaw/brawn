# pylint: disable=missing-function-docstring,missing-module-docstring,missing-class-docstring
# pylint: disable=protected-access

from io import StringIO

from brawn import Alignment

from .shared import SMALL_DB, SMALL_QUERY, dict_to_fasta, near_equal


def test_from_local_dict():
    # an alternate constructor, but the rest of the tests are fine to reuse
    test_from_local_trivial(Alignment(SMALL_QUERY))


def test_from_local_trivial(alignment=None):
    if not alignment:
        alignment = Alignment.from_file(StringIO(dict_to_fasta(SMALL_QUERY)))
    profile_positions = alignment.positions

    assert len(profile_positions) == 4
    for prof in profile_positions:
        assert sum(prof.base_counts) == 1
        assert prof.ungapped_weight == 1
        assert prof.gap_opens == 1
        assert prof.gap_closes == 1
        assert prof.score_gap_open == -1.45
        assert prof.score_gap_close == -1.45
        assert prof.base_counts[prof.sort_order[0]] == 1
        assert prof.base_counts[prof.sort_order[1]] == 0

    position = profile_positions[0]
    assert position.sort_order[:2] == [5, 0]
    assert near_equal(position.scores, [
        1.09860003, 0.655409992, 0.909370005, 0.771830022, 0.301400006, 5.62829018,
        0.641910017, 0.284319997, 0.678740025, 0.305489987, 0.377389997, 1.01012003,
        0.608510017, 0.659959972, 0.636600018, 1.03447998, 0.684350014, 0.407279998,
        0.360339999, 0.356790006
    ])

    position = profile_positions[1]
    assert position.sort_order[:2] == [16, 0]
    assert near_equal(position.scores, [
        1.20372999, 0.990689993, 0.881910026, 0.88132, 0.571720004, 0.684350014,
        0.835430026, 0.881680012, 0.927299976, 0.699240029, 0.877590001, 1.11251998,
        0.830139995, 0.931029975, 0.806729972, 1.52910995, 2.58221006, 0.987020016,
        0.315409988, 0.579540014
    ])

    position = profile_positions[2]
    assert position.sort_order[:2] == [7, 0]
    assert near_equal(position.scores, [
        0.808049977, 0.971849978, 0.328599989, 0.431510001, 1.10968995, 0.284319997,
        0.507830024, 3.03765988, 0.493099988, 1.88885999, 1.75039005, 0.442460001,
        0.444310009, 0.532130003, 0.481530011, 0.556029975, 0.881680012, 2.3736701,
        0.684939981, 0.700349987
    ])

    position = profile_positions[3]
    assert position.sort_order[:2] == [17, 0]
    assert near_equal(position.scores, [
        1.05955994, 1.21604002, 0.435570002, 0.540470004, 0.912559986, 0.407279998,
        0.548169971, 2.3736701, 0.554669976, 1.50372005, 1.42742002, 0.506030023,
        0.56795001, 0.610849977, 0.514219999, 0.677670002, 0.987020016, 2.6558001,
        0.434190005, 0.63805002
    ])


def test_profiles_less_trivial():
    msa = Alignment.from_file(StringIO(dict_to_fasta(SMALL_DB)))
    profile_positions = msa.positions

    assert len(profile_positions) == 6
    for prof in profile_positions:
        assert sum(prof.base_counts) == 1
        assert prof.base_counts[prof.sort_order[0]] == 1
        assert prof.base_counts[prof.sort_order[1]] == 0

    assert [p.gap_opens for p in profile_positions] == [1, 1, .5, .5, 1, 1]
    assert [p.gap_closes for p in profile_positions] == [1, 1, .5, .5, 1, 1]
    assert [p.ungapped_weight for p in profile_positions] == [1, 1, .5, .5, 1, 1]

    position = profile_positions[0]
    assert position.sort_order[:2] == [5, 0]
    assert near_equal(position.scores, [
        1.09860003, 0.655409992, 0.909370005, 0.771830022, 0.301400006, 5.62829018,
        0.641910017, 0.284319997, 0.678740025, 0.305489987, 0.377389997, 1.01012003,
        0.608510017, 0.659959972, 0.636600018, 1.03447998, 0.684350014, 0.407279998,
        0.360339999, 0.356790006
    ])
    assert position.score_gap_open == -1.45
    assert position.score_gap_close == -1.45

    position = profile_positions[1]
    assert position.sort_order[:2] == [16, 0]
    assert sum(position.base_counts) == 1
    assert near_equal(position.scores, [
        1.20372999, 0.990689993, 0.881910026, 0.88132, 0.571720004, 0.684350014,
        0.835430026, 0.881680012, 0.927299976, 0.699240029, 0.877590001, 1.11251998,
        0.830139995, 0.931029975, 0.806729972, 1.52910995, 2.58221006, 0.987020016,
        0.315409988, 0.579540014
    ])
    assert position.score_gap_open == -1.45
    assert position.score_gap_close == -1.45

    position = profile_positions[2]
    assert position.sort_order[:2] == [8, 0]
    assert near_equal(position.scores, [
        0.812129974, 0.464139998, 1.03391004, 1.35988998, 0.370689988, 0.678740025,
        1.03822005, 0.493099988, 2.7288301, 0.527390003, 0.682439983, 1.15671003,
        0.829110026, 1.51332998, 2.33521008, 0.938579977, 0.927299976, 0.554669976,
        0.399439991, 0.525489986
    ])
    assert position.score_gap_open == -0.725
    assert position.score_gap_close == -0.725

    position = profile_positions[3]
    assert position.sort_order[:2] == [2, 0]
    assert near_equal(position.scores, [
        0.827040017, 0.398620009, 4.18833017, 2.06850004, 0.251940012, 0.909370005,
        1.01617002, 0.328599989, 1.03391004, 0.312999994, 0.424980015, 1.80887997,
        0.813069999, 1.20043004, 0.637120008, 1.03000998, 0.881910026, 0.435570002,
        0.263130009, 0.379469991
    ])
    assert position.score_gap_open == -0.725
    assert position.score_gap_close == -0.725

    position = profile_positions[4]
    assert position.sort_order[:2] == [17, 0]
    assert near_equal(position.scores, [
        1.05955994, 1.21604002, 0.435570002, 0.540470004, 0.912559986, 0.407279998,
        0.548169971, 2.3736701, 0.554669976, 1.50372005, 1.42742002, 0.506030023,
        0.56795001, 0.610849977, 0.514219999, 0.677670002, 0.987020016, 2.6558001,
        0.434190005, 0.63805002
    ])
    assert position.score_gap_open == -1.45
    assert position.score_gap_close == -1.45

    position = profile_positions[5]
    assert position.sort_order[:2] == [5, 0]
    assert near_equal(position.scores, [
        1.09860003, 0.655409992, 0.909370005, 0.771830022, 0.301400006, 5.62829018,
        0.641910017, 0.284319997, 0.678740025, 0.305489987, 0.377389997, 1.01012003,
        0.608510017, 0.659959972, 0.636600018, 1.03447998, 0.684350014, 0.407279998,
        0.360339999, 0.356790006
    ])
    assert position.score_gap_open == -1.45
    assert position.score_gap_close == -1.45
