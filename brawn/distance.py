""" Distance calculation methods """

from math import floor, log

# a precalculated table straight from MUSCLE covering regions
# of the domain where kimura would attempt a log of a negative
_KIMURA_TABLE = [
    1.95, 1.96, 1.97, 1.98, 1.99, 2.00, 2.00, 2.01, 2.02, 2.03,
    2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.09, 2.10, 2.11, 2.12,
    2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22,
    2.23, 2.24, 2.26, 2.27, 2.28, 2.29, 2.30, 2.31, 2.32, 2.33,
    2.34, 2.36, 2.37, 2.38, 2.39, 2.40, 2.41, 2.43, 2.44, 2.45,
    2.46, 2.48, 2.49, 2.50, 2.52, 2.53, 2.54, 2.55, 2.57, 2.58,
    2.60, 2.61, 2.62, 2.64, 2.65, 2.67, 2.68, 2.70, 2.71, 2.73,
    2.74, 2.76, 2.77, 2.79, 2.81, 2.82, 2.84, 2.85, 2.87, 2.89,
    2.91, 2.92, 2.94, 2.96, 2.98, 2.99, 3.01, 3.03, 3.05, 3.07,
    3.09, 3.11, 3.13, 3.15, 3.17, 3.19, 3.21, 3.23, 3.25, 3.28,
    3.30, 3.32, 3.35, 3.37, 3.39, 3.42, 3.44, 3.47, 3.49, 3.52,
    3.54, 3.57, 3.60, 3.62, 3.65, 3.68, 3.71, 3.74, 3.77, 3.80,
    3.83, 3.86, 3.89, 3.93, 3.96, 3.99, 4.03, 4.07, 4.10, 4.14,
    4.18, 4.22, 4.26, 4.30, 4.34, 4.38, 4.42, 4.47, 4.51, 4.56,
    4.61, 4.66, 4.71, 4.76, 4.82, 4.87, 4.93, 4.98, 5.04, 5.11,
    5.17, 5.24, 5.31, 5.38, 5.45, 5.53, 5.60, 5.69, 5.77, 5.86,
    5.95, 6.05, 6.15, 6.26, 6.37, 6.49, 6.61, 6.75, 6.88, 7.03,
    7.19, 7.36, 7.54, 7.75, 7.96, 8.19, 8.45, 8.74, 9.07, 9.45,
    9.88
]


# straight from MUSCLE
def calculate_distance(similarity: float) -> float:
    """ Calculates a distance from the given similarity, via different methods
        depending on where in the range 0-1 the similarity falls.
    """
    diff = 1 - similarity
    # use the full Kimura protein distance estimation if it wouldn't
    # cause a log domain error by being negative
    # (this should manage up to 0.85, but 0.75 matches MUSCLE and the values differ)
    if diff < 0.75:  # kimura proper (1983)
        return -log(1 - diff - diff**2 / 5)
    # again, match MUSCLE for anything above 0.93
    if diff > 0.93:
        return 10.
    # if it's in the awkward range, use MUSCLE's lookup table
    index = floor((diff - 0.75) * 1000 + .5)
    return _KIMURA_TABLE[index]
