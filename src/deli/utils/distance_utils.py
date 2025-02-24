"""utilities for calculating distance"""

from typing import Union

from Levenshtein import distance as levenshtein_distance


def get_min_levenshtein_distance(items: list[str]) -> int:
    """
    get the minimum Levenshtein distance between any pair of items

    Parameters
    ----------
    items: list[str]
        items to check distance for

    Returns
    -------
    min_distance: int
    """
    min_distance: Union[int, float] = float("inf")
    for i, item1 in enumerate(items):
        for item2 in items[i + 1 :]:
            min_distance = min(min_distance, levenshtein_distance(item1, item2))
    return int(min_distance)
