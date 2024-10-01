import json

from deli.constants import MAX_INDEX_RISK_DIST_THRESHOLD
from Levenshtein import distance


class Index:
    def __init__(
            self,
            index_id: str,
            dna_tag: str
    ):
        self.index_id = index_id
        self.dna_tag = dna_tag

    @classmethod
    def from_json(cls, path: str):
        data = json.load(open(path))
        return cls(**data)


def get_min_index_distance(included_index: list[Index]) -> int:
    """
    Given a list of index IDs, determine the minimum Levenshtein distance between all of them

    Parameters
    ----------
    included_index: list[Index]
        the index_ids to use

    Returns
    -------
    min_distance: int
        the minimum Levenshtein distance between the passed indexes
    """
    # if only one index just give it a constant threshold
    if len(included_index) == 1:
        return MAX_INDEX_RISK_DIST_THRESHOLD
    if included_index is None:
        return 0
    _index_sequences = [i.dna_tag for i in included_index]
    return min([distance(s1, s2) for s1 in _index_sequences for s2 in _index_sequences if s1 != s2])
