from dataclasses import dataclass
from collections import Counter

from Levenshtein import distance

from .constants import UMI_CLUSTER_THRESHOLD


@dataclass
class UMICluster:
    umi_set: set[str]
    umi_cluster_id: int

    def includes_umi(self, umi: str) -> bool:
        min_dist = 100
        for _cluster_umi in self.umi_set:
            _dist = distance(umi, _cluster_umi)
            if _dist < min_dist:
                min_dist = _dist
            if min_dist <= UMI_CLUSTER_THRESHOLD:  # early stopping
                return True
        return min_dist <= UMI_CLUSTER_THRESHOLD

    def add_umi(self, umi: str):
        self.umi_set.add(umi)

    def __contains__(self, item):
        if isinstance(item, str):
            return self.includes_umi(item)
        return False


def cluster_umis(umis: list[str], return_count: bool = True):
    umi_set = list(dict(sorted(Counter(umis).items(), key=lambda x: x[1])).keys())
    _clusters = []
    _cluster_ids = []
    for umi in umi_set:
        _found_cluster = False
        for _cluster in _clusters:
            if _cluster.includes_umi(umi):
                _found_cluster = True
                break
        if not _found_cluster:
            _clusters.append(UMICluster({umi}, len(_clusters)))

    for umi in umis:
        found_clusters = []
        for _cluster in _clusters:
            if _cluster.includes_umi(umi):
                found_clusters.append(_cluster.umi_cluster_id)
        _cluster_ids.append(tuple(found_clusters))

    if return_count:
        return max([__ for _ in _cluster_ids for __ in _])
    else:
        return _cluster_ids
