import json
import pickle
import sys

from deli.match import SequenceMatch
from deli.called import call_hit
from deli.del_make import DELMake
from deli.utils import get_min_index_distance, get_min_lib_distance
from deli.constants import FILE_HEADER

data = pickle.load(open(sys.argv[1], "rb"))

libs = sys.argv[2].split(",")
idxs = sys.argv[3].split(",")
strict = bool(int(sys.argv[4]))

outfile = open(sys.argv[5], "w")
outfile.write(FILE_HEADER + "\n")

make = DELMake(sys.argv[6])

_num_called = 0
_num_error = 0
_num_failed = 0
_num_passed = 0

_min_lib_dist = get_min_lib_distance(libs)
_min_idx_dist = get_min_index_distance(idxs)

_count_file = sys.argv[5].replace(".csv", ".count")
open(_count_file, "w").write(str(0) + "\n")

try:
    data_number = int(sys.argv[1].split("_")[-1].split(".")[0])
except ValueError:
    data_number = 1

for i, hit in enumerate(data):
    _hit = SequenceMatch(hit[0], hit[1], hit[4], hit[5], hit[2][0], hit[3][0],
                         hit[2][1], hit[3][1], hit[2][2], hit[3][2], data_number*i)
    call = call_hit(_hit, libs, idxs, make, strict, _min_lib_dist, _min_idx_dist)

    if call is None:
        _num_error += 0
        continue
    else:
        _num_called += 1
        if call.del_id == "FAILED":
            _num_failed += 1
        else:
            _num_passed += 1
        outfile.write(call.to_row())

    if ((i + 1) % 1000) == 0:
        open(_count_file, "w").write(str(i + 1) + "\n")

data = {
    "num_called": _num_called,
    "num_error": _num_error,
    "num_failed": _num_failed,
    "num_passed": _num_passed
}

outfile.close()
json.dump(data, open(sys.argv[5].replace(".csv", ".json"), "w"))
open(_count_file, "w").write(str(_num_called) + "\n")
