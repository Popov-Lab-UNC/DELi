import os
import pickle
import sys
import json

from deli.match import find_matches

_pattern = sys.argv[1]

_num_no_matches = 0
_num_single_matches = 0
_num_multiple_matches = 0

matches = []
_count_file = sys.argv[3].replace(".pkl", ".count")
open(_count_file, "w").write(str(0) + "\n")
with open(sys.argv[2], 'r') as f:
    for i, line in enumerate(f):
        line = line.strip().split("\t")
        seq = line[1]
        _id = line[0]
        match = find_matches(seq, _pattern, _id)

        if len(match) == 0:
            _num_no_matches += 1
        elif len(match) == 1:
            _num_single_matches += 1
        else:
            _num_multiple_matches += 1
        matches.extend(match)

        if ((i + 1) % 10000) == 0:
            open(_count_file, "w").write(str(i+1) + "\n")

data = {
    '_num_no_matches': _num_no_matches,
    '_num_single_matches': _num_single_matches,
    '_num_multiple_matches': _num_multiple_matches,
}

pickle.dump(matches, open(sys.argv[3], "wb"))
json.dump(data, open(sys.argv[3].replace(".pkl", ".json"), "w"))
open(_count_file, "w").write(str(i+1) + "\n")
