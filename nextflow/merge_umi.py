"""nextflow process script for merging decode runs into single output when using UMI degen"""

import argparse
import gzip
import json
import os
from collections import defaultdict
from itertools import zip_longest

from deli.decode import DecodeStatistics, build_decoding_report
from deli.dels import Selection


parser = argparse.ArgumentParser()
parser.add_argument("--decode_file")
parser.add_argument("--include_bb_smi", action="store_true")
parser.add_argument("--enumerate_smiles", action="store_true")
parser.add_argument(
    "--save_counter",
    action="store_true",
)
parser.add_argument("--counters", nargs="+")
args = parser.parse_args()

# merge the counter files into a single count
counter_files = args.counters
print(counter_files)
data = json.load(gzip.open(counter_files[0], "rt", encoding="utf-8"))
for f in counter_files[1:]:
    new_data = json.load(gzip.open(f, "rt", encoding="utf-8"))
    for lib_id, counter_data in new_data.items():
        if lib_id not in data.keys():
            data[lib_id] = {}
        for del_id, count_data in counter_data.items():
            if del_id not in data[lib_id]:
                data[lib_id][del_id] = {
                    "lib_id": count_data["lib_id"],
                    "bb_ids": count_data["bb_ids"],
                    "raw_count": 0,
                    "umis": [],
                }
            data[lib_id][del_id]["raw_count"] += count_data["raw_count"]
            data[lib_id][del_id]["umis"].extend(count_data["umis"])

print(len(data["DEL004"]), len(data["DEL005"]))

# merge the relevant statistics from the stats files
statistic_files = [f for f in os.listdir() if f.endswith("_decode_statistics_subjob.json")]
merged_stats = DecodeStatistics.from_file(statistic_files[0])
for statistic_file in statistic_files[1:]:
    loaded_statistics = DecodeStatistics.from_file(statistic_file)
    merged_stats += loaded_statistics

# build the cube file
selection = Selection.from_yaml(args.decode_file)
_max_cycle_size = selection.library_collection.max_cycle_size()

# build cube file header
_header = "DEL_ID"
_header += ",LIBRARY_ID"
for i in range(_max_cycle_size):
    _header += f",BB{i + 1}_ID"
    if args.include_bb_smi:
        _header += f",BB{i + 1}_SMILES"
if args.enumerate_smiles:
    _header += ",SMILES"
_header += ",RAW_COUNT"
_header += ",UMI_CORRECTED_COUNT\n"

updated_num_seqs_degen_per_lib = defaultdict(
    int, {lib.library_id: 0 for lib in selection.library_collection.libraries}
)

# write the new cube_file
with open(f"./{selection.selection_id}_cube.csv", "w") as f:
    f.write(_header)
    for lib_id, lib_data in data.items():
        for del_id, compound_data in lib_data.items():
            degen_count = len(set(compound_data["umis"]))
            updated_num_seqs_degen_per_lib[lib_id] += degen_count
            _row = f"{del_id},{lib_id}"
            for _, bb in zip_longest(range(_max_cycle_size), compound_data["bb_ids"]):
                _row += f",{bb if bb is not None else 'null'}"
                if args.include_bb_smi:
                    _row += f",{bb.smiles if (bb is not None and bb.smiles) else 'null'}"
                if args.enumerate_smiles:
                    try:
                        _smi = (
                            selection.library_collection.get_library(lib_id)
                            .enumerate_by_bb_ids(compound_data["bb_ids"])
                            .smi
                        )
                    except Exception:
                        _smi = "null"
                    _row += f",{_smi}"
            _row += f",{compound_data['raw_count']}"
            _row += f",{degen_count}\n"

# reset the umi counts in the stats object
merged_stats.num_seqs_degen_per_lib = updated_num_seqs_degen_per_lib

# save counter if requested
if args.save_counter:
    with gzip.open(f"./{selection.selection_id}_counter.json", "wt", encoding="utf-8") as f:
        json.dump(data, f)

# write the report at the end
build_decoding_report(
    selection=selection,
    stats=merged_stats,
    out_path=f"./{selection.selection_id}_decode_report.html",
)

# write the statistics file
statistics_out_path = f"./{selection.selection_id}_decode_statistics.json"
merged_stats.to_file(statistics_out_path)
