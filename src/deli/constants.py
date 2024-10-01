# TODO this should load these constants from a configuration file?

import os


BB_MASK = "###"
BB_NULL = "NUL"

DELI_DATA_DIR = os.getenv("DELI_DATA_DIR")

MAX_INDEX_RISK_DIST_THRESHOLD = 4
MAX_LIBRARY_RISK_DIST_THRESHOLD = 3

# PACKAGE_DIR = os.path.dirname(os.path.realpath(__file__))
# MODULE_DIR = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
#
# FILE_HEADER = ("MATCH_ID,READ_ID,OVERALL_PASS,HAMMING_FIXABLE,TOTAL_RISK,MATCH_RISK,CALLED_LIB,LIB_RISK,CALLED_INDEX,"
#                "INDEX_RISK,DEL_ID,BBA_ID,BBA_RISK,BBB_ID,BBB_RISK,BBC_ID,BBC_RISK,UMI,MATCH_SEQ")
#
# # risk constants are exclusive threshold (> or <)
# MATCH_RISK_THRESHOLD = 3
# INDEX_RISK_NORMAL_THRESHOLD = 3
# LIBRARY_RISK_NORMAL_THRESHOLD = 3
# UMI_CLUSTER_THRESHOLD = 2  # inclusive, so if the distance is 2 then it will be in the cluster
#
# BB_CODES = {}
# for _dir in os.listdir(os.path.join(os.path.dirname(__file__), '../data/bb_codes')):
#     BB_CODES[_dir] = {}
#     for _bb_file in os.listdir(os.path.join(os.path.dirname(__file__), '../data/bb_codes', _dir)):
#         _bb_cycle = _bb_file.split("_")[-1].split(".")[0]
#         _path = os.path.join(os.path.dirname(__file__), '../data/bb_codes', _dir, _bb_file)
#         BB_CODES[_dir][_bb_cycle] = json.load(open(_path, "r"))
#
# INDEX = json.load(open(os.path.join(os.path.dirname(__file__), '../data/experiment_index.json'), 'r'))
# LIBS = json.load(open(os.path.join(os.path.dirname(__file__), '../data/libs.json'), 'r'))
# MAKES = json.load(open(os.path.join(os.path.dirname(__file__), '../data/barcode_makes.json'), 'r'))
# LIBS_TO_MAKE = json.load(open(os.path.join(os.path.dirname(__file__), '../data/libs_to_make.json'), 'r'))
# MAKE_TO_LIBS = {}
# for lib, make in LIBS_TO_MAKE.items():
#     MAKE_TO_LIBS.setdefault(make, []).append(lib)
#
#
#
# class DeliConfigError(Exception):
#     """
#     Exception class to be thrown whenever there is an issue with data missing from the DELi configs
#     """
#     pass
