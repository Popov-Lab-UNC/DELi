from deli.constants import MAKES, BB_CODES, MAKE_TO_LIBS


class DELMake:
    def __init__(self, make_id: str):
        self.make_id = make_id
        self._make_info: dict = MAKES[make_id]

        BB_CODES_NO_TOE = {key1: {key2: {key3[:-self._make_info["TOE_SIZE"]]: val3
                                         for key3, val3 in val2.items()}
                                  for key2, val2 in val1.items()}
                           for key1, val1 in BB_CODES.items() if key1 in MAKE_TO_LIBS[self.make_id]}

        BB_LENGTHS = {key: len(list(val["BBA"].keys())[0]) - self._make_info["TOE_SIZE"]
                      for key, val in BB_CODES.items() if key in MAKE_TO_LIBS[self.make_id]}

        self._make_info["BB_CODES_NO_TOE"] = BB_CODES_NO_TOE
        self._make_info["BB_LENGTHS"] = BB_LENGTHS

    def __getitem__(self, item):
        return self._make_info.get(item, None)
