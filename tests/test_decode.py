import time

from deli.dels import DELibrary, DELibraryPool
from deli.decode import DecodingExperimentRunner, DecodingExperiment, DecodingSettings
from deli.dna.io import SingleSequenceReader

lib4 = DELibrary.load("DEL004")
lib5 = DELibrary.load("DEL005")

pool = DELibraryPool([lib4, lib5])
reader = SingleSequenceReader("TEST_SEQS")

decode_settings = DecodingSettings(revcomp=True, bb_calling_approach="bio")
experiment = DecodingExperiment(pool, target_id="TEST", decode_settings=decode_settings, experiment_id="TEST")
runner = DecodingExperimentRunner(experiment, debug=False, disable_logging=True)

t0 = time.time()
runner.run(reader, use_tqdm=True)
print(time.time() - t0)
