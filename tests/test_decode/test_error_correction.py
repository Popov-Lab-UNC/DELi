"""tests for building block error correction classes/functions"""

import pytest

from deli.decode.barcode_calling import HammingDistBarcodeCaller, HashMapCollisionError, LevenshteinDistBarcodeCaller
from deli.dels.building_block import TaggedBuildingBlock, TaggedBuildingBlockSet


class TestHammingHashCorrectors:
    """test the accuracy of the HammingDistBarcodeCaller correctors"""

    @pytest.fixture
    def tagged_hamming3_bb_set(self) -> TaggedBuildingBlockSet:
        """A TaggedBuildingBlockSet fixture for testing"""
        return TaggedBuildingBlockSet(
            bb_set_id="TEST_BB_SET",
            building_blocks=[
                TaggedBuildingBlock(tags="TGGGCGTT", bb_id="BB1"),
                TaggedBuildingBlock(tags="GCTTGCCA", bb_id="BB2"),
                TaggedBuildingBlock(tags="TCCCAGAC", bb_id="BB3"),
                TaggedBuildingBlock(tags="GTGTTCCC", bb_id="BB4"),
                TaggedBuildingBlock(tags="AGTCATAG", bb_id="BB5"),
                TaggedBuildingBlock(tags="TAGGGAGG", bb_id="BB6"),
                TaggedBuildingBlock(tags="CTCGTCCC", bb_id="BB7"),
                TaggedBuildingBlock(tags="CGTACCAC", bb_id="BB8"),
                TaggedBuildingBlock(tags="GGGGCAAG", bb_id="BB9"),
                TaggedBuildingBlock(tags="GGCCTGTT", bb_id="BB10"),
                TaggedBuildingBlock(tags="AGCGGCGA", bb_id="BB11"),
                TaggedBuildingBlock(tags="TCGTAAGA", bb_id="BB12"),
                TaggedBuildingBlock(tags="ACTCGCTT", bb_id="BB13"),
                TaggedBuildingBlock(tags="TTTAGCTA", bb_id="BB14"),
                TaggedBuildingBlock(tags="CGAGACGG", bb_id="BB15"),
                TaggedBuildingBlock(tags="GCTCCCTA", bb_id="BB16"),
                TaggedBuildingBlock(tags="TGCACAAA", bb_id="BB17"),
                TaggedBuildingBlock(tags="CTTGTTGT", bb_id="BB18"),
                TaggedBuildingBlock(tags="CTGCTGCT", bb_id="BB19"),
                TaggedBuildingBlock(tags="TCCAGGTC", bb_id="BB20"),
            ],
        )

    @pytest.fixture
    def tagged_hamming5_bb_set(self) -> TaggedBuildingBlockSet:
        """A TaggedBuildingBlockSet fixture for testing"""
        return TaggedBuildingBlockSet(
            bb_set_id="TEST_BB_SET",
            building_blocks=[
                TaggedBuildingBlock(tags="ACCTTCAACCGT", bb_id="BB1"),
                TaggedBuildingBlock(tags="TGTGTCACGCTG", bb_id="BB2"),
                TaggedBuildingBlock(tags="GACTTTGGCTTT", bb_id="BB3"),
                TaggedBuildingBlock(tags="GGCTAATCATAG", bb_id="BB4"),
                TaggedBuildingBlock(tags="ATTTCGGGAAAC", bb_id="BB5"),
                TaggedBuildingBlock(tags="GCAAAGCCCCGT", bb_id="BB6"),
                TaggedBuildingBlock(tags="GTCAGGCTGCCA", bb_id="BB7"),
                TaggedBuildingBlock(tags="AACTGACATTGA", bb_id="BB8"),
                TaggedBuildingBlock(tags="CACGATCTTGGA", bb_id="BB9"),
                TaggedBuildingBlock(tags="CCTCCTTGATTG", bb_id="BB10"),
                TaggedBuildingBlock(tags="TCGTTCGGGGGG", bb_id="BB11"),
                TaggedBuildingBlock(tags="TCCTGCAGGCTC", bb_id="BB12"),
                TaggedBuildingBlock(tags="CTGACTGCCAAG", bb_id="BB13"),
                TaggedBuildingBlock(tags="AGGCTAGAGAAC", bb_id="BB14"),
                TaggedBuildingBlock(tags="TCCTGTTGGTAA", bb_id="BB15"),
                TaggedBuildingBlock(tags="GAAGAGTCACTT", bb_id="BB16"),
                TaggedBuildingBlock(tags="CAGTAGTCGTAT", bb_id="BB17"),
                TaggedBuildingBlock(tags="TCGAACAGCTAC", bb_id="BB18"),
                TaggedBuildingBlock(tags="CCCGACGTCTAA", bb_id="BB19"),
                TaggedBuildingBlock(tags="AGATATTTGACG", bb_id="BB20"),
            ],
        )

    def test_hamming3_1_hashmap_corrector(self, tagged_hamming3_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 3 set with non-asym dist cutoff 3"""
        # set up the hamming 3 correctors
        hamming3_1_corrector = HammingDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_hamming3_bb_set.building_blocks for tag in bb.tags}, distance_cutoff=1
        )

        # test the correctors
        bb = tagged_hamming3_bb_set.building_blocks[0]
        res = hamming3_1_corrector._hash_map["TGGGCGTT"]  # direct lookup should pass
        assert res == (bb, 0)
        res = hamming3_1_corrector._hash_map["TGGGCGTA"]  # 1 error should pass
        assert res == (bb, 1)
        with pytest.raises(KeyError):
            hamming3_1_corrector._hash_map["TGGGCGAA"]  # 2 errors should fail
        with pytest.raises(KeyError):
            hamming3_1_corrector._hash_map["AAAAAAAA"]  # not in set should fail

    def test_hamming3_2_hashmap_corrector(self, tagged_hamming3_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 3 set with asym dist cutoff 2"""
        # set up the hamming 3 correctors
        with pytest.raises(HashMapCollisionError):
            HammingDistBarcodeCaller(
                tag_map={tag: bb for bb in tagged_hamming3_bb_set.building_blocks for tag in bb.tags},
                distance_cutoff=2,
                asymmetrical=False,
            )

        hamming3_2_corrector = HammingDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_hamming3_bb_set.building_blocks for tag in bb.tags},
            distance_cutoff=2,
            asymmetrical=True,
        )

        # test the correctors
        bb = tagged_hamming3_bb_set.building_blocks[0]
        res = hamming3_2_corrector._hash_map["TGGGCGTT"]  # direct lookup should pass
        assert res == (bb, 0)
        res = hamming3_2_corrector._hash_map["TGGGCGTA"]  # 1 error should pass
        assert res == (bb, 1)
        res = hamming3_2_corrector._hash_map["TGGGCGAA"]  # 2 errors should pass
        assert res == (bb, 2)
        with pytest.raises(KeyError):
            hamming3_2_corrector._hash_map["AAAAAAAA"]  # not in set should fail

    def test_hamming5_1_hashmap_corrector(self, tagged_hamming5_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 5 set with non-asym dist cutoff 1"""
        # set up the hamming 5 correctors
        hamming5_1_corrector = HammingDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_hamming5_bb_set.building_blocks for tag in bb.tags}, distance_cutoff=1
        )

        # test the correctors
        bb = tagged_hamming5_bb_set.building_blocks[0]
        res = hamming5_1_corrector._hash_map["ACCTTCAACCGT"]  # direct lookup should pass
        assert res == (bb, 0)
        res = hamming5_1_corrector._hash_map["ACCTTCAACCCT"]  # 1 error should pass
        assert res == (bb, 1)
        with pytest.raises(KeyError):
            hamming5_1_corrector._hash_map["ACCTTCAAGGGT"]  # 2 errors should fail
        with pytest.raises(KeyError):
            hamming5_1_corrector._hash_map["AAAAAAAAAAA"]  # not in set should fail

    def test_hamming5_2_hashmap_corrector(self, tagged_hamming5_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 5 set with non-asym dist cutoff 2"""
        # set up the hamming 5 correctors
        with pytest.raises(HashMapCollisionError):
            HammingDistBarcodeCaller(
                tag_map={tag: bb for bb in tagged_hamming5_bb_set.building_blocks for tag in bb.tags},
                distance_cutoff=3,
                asymmetrical=False,
            )

        hamming5_2_corrector = HammingDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_hamming5_bb_set.building_blocks for tag in bb.tags},
            distance_cutoff=2,
            asymmetrical=True,
        )

        # test the correctors
        bb = tagged_hamming5_bb_set.building_blocks[0]
        res = hamming5_2_corrector._hash_map["ACCTTCAACCGT"]  # direct lookup should pass
        assert res == (bb, 0)
        res = hamming5_2_corrector._hash_map["ACCTTCAACCCT"]  # 1 error should pass
        assert res == (bb, 1)
        res = hamming5_2_corrector._hash_map["ACCTTCAAGGGT"]  # 2 errors should pass
        assert res == (bb, 2)
        with pytest.raises(KeyError):
            hamming5_2_corrector._hash_map["AAAAAAAAAAA"]  # not in set should fail


class TestLevenshteinHashCorrectors:
    """Test the accuracy of the LevenshteinDistBarcodeCaller correctors"""

    @pytest.fixture
    def tagged_levenshtein3_bb_set(self) -> TaggedBuildingBlockSet:
        """A TaggedBuildingBlockSet fixture for testing"""
        return TaggedBuildingBlockSet(
            bb_set_id="TEST_BB_SET",
            building_blocks=[
                TaggedBuildingBlock(tags="TCGGTGGA", bb_id="BB1"),
                TaggedBuildingBlock(tags="AGCCTCAC", bb_id="BB2"),
                TaggedBuildingBlock(tags="CTTTAATA", bb_id="BB3"),
                TaggedBuildingBlock(tags="CAATAAAT", bb_id="BB4"),
                TaggedBuildingBlock(tags="TTTGCCGG", bb_id="BB5"),
                TaggedBuildingBlock(tags="TGCAGTCA", bb_id="BB6"),
                TaggedBuildingBlock(tags="CTTTGTAG", bb_id="BB7"),
                TaggedBuildingBlock(tags="CAACGAAA", bb_id="BB8"),
                TaggedBuildingBlock(tags="CTGGAATC", bb_id="BB9"),
                TaggedBuildingBlock(tags="GAAGAGCT", bb_id="BB10"),
                TaggedBuildingBlock(tags="CTGTTAGA", bb_id="BB11"),
                TaggedBuildingBlock(tags="GTGAGGAC", bb_id="BB12"),
                TaggedBuildingBlock(tags="TCGCGCGC", bb_id="BB13"),
                TaggedBuildingBlock(tags="TGTCTCTA", bb_id="BB14"),
                TaggedBuildingBlock(tags="GGAGGCAA", bb_id="BB15"),
                TaggedBuildingBlock(tags="GAATACCG", bb_id="BB16"),
                TaggedBuildingBlock(tags="ATAATATC", bb_id="BB17"),
                TaggedBuildingBlock(tags="AGCTGCGC", bb_id="BB18"),
                TaggedBuildingBlock(tags="TTCACCCT", bb_id="BB19"),
                TaggedBuildingBlock(tags="GCTTAATT", bb_id="BB20"),
            ],
        )

    @pytest.fixture
    def tagged_levenshtein5_bb_set(self) -> TaggedBuildingBlockSet:
        """A TaggedBuildingBlockSet fixture for testing"""
        return TaggedBuildingBlockSet(
            bb_set_id="TEST_BB_SET",
            building_blocks=[
                TaggedBuildingBlock(tags="GATCACGGTTTG", bb_id="BB1"),
                TaggedBuildingBlock(tags="AGCGTGCTGTCA", bb_id="BB2"),
                TaggedBuildingBlock(tags="TTTCCTCTGAAT", bb_id="BB3"),
                TaggedBuildingBlock(tags="CTCCTGGGGCGT", bb_id="BB4"),
                TaggedBuildingBlock(tags="CTAATTGGGACC", bb_id="BB5"),
            ],
        )

    def test_levenshtein3_1_hashmap_corrector(self, tagged_levenshtein3_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 3 set with non-asym dist cutoff 3"""
        # set up the levenshtein 3 correctors
        leven3_1_corrector = LevenshteinDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_levenshtein3_bb_set.building_blocks for tag in bb.tags},
            distance_cutoff=1,
            asymmetrical=False,
        )

        # test the correctors
        bb = tagged_levenshtein3_bb_set.building_blocks[0]
        res = leven3_1_corrector._hash_map["TCGGTGGA"]  # direct lookup should pass
        assert res == (bb, 0)
        res = leven3_1_corrector._hash_map["TCGGTGGG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven3_1_corrector._hash_map["TCGGTGG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven3_1_corrector._hash_map["TCGGTGGAA"]  # 1 error should pass
        assert res == (bb, 1)
        with pytest.raises(KeyError):
            res = leven3_1_corrector._hash_map["TCGTGGG"]  # 2 errors should fail
        with pytest.raises(KeyError):
            res = leven3_1_corrector._hash_map["AAAAAAAA"]  # not in set should fail

    def test_levenshtein3_2_hashmap_corrector(self, tagged_levenshtein3_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 3 set with asym dist cutoff 2"""
        # set up the levenshtein 3 correctors
        with pytest.raises(HashMapCollisionError):
            LevenshteinDistBarcodeCaller(
                tag_map={tag: bb for bb in tagged_levenshtein3_bb_set.building_blocks for tag in bb.tags},
                distance_cutoff=2,
                asymmetrical=False,
            )

        leven3_2_corrector = LevenshteinDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_levenshtein3_bb_set.building_blocks for tag in bb.tags},
            distance_cutoff=2,
            asymmetrical=True,
        )

        bb = tagged_levenshtein3_bb_set.building_blocks[0]
        res = leven3_2_corrector._hash_map["TCGGTGGA"]  # direct lookup should pass
        assert res == (bb, 0)
        res = leven3_2_corrector._hash_map["TCGGTGGG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven3_2_corrector._hash_map["TCGGTGG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven3_2_corrector._hash_map["TCGGTGGAA"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven3_2_corrector._hash_map["TCGTGGG"]  # 2 errors valid asym should pass
        assert res == (bb, 2)
        res = leven3_2_corrector._hash_map["TCTGGA"]  # 2 errors valid asym should pass
        assert res == (bb, 2)
        res = leven3_2_corrector._hash_map["TCGGCCGA"]  # 2 errors valid asym should pass
        assert res == (bb, 2)
        res = leven3_2_corrector._hash_map["TCGGTGGAGG"]  # 2 errors valid asym should pass
        assert res == (bb, 2)

        with pytest.raises(KeyError):
            leven3_2_corrector._hash_map["TCGGT"]  # 3 errors should fail
        with pytest.raises(KeyError):
            leven3_2_corrector._hash_map["TCGGTAAAC"]  # 3 errors should fail
        with pytest.raises(KeyError):
            leven3_2_corrector._hash_map["TCCCCGGA"]  # 3 errors should fail
        with pytest.raises(KeyError):
            leven3_2_corrector._hash_map["AAAAAAAA"]  # not in set should fail

    def test_levenshtein5_1_hashmap_corrector(self, tagged_levenshtein5_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 5 set with non-asym dist cutoff 1"""
        # set up the levenshtein 5 correctors
        leven5_1_corrector = LevenshteinDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_levenshtein5_bb_set.building_blocks for tag in bb.tags},
            distance_cutoff=1,
            asymmetrical=False,
        )

        bb = tagged_levenshtein5_bb_set.building_blocks[0]
        res = leven5_1_corrector._hash_map["GATCACGGTTTG"]  # direct lookup should pass
        assert res == (bb, 0)
        res = leven5_1_corrector._hash_map["GATCACGTTTG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven5_1_corrector._hash_map["GATCACGGTCTG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven5_1_corrector._hash_map["GATCACGGTTTGA"]  # 1 error should pass
        assert res == (bb, 1)

        with pytest.raises(KeyError):
            leven5_1_corrector._hash_map["GATCACTTTG"]  # 2 errors should fail
        with pytest.raises(KeyError):
            leven5_1_corrector._hash_map["GATCACGGTTC"]  # 2 errors should fail
        with pytest.raises(KeyError):
            leven5_1_corrector._hash_map["GATCACGGTCCG"]  # 2 errors should fail
        with pytest.raises(KeyError):
            leven5_1_corrector._hash_map["AAAAAAAAAAAA"]  # not in set should fail

    def test_levenshtein5_2_hashmap_corrector(self, tagged_levenshtein5_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 5 set with non-asym dist cutoff 2"""
        # set up the levenshtein 5 correctors
        with pytest.raises(HashMapCollisionError):
            LevenshteinDistBarcodeCaller(
                tag_map={tag: bb for bb in tagged_levenshtein5_bb_set.building_blocks for tag in bb.tags},
                distance_cutoff=3,
                asymmetrical=False,
            )

        leven5_2_corrector = LevenshteinDistBarcodeCaller(
            tag_map={tag: bb for bb in tagged_levenshtein5_bb_set.building_blocks for tag in bb.tags},
            distance_cutoff=2,
            asymmetrical=True,
        )

        bb = tagged_levenshtein5_bb_set.building_blocks[0]
        res = leven5_2_corrector._hash_map["GATCACGGTTTG"]  # direct lookup should pass
        assert res == (bb, 0)
        res = leven5_2_corrector._hash_map["GATCACGTTTG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven5_2_corrector._hash_map["GATCACGGTCTG"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven5_2_corrector._hash_map["GATCACGGTTTGA"]  # 1 error should pass
        assert res == (bb, 1)
        res = leven5_2_corrector._hash_map["GATCACTTTG"]  # 2 errors should pass
        assert res == (bb, 2)
        res = leven5_2_corrector._hash_map["GATCACGGTTC"]  # 2 errors should pass
        assert res == (bb, 2)
        res = leven5_2_corrector._hash_map["GATCACGGTCCG"]  # 2 errors should pass
        assert res == (bb, 2)

        with pytest.raises(KeyError):
            leven5_2_corrector._hash_map["GATCACGGG"]  # 3 errors should fail
        with pytest.raises(KeyError):
            leven5_2_corrector._hash_map["GATCACGGTAAGC"]  # 3 errors should fail
        with pytest.raises(KeyError):
            leven5_2_corrector._hash_map["GATCACGGAAAG"]  # 3 errors should fail
        with pytest.raises(KeyError):
            leven5_2_corrector._hash_map["GATCCGGATTGC"]  # 3 errors should fail
        with pytest.raises(KeyError):
            leven5_2_corrector._hash_map["AAAAAAAAAAAA"]  # not in set should fail
