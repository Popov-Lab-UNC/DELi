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
                TaggedBuildingBlock(tag="TGGGCGTT", bb_id="BB1"),
                TaggedBuildingBlock(tag="GCTTGCCA", bb_id="BB2"),
                TaggedBuildingBlock(tag="TCCCAGAC", bb_id="BB3"),
                TaggedBuildingBlock(tag="GTGTTCCC", bb_id="BB4"),
                TaggedBuildingBlock(tag="AGTCATAG", bb_id="BB5"),
                TaggedBuildingBlock(tag="TAGGGAGG", bb_id="BB6"),
                TaggedBuildingBlock(tag="CTCGTCCC", bb_id="BB7"),
                TaggedBuildingBlock(tag="CGTACCAC", bb_id="BB8"),
                TaggedBuildingBlock(tag="GGGGCAAG", bb_id="BB9"),
                TaggedBuildingBlock(tag="GGCCTGTT", bb_id="BB10"),
                TaggedBuildingBlock(tag="AGCGGCGA", bb_id="BB11"),
                TaggedBuildingBlock(tag="TCGTAAGA", bb_id="BB12"),
                TaggedBuildingBlock(tag="ACTCGCTT", bb_id="BB13"),
                TaggedBuildingBlock(tag="TTTAGCTA", bb_id="BB14"),
                TaggedBuildingBlock(tag="CGAGACGG", bb_id="BB15"),
                TaggedBuildingBlock(tag="GCTCCCTA", bb_id="BB16"),
                TaggedBuildingBlock(tag="TGCACAAA", bb_id="BB17"),
                TaggedBuildingBlock(tag="CTTGTTGT", bb_id="BB18"),
                TaggedBuildingBlock(tag="CTGCTGCT", bb_id="BB19"),
                TaggedBuildingBlock(tag="TCCAGGTC", bb_id="BB20"),
            ],
        )

    @pytest.fixture
    def tagged_hamming5_bb_set(self) -> TaggedBuildingBlockSet:
        """A TaggedBuildingBlockSet fixture for testing"""
        return TaggedBuildingBlockSet(
            bb_set_id="TEST_BB_SET",
            building_blocks=[
                TaggedBuildingBlock(tag="ACCTTCAACCGT", bb_id="BB1"),
                TaggedBuildingBlock(tag="TGTGTCACGCTG", bb_id="BB2"),
                TaggedBuildingBlock(tag="GACTTTGGCTTT", bb_id="BB3"),
                TaggedBuildingBlock(tag="GGCTAATCATAG", bb_id="BB4"),
                TaggedBuildingBlock(tag="ATTTCGGGAAAC", bb_id="BB5"),
                TaggedBuildingBlock(tag="GCAAAGCCCCGT", bb_id="BB6"),
                TaggedBuildingBlock(tag="GTCAGGCTGCCA", bb_id="BB7"),
                TaggedBuildingBlock(tag="AACTGACATTGA", bb_id="BB8"),
                TaggedBuildingBlock(tag="CACGATCTTGGA", bb_id="BB9"),
                TaggedBuildingBlock(tag="CCTCCTTGATTG", bb_id="BB10"),
                TaggedBuildingBlock(tag="TCGTTCGGGGGG", bb_id="BB11"),
                TaggedBuildingBlock(tag="TCCTGCAGGCTC", bb_id="BB12"),
                TaggedBuildingBlock(tag="CTGACTGCCAAG", bb_id="BB13"),
                TaggedBuildingBlock(tag="AGGCTAGAGAAC", bb_id="BB14"),
                TaggedBuildingBlock(tag="TCCTGTTGGTAA", bb_id="BB15"),
                TaggedBuildingBlock(tag="GAAGAGTCACTT", bb_id="BB16"),
                TaggedBuildingBlock(tag="CAGTAGTCGTAT", bb_id="BB17"),
                TaggedBuildingBlock(tag="TCGAACAGCTAC", bb_id="BB18"),
                TaggedBuildingBlock(tag="CCCGACGTCTAA", bb_id="BB19"),
                TaggedBuildingBlock(tag="AGATATTTGACG", bb_id="BB20"),
            ],
        )

    @pytest.mark.unit
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

    @pytest.mark.unit
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

    @pytest.mark.unit
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

    @pytest.mark.unit
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
                TaggedBuildingBlock(tag="TCGGTGGA", bb_id="BB1"),
                TaggedBuildingBlock(tag="AGCCTCAC", bb_id="BB2"),
                TaggedBuildingBlock(tag="CTTTAATA", bb_id="BB3"),
                TaggedBuildingBlock(tag="CAATAAAT", bb_id="BB4"),
                TaggedBuildingBlock(tag="TTTGCCGG", bb_id="BB5"),
                TaggedBuildingBlock(tag="TGCAGTCA", bb_id="BB6"),
                TaggedBuildingBlock(tag="CTTTGTAG", bb_id="BB7"),
                TaggedBuildingBlock(tag="CAACGAAA", bb_id="BB8"),
                TaggedBuildingBlock(tag="CTGGAATC", bb_id="BB9"),
                TaggedBuildingBlock(tag="GAAGAGCT", bb_id="BB10"),
                TaggedBuildingBlock(tag="CTGTTAGA", bb_id="BB11"),
                TaggedBuildingBlock(tag="GTGAGGAC", bb_id="BB12"),
                TaggedBuildingBlock(tag="TCGCGCGC", bb_id="BB13"),
                TaggedBuildingBlock(tag="TGTCTCTA", bb_id="BB14"),
                TaggedBuildingBlock(tag="GGAGGCAA", bb_id="BB15"),
                TaggedBuildingBlock(tag="GAATACCG", bb_id="BB16"),
                TaggedBuildingBlock(tag="ATAATATC", bb_id="BB17"),
                TaggedBuildingBlock(tag="AGCTGCGC", bb_id="BB18"),
                TaggedBuildingBlock(tag="TTCACCCT", bb_id="BB19"),
                TaggedBuildingBlock(tag="GCTTAATT", bb_id="BB20"),
            ],
        )

    @pytest.fixture
    def tagged_levenshtein5_bb_set(self) -> TaggedBuildingBlockSet:
        """A TaggedBuildingBlockSet fixture for testing"""
        return TaggedBuildingBlockSet(
            bb_set_id="TEST_BB_SET",
            building_blocks=[
                TaggedBuildingBlock(tag="GATCACGGTTTG", bb_id="BB1"),
                TaggedBuildingBlock(tag="AGCGTGCTGTCA", bb_id="BB2"),
                TaggedBuildingBlock(tag="TTTCCTCTGAAT", bb_id="BB3"),
                TaggedBuildingBlock(tag="CTCCTGGGGCGT", bb_id="BB4"),
                TaggedBuildingBlock(tag="CTAATTGGGACC", bb_id="BB5"),
            ],
        )

    @pytest.mark.unit
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

    @pytest.mark.unit
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

    @pytest.mark.unit
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

    @pytest.mark.unit
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
