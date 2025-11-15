"""tests for building block error correction classes/functions"""

import pytest

from deli.decode.barcode_calling import HashMapCollisionError, LevenshteinDistBarcodeCaller, HammingDistBarcodeCaller
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
        hamming3_1_corrector = HammingDistBarcodeCaller(bb_set=tagged_hamming3_bb_set, distance_cutoff=1)

        # test the correctors
        res = hamming3_1_corrector.correct_sequence("TGGGCGTT")  # direct lookup should pass
        assert res == "TGGGCGTT"
        res = hamming3_1_corrector.correct_sequence("TGGGCGTA")  # 1 error should pass
        assert res == "TGGGCGTT"
        res = hamming3_1_corrector.correct_sequence("TGGGCGAA")  # 2 errors should fail
        assert res is None
        res = hamming3_1_corrector.correct_sequence("AAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_hamming3_2_hashmap_corrector(self, tagged_hamming3_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 3 set with asym dist cutoff 2"""
        # set up the hamming 3 correctors
        with pytest.raises(HashMapCollisionError):
            HammingDistBarcodeCaller(bb_set=tagged_hamming3_bb_set, distance_cutoff=2)
        hamming3_2_corrector = HammingDistBarcodeCaller(
            bb_set=tagged_hamming3_bb_set, distance_cutoff=2, asymmetrical=True
        )

        res = hamming3_2_corrector.correct_sequence("TGGGCGTT")  # direct lookup should pass
        assert res == "TGGGCGTT"
        res = hamming3_2_corrector.correct_sequence("TGGGCGTA")  # 1 error should pass
        assert res == "TGGGCGTT"
        res = hamming3_2_corrector.correct_sequence(
            "TGGGCGAA"
        )  # 2 errors but valid asymmetrical should pass
        assert res == "TGGGCGTT"
        res = hamming3_2_corrector.correct_sequence("TGGGCAAA")  # 3 errors should fail
        assert res is None
        res = hamming3_2_corrector.correct_sequence("AAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_hamming5_1_hashmap_corrector(self, tagged_hamming5_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 5 set with non-asym dist cutoff 1"""
        # set up the hamming 5 correctors
        hamming5_1_corrector = HammingDistBarcodeCaller(bb_set=tagged_hamming5_bb_set, distance_cutoff=1)

        res = hamming5_1_corrector.correct_sequence("ACCTTCAACCGT")  # direct lookup should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_1_corrector.correct_sequence("ACCTTCAACCGA")  # 1 error should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_1_corrector.correct_sequence("ACCTTCAACCAA")  # 2 errors should fail
        assert res is None
        res = hamming5_1_corrector.correct_sequence("AAAAAAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_hamming5_2_hashmap_corrector(self, tagged_hamming5_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 5 set with non-asym dist cutoff 2"""
        # set up the hamming 5 correctors
        hamming5_2_corrector = HammingDistBarcodeCaller(bb_set=tagged_hamming5_bb_set, distance_cutoff=2)

        res = hamming5_2_corrector.correct_sequence("ACCTTCAACCGT")  # direct lookup should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_2_corrector.correct_sequence("ACCTTCAACCGA")  # 1 error should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_2_corrector.correct_sequence("ACCTTCAACCAA")  # 2 errors should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_2_corrector.correct_sequence("ACCTTCAACAAA")  # 3 errors should fail
        assert res is None
        res = hamming5_2_corrector.correct_sequence("AAAAAAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_hamming5_3_hashmap_corrector(self, tagged_hamming5_bb_set):
        """Test HammingDistBarcodeCaller on a hamming 5 set with asym dist cutoff 3"""
        # set up the hamming 5 correctors
        with pytest.raises(HashMapCollisionError):
            HammingDistBarcodeCaller(bb_set=tagged_hamming5_bb_set, distance_cutoff=3)
        hamming5_3_corrector = HammingDistBarcodeCaller(
            bb_set=tagged_hamming5_bb_set, distance_cutoff=3, asymmetrical=True
        )

        res = hamming5_3_corrector.correct_sequence("ACCTTCAACCGT")  # direct lookup should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_3_corrector.correct_sequence("ACCTTCAACCGA")  # 1 error should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_3_corrector.correct_sequence("ACCTTCAACCAA")  # 2 errors should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_3_corrector.correct_sequence(
            "ACCTTCAACAAA"
        )  # 3 errors but valid asym should pass
        assert res == "ACCTTCAACCGT"
        res = hamming5_3_corrector.correct_sequence("ACCTTCAAAAAA")  # 4 errors should fail
        assert res is None
        res = hamming5_3_corrector.correct_sequence("AAAAAAAAAAAA")  # not in set should fail
        assert res is None


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
            bb_set=tagged_levenshtein3_bb_set, distance_cutoff=1
        )

        # test the correctors
        res = leven3_1_corrector.correct_sequence("TCGGTGGA")  # direct lookup should pass
        assert res == "TCGGTGGA"
        res = leven3_1_corrector.correct_sequence("TCGGTGGG")  # 1 error should pass
        assert res == "TCGGTGGA"
        res = leven3_1_corrector.correct_sequence("TCGGTGG")  # 1 error should pass
        assert res == "TCGGTGGA"
        res = leven3_1_corrector.correct_sequence("TCGGTGGAA")  # 1 error should pass
        assert res == "TCGGTGGA"
        res = leven3_1_corrector.correct_sequence("TCGTGGG")  # 2 errors should fail
        assert res is None
        res = leven3_1_corrector.correct_sequence("AAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_levenshtein3_2_hashmap_corrector(self, tagged_levenshtein3_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 3 set with asym dist cutoff 2"""
        # set up the levenshtein 3 correctors
        with pytest.raises(HashMapCollisionError):
            LevenshteinDistBarcodeCaller(bb_set=tagged_levenshtein3_bb_set, distance_cutoff=2)
        leven3_2_corrector = LevenshteinDistBarcodeCaller(
            bb_set=tagged_levenshtein3_bb_set, distance_cutoff=2, asymmetrical=True
        )

        res = leven3_2_corrector.correct_sequence("TCGGTGGA")  # direct lookup should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCGGTGGG")  # 1 error should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCGGTGG")  # 1 error should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCGGTGGAA")  # 1 error should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCGTGGG")  # 2 errors valid asym should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCTGGA")  # 2 errors valid asym should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCGGCCGA")  # 2 errors valid asym should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCGGTGGAGG")  # 2 errors valid asym should pass
        assert res == "TCGGTGGA"
        res = leven3_2_corrector.correct_sequence("TCGGT")  # 3 errors should fail
        assert res is None
        res = leven3_2_corrector.correct_sequence("TCGGTAAAC")  # 3 errors should fail
        assert res is None
        res = leven3_2_corrector.correct_sequence("TCCCCGGA")  # 3 errors should fail
        assert res is None
        res = leven3_2_corrector.correct_sequence("AAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_levenshtein5_1_hashmap_corrector(self, tagged_levenshtein5_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 5 set with non-asym dist cutoff 1"""
        # set up the levenshtein 5 correctors
        leven5_1_corrector = LevenshteinDistBarcodeCaller(
            bb_set=tagged_levenshtein5_bb_set, distance_cutoff=1
        )

        res = leven5_1_corrector.correct_sequence("GATCACGGTTTG")  # direct lookup should pass
        assert res == "GATCACGGTTTG"
        res = leven5_1_corrector.correct_sequence("GATCACGTTTG")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_1_corrector.correct_sequence("GATCACGGTCTG")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_1_corrector.correct_sequence("GATCACGGTTTGA")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_1_corrector.correct_sequence("GATCACTTTG")  # 2 errors should fail
        assert res is None
        res = leven5_1_corrector.correct_sequence("GATCACGGTTC")  # 2 errors should fail
        assert res is None
        res = leven5_1_corrector.correct_sequence("GATCACGGTCCG")  # 2 errors should fail
        assert res is None
        res = leven5_1_corrector.correct_sequence("AAAAAAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_levenshtein5_2_hashmap_corrector(self, tagged_levenshtein5_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 5 set with non-asym dist cutoff 2"""
        # set up the levenshtein 5 correctors
        leven5_2_corrector = LevenshteinDistBarcodeCaller(
            bb_set=tagged_levenshtein5_bb_set, distance_cutoff=2
        )

        res = leven5_2_corrector.correct_sequence("GATCACGGTTTG")  # direct lookup should pass
        assert res == "GATCACGGTTTG"
        res = leven5_2_corrector.correct_sequence("GATCACGTTTG")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_2_corrector.correct_sequence("GATCACGGTCTG")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_2_corrector.correct_sequence("GATCACGGTTTGA")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_2_corrector.correct_sequence("GATCACTTTG")  # 2 errors should pass
        assert res == "GATCACGGTTTG"
        res = leven5_2_corrector.correct_sequence("GATCACGGTTC")  # 2 errors should pass
        assert res == "GATCACGGTTTG"
        res = leven5_2_corrector.correct_sequence("GATCACGGTCCG")  # 2 errors should pass
        assert res == "GATCACGGTTTG"
        res = leven5_2_corrector.correct_sequence("GATCACGGG")  # 3 errors should fail
        assert res is None
        res = leven5_2_corrector.correct_sequence("GATCACGGTAAGC")  # 3 errors should fail
        assert res is None
        res = leven5_2_corrector.correct_sequence("GATCACGGAAAG")  # 3 errors should fail
        assert res is None
        res = leven5_2_corrector.correct_sequence("GATCCGGATTGC")  # 3 errors should fail
        assert res is None
        res = leven5_2_corrector.correct_sequence("AAAAAAAAAAAA")  # not in set should fail
        assert res is None

    @pytest.mark.unit
    def test_levenshtein5_3_hashmap_corrector(self, tagged_levenshtein5_bb_set):
        """Test LevenshteinDistBarcodeCaller on a Levenshtein 5 set with asym dist cutoff 3"""
        # set up the levenshtein 5 correctors
        with pytest.raises(HashMapCollisionError):
            LevenshteinDistBarcodeCaller(bb_set=tagged_levenshtein5_bb_set, distance_cutoff=3)
        leven5_3_corrector = LevenshteinDistBarcodeCaller(
            bb_set=tagged_levenshtein5_bb_set, distance_cutoff=3, asymmetrical=True
        )

        res = leven5_3_corrector.correct_sequence("GATCACGGTTTG")  # direct lookup should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACGTTTG")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACGGTCTG")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACGGTTTGA")  # 1 error should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACTTTG")  # 2 errors should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACGGTTC")  # 2 errors should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACGGTCCG")  # 2 errors should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACGGG")  # 3 errors valid asym should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence(
            "GATCACGGTAAGC"
        )  # 3 errors valid asym should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence(
            "GATCACGGAAAG"
        )  # 3 errors valid asym should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence(
            "GATCCGGATTGC"
        )  # 3 errors valid asym should pass
        assert res == "GATCACGGTTTG"
        res = leven5_3_corrector.correct_sequence("GATCACGG")  # 4 errors should fail
        assert res is None
        res = leven5_3_corrector.correct_sequence("GATCACGGTTTGAAAA")  # 4 errors should fail
        assert res is None
        res = leven5_3_corrector.correct_sequence("GATTTTTGTTTG")  # 4 errors should fail
        assert res is None
        res = leven5_3_corrector.correct_sequence("AGACACCGTTG")  # 4 errors should fail
        assert res is None
        res = leven5_3_corrector.correct_sequence("AAAAAAAAAAAA")  # not in set should fail
        assert res is None
