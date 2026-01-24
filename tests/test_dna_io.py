"""
Tests for deli.dna.io module.
"""

import gzip
from pathlib import Path

import dnaio
import pytest
from dnaio.exceptions import UnknownFileFormat

from deli.dna import io as dna_io


def write_fastq(path: Path, records: list[tuple[str, str, str]]):
    """Write a simple FASTQ file. records: list of (name, seq, qual)"""
    content_lines = []
    for name, seq, qual in records:
        content_lines.extend([f"@{name}", seq, "+", qual])
    path.write_text("\n".join(content_lines) + "\n")


def write_fasta(path: Path, records: list[tuple[str, str]]):
    """Write a simple FASTA file. records: list of (name, seq)"""
    content_lines = []
    for name, seq in records:
        content_lines.extend([f">{name}", seq])
    path.write_text("\n".join(content_lines) + "\n")


def write_gzip(src: Path, dst: Path):
    """Write a gzip compressed version of src to dst."""
    data = src.read_bytes()
    with gzip.open(dst, "wb") as f:
        f.write(data)


def test__check_if_seq_file_positive_and_case_insensitive():
    """Test various positive sequence file names."""
    positives = [
        "reads.fastq",
        "reads.FQ.gz",
        "sample.fasta.xz",
        "abc.fa",
        "something.CSFASTA",
        "longname.fna.bz2",
        "weird.fastq.zst",
    ]
    for name in positives:
        assert dna_io._check_if_seq_file(name) is True


def test__check_if_seq_file_negative():
    """Test various negative sequence file names."""
    negatives = ["notes.txt", "table.csv", "image.png", "archive.tar.gz"]
    for name in negatives:
        assert dna_io._check_if_seq_file(name) is False


def test_singlefile_fastq_read(tmp_path: Path):
    """Test reading a FASTQ file with SingleFileSequenceReader."""
    p = tmp_path / "reads.fastq"
    recs = [("r1", "ATGC", "IIII"), ("r2", "GGTT", "JJJJ")]
    write_fastq(p, recs)

    reader = dna_io.SingleFileSequenceReader(p)
    got = list(reader.iter_seqs())
    assert len(got) == 2
    assert got[0].name == "r1"
    assert str(got[0].sequence) == "ATGC"
    assert got[0].qualities is not None
    assert reader.get_sequence_files() == (p,)


def test_singlefile_fasta_read(tmp_path: Path):
    """Test reading a FASTA file with SingleFileSequenceReader."""
    p = tmp_path / "seqs.fasta"
    recs = [("f1", "AAA"), ("f2", "TTT")]
    write_fasta(p, recs)

    reader = dna_io.SingleFileSequenceReader(p)
    got = list(reader.iter_seqs())
    assert len(got) == 2
    assert got[1].name == "f2"
    assert str(got[1].sequence) == "TTT"
    # fasta has no qualities
    assert got[1].qualities is None


def test_iter_seqs_with_filenames(tmp_path: Path):
    """Test iter_seqs_with_filenames method of SingleFileSequenceReader."""
    p = tmp_path / "seqs.fasta"
    recs = [("x", "A"), ("y", "C")]
    write_fasta(p, recs)

    reader = dna_io.SingleFileSequenceReader(p)
    pairs = list(reader.iter_seqs_with_filenames())
    assert all(isinstance(pair[0], Path) for pair in pairs)
    assert all(isinstance(pair[1], dnaio.SequenceRecord) for pair in pairs)
    assert all(pair[0] == p for pair in pairs)


def test_constructor_raises_on_unknown_fileformat(monkeypatch, tmp_path: Path):
    """Test that SingleFileSequenceReader raises SequenceIOError on unknown format."""
    p = tmp_path / "bad.fastq"
    p.write_text("not a valid fastq")

    def fake_open(_):
        raise UnknownFileFormat("boom")

    monkeypatch.setattr(dnaio, "open", fake_open)

    with pytest.raises(dna_io.SequenceIOError):
        dna_io.SingleFileSequenceReader(p)


def test_multifile_sequence_reader_order(tmp_path: Path):
    """Test that MultiFileSequenceReader reads files in the given order."""
    a = tmp_path / "a.fastq"
    b = tmp_path / "b.fastq"
    write_fastq(a, [("a1", "A", "I")])
    write_fastq(b, [("b1", "T", "I")])

    m = dna_io.MultiFileSequenceReader([a, b])
    files = m.get_sequence_files()
    assert files == (a, b)

    names = [r.name for r in m.iter_seqs()]
    assert names == ["a1", "b1"]


def test_directory_and_glob_reader(tmp_path: Path):
    """Test SequenceDirectoryReader and SequenceGlobReader."""
    seq1 = tmp_path / "one.fastq"
    seq2 = tmp_path / "two.fasta"
    write_fastq(seq1, [("o1", "AT", "II")])
    write_fasta(seq2, [("t1", "GG")])

    # add a non-seq file
    (tmp_path / "notes.txt").write_text("hello")

    # add gz compressed version of a fastq
    seq1_gz = tmp_path / "one.fastq.gz"
    write_gzip(seq1, seq1_gz)

    dreader = dna_io.SequenceDirectoryReader(tmp_path)
    files = set(dreader.get_sequence_files())
    # should detect at least the sequence files (gz and uncompressed)
    assert any(str(f).endswith("one.fastq") or str(f).endswith("one.fastq.gz") for f in files)

    greader = dna_io.SequenceGlobReader(tmp_path / "*")
    gfiles = set(greader.get_sequence_files())
    assert len(gfiles) >= 2


def test_get_reader_various_inputs(tmp_path: Path):
    """Test dna_io.get_reader with different input types."""
    f = tmp_path / "s.fastq"
    write_fastq(f, [("s1", "A", "I")])

    # single path
    r1 = dna_io.get_reader(f)
    assert isinstance(r1, dna_io.MultiFileSequenceReader)
    assert r1.get_sequence_files() == (f,)

    # list of paths
    r2 = dna_io.get_reader([f])
    assert r2.get_sequence_files() == (f,)

    # directory
    r3 = dna_io.get_reader(tmp_path)
    assert isinstance(r3, dna_io.MultiFileSequenceReader)
    assert f in r3.get_sequence_files()

    # glob (pass a Path so it's not treated as a character sequence)
    r4 = dna_io.get_reader(tmp_path / "*.fastq")
    assert f in r4.get_sequence_files()


def test_sequence_directory_reader_empty_dir(tmp_path: Path):
    """Test that SequenceDirectoryReader raises on empty directory."""
    # create an empty dir
    empty = tmp_path / "emptydir"
    empty.mkdir()
    with pytest.raises(dna_io.SequenceIOError):
        dna_io.SequenceDirectoryReader(empty)


def test_compound_suffix_handling():
    """Test that _check_if_seq_file handles compound suffixes correctly."""
    assert dna_io._check_if_seq_file("reads.fastq.gz") is True
    assert dna_io._check_if_seq_file("reads.fasta.xz") is True
    # unknown compressor should not falsely identify
    assert dna_io._check_if_seq_file("reads.fastq.tar.gz") is False
