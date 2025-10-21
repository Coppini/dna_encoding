#!/usr/bin/env python3

import random
import string

import pytest

from encoding_2bit import ENCODING_TO_BASES, Encoded2bitSequence, Encoding
from encoding_2bit import decode_2bit_sequence, EncodingError, DecodingError
from encoding_2bit import bits_to_bytes, decode_2bit_sequence
from encoding_2bit import bits_to_bytes, decode_2bit_sequence


@pytest.mark.parametrize(
    "dna_sequence, expected_sequence",
    [
        ("atcg", "ATCG"),
        ("AGTC", "AGTC"),
    ],
)
def test_encode_decode_known_examples(dna_sequence, expected_sequence):
    # Act
    encoded_dna = Encoded2bitSequence.from_sequence(dna_sequence)

    # Assert
    assert encoded_dna.sequence == expected_sequence

@pytest.mark.parametrize(
    "dna_sequence",
    [
        "ATCGN",
        "ATCG-",
        "ATCGR",
        "ATCGY",
        "ATCGZ",
    ],
)
def test_encode_invalid_bases(dna_sequence):
    # Act / Assert
    with pytest.raises(EncodingError):
        Encoded2bitSequence.from_sequence(dna_sequence)


def test_encode_quality_score_length_mismatch():
    # Act / Assert
    with pytest.raises(ValueError):
        Encoded2bitSequence.from_sequence("ATCG", "!#$%&'()*+,-./01234")


@pytest.mark.parametrize("header", [True, False])
def test_encode_decode_random_fasta(header: bool):
    # Arrange
    length = random.randint(10, 100)
    available_bases = sorted(ENCODING_TO_BASES[Encoding.BIT2_ATCG])
    dna_sequence = "".join(random.choice(available_bases) for _ in range(length))
    seq_header = random.choice(string.ascii_letters) if header else None
    quality = None

    # Act
    encoded_dna = Encoded2bitSequence.from_sequence(
        dna_sequence=dna_sequence, header=seq_header
    )

    # Assert
    assert encoded_dna.sequence == dna_sequence
    assert encoded_dna.quality == quality
    # assert average_quality is float("nan")
    assert isinstance(encoded_dna.average_quality, float)
    assert not any([encoded_dna.average_quality <= 0, encoded_dna.average_quality >= 0])
    fasta_header, fasta_sequence = encoded_dna.fasta.split("\n")
    if header:
        assert fasta_header == f">{seq_header}"
    else:
        assert fasta_header.startswith(">")
    assert fasta_sequence == dna_sequence
    with pytest.raises(ValueError):
        encoded_dna.fastq


@pytest.mark.parametrize("quality1", [0, 31, 32, 64])
@pytest.mark.parametrize(
    "quality2", [0, 1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 63, 64, 127, 128, 255]
)
@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("header", [True, False])
def test_encode_decode_random_fastq(
    quality1: int, quality2: int, sort: bool, header: bool
):
    # Arrange
    length = random.randint(100, 1000)
    available_bases = sorted(ENCODING_TO_BASES[Encoding.BIT2_ATCG])
    dna_sequence = "".join(random.choice(available_bases) for _ in range(length))
    min_quality = min([quality1, quality2])
    max_quality = max([quality1, quality2])
    seq_header = random.choice(string.ascii_letters) if header else None
    quality_scores = (
        [min_quality]
        + [random.randint(min_quality, max_quality) for _ in range(length - 2)]
        + [max_quality]
    )
    if sort:
        quality_scores.sort()
    else:
        random.shuffle(quality_scores)
    average_quality = sum(quality_scores) / length
    quality = "".join(chr(score + 33) for score in quality_scores)

    # Act
    encoded_dna = Encoded2bitSequence.from_sequence(
        dna_sequence=dna_sequence, quality_score_str=quality, header=seq_header
    )

    # Assert
    assert encoded_dna.sequence == dna_sequence
    assert encoded_dna.quality == quality
    assert encoded_dna.average_quality == average_quality
    fasta_header, fasta_sequence = encoded_dna.fasta.split("\n")
    assert fasta_sequence == dna_sequence
    fastq_header, fastq_sequence, plus_signal, fastq_quality = encoded_dna.fastq.split(
        "\n"
    )
    assert fastq_sequence == dna_sequence
    assert plus_signal == "+"
    assert fastq_quality == quality
    if header:
        assert fasta_header == f">{seq_header}"
        assert fastq_header == f"@{seq_header}"
    else:
        assert fasta_header.startswith(">")
        assert fastq_header.startswith("@")

@pytest.mark.parametrize(
    "sequence, invalid_bases",
    [
        ("ATCGN", "N"),
        ("ATCGY", "Y"),
        ("ATCGZ", "Z"),
    ],
)
def test_encoding_error_invalid_bases(sequence, invalid_bases):
    with pytest.raises(EncodingError) as exc_info:
        Encoded2bitSequence.from_sequence(sequence)
    assert exc_info.value.encoding == Encoding.BIT2_ATCG
    assert f"['{invalid_bases}']" in str(exc_info.value)
    assert "Unsupported symbols in sequence" in str(exc_info.value)

def test_decode_invalid_tag():
    # Arrange: create an invalid tag (should be '01')
    invalid_bytes = bytes([0b11000000])  # '11' as tag (BIT4)
    # Act / Assert
    with pytest.raises(DecodingError):
        decode_2bit_sequence(invalid_bytes)

def test_decode_nonzero_padding():
    # Arrange: valid tag, but padding bits are not zero
    # Header: '01' + '01' (pad=2), padding bits: '11' (should be '00')
    # Sequence bits: '00' (A)
    bits = '01' + '01' + '11' + '00'
    encoded = bits_to_bytes(bits)
    with pytest.raises(DecodingError):
        decode_2bit_sequence(encoded)

def test_encode_empty_sequence():
    # Act / Assert
    encoded = Encoded2bitSequence.from_sequence("")
    assert encoded.sequence == ""
    assert encoded.quality is None
    assert isinstance(encoded.average_quality, float)
    assert not any([encoded.average_quality <= 0, encoded.average_quality >= 0])

def test_encode_decode_long_sequence():
    seq = "ATCG" * 1000
    encoded = Encoded2bitSequence.from_sequence(seq)
    assert encoded.sequence == seq
    decoded = Encoded2bitSequence.decode_sequence(encoded.encoded_sequence)
    assert decoded == seq

def test_encode_decode_lowercase():
    seq = "atcg"
    encoded = Encoded2bitSequence.from_sequence(seq)
    assert encoded.sequence == "ATCG"
    decoded = Encoded2bitSequence.decode_sequence(encoded.encoded_sequence)
    assert decoded == "ATCG"

def test_encode_decode_with_header_and_quality():
    seq = "ATCG"
    quality = "!I#$"
    header = "testheader"
    encoded = Encoded2bitSequence.from_sequence(seq, quality_score_str=quality, header=header)
    assert encoded.sequence == seq
    assert encoded.quality == quality
    assert encoded.header == header
    assert encoded.fasta.startswith(f">{header}")
    assert encoded.fastq.startswith(f"@{header}")

def test_encode_quality_score_invalid_characters():
    seq = "ATCG"
    quality = " III"  # Not valid ASCII Phred+33
    with pytest.raises(ValueError):
        Encoded2bitSequence.from_sequence(seq, quality_score_str=quality)

def test_encode_quality_score_length_mismatch_empty_quality():
    seq = "ATCG"
    quality = ""
    with pytest.raises(ValueError):
        Encoded2bitSequence.from_sequence(seq, quality_score_str=quality)

if __name__ == "__main__":
    pytest.main()
