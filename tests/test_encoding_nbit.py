#!/usr/bin/env python3

import random
import string

import pytest

from encoding_nbit import ENCODING_TO_BASES, EncodedNbitSequence, Encoding, decode_Nbit_sequence

from generic_encoding import DecodingError, EncodingError, bits_to_bytes


@pytest.mark.parametrize(
    "dna_sequence, expected_sequence, encoding_type",
    [
        ("atcg", "ATCG", Encoding.BIT2_ATCG),
        ("ATCGN-A", "ATCGN-A", Encoding.BIT3_Ns_and_GAPs),
        ("ATCGRymKSWHBvd", "ATCGRYMKSWHBVD", Encoding.BIT4_FULL_IUPAC),
    ],
)
def test_encode_decode_known_examples(dna_sequence, expected_sequence, encoding_type):
    # Act
    encoded_dna = EncodedNbitSequence.from_sequence(dna_sequence)

    # Assert
    assert encoded_dna.encoding == encoding_type
    assert encoded_dna.sequence == expected_sequence


def test_encode_invalid_bases():
    # Act / Assert
    with pytest.raises(ValueError):
        EncodedNbitSequence.from_sequence("ATCGXYZ")


def test_encode_quality_score_length_mismatch():
    # Act / Assert
    with pytest.raises(ValueError):
        EncodedNbitSequence.from_sequence("ATCG", "!#$%&'()*+,-./01234")


@pytest.mark.parametrize("encoding_type", list(Encoding))
@pytest.mark.parametrize("header", [True, False])
def test_encode_decode_random_fasta(encoding_type: Encoding, header: bool):
    # Arrange
    length = random.randint(10, 100)
    available_bases = sorted(ENCODING_TO_BASES[encoding_type])
    dna_sequence = "".join(random.choice(available_bases) for _ in range(length))
    expected_dna_sequence = dna_sequence.replace(".", "-")
    seq_header = random.choice(string.ascii_letters) if header else None
    quality = None

    # Act
    encoded_dna = EncodedNbitSequence.from_sequence(
        dna_sequence=dna_sequence, header=seq_header
    )

    # Assert
    assert encoded_dna.encoding == encoding_type
    assert encoded_dna.sequence == expected_dna_sequence
    assert encoded_dna.quality == quality
    # assert average_quality is float("nan")
    assert isinstance(encoded_dna.average_quality, float)
    assert not any([encoded_dna.average_quality <= 0, encoded_dna.average_quality >= 0])
    fasta_header, fasta_sequence = encoded_dna.fasta.split("\n")
    if header:
        assert fasta_header == f">{seq_header}"
    else:
        assert fasta_header.startswith(">")
    assert fasta_sequence == expected_dna_sequence
    with pytest.raises(ValueError):
        encoded_dna.fastq


@pytest.mark.parametrize("encoding_type", list(Encoding))
@pytest.mark.parametrize("quality1", [0, 31, 32, 64])
@pytest.mark.parametrize(
    "quality2", [0, 1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 63, 64, 127, 128, 255]
)
@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("header", [True, False])
def test_encode_decode_random_fastq(
    encoding_type: Encoding, quality1: int, quality2: int, sort: bool, header: bool
):
    # Arrange
    length = random.randint(100, 1000)
    available_bases = sorted(ENCODING_TO_BASES[encoding_type])
    dna_sequence = "".join(random.choice(available_bases) for _ in range(length))
    expected_dna_sequence = dna_sequence.replace(".", "-")
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
    encoded_dna = EncodedNbitSequence.from_sequence(
        dna_sequence=dna_sequence, quality_score_str=quality, header=seq_header
    )

    # Assert
    assert encoded_dna.encoding == encoding_type
    assert encoded_dna.sequence == expected_dna_sequence
    assert encoded_dna.quality == quality
    assert encoded_dna.average_quality == average_quality
    fasta_header, fasta_sequence = encoded_dna.fasta.split("\n")
    assert fasta_sequence == expected_dna_sequence
    fastq_header, fastq_sequence, plus_signal, fastq_quality = encoded_dna.fastq.split(
        "\n"
    )
    assert fastq_sequence == expected_dna_sequence
    assert plus_signal == "+"
    assert fastq_quality == quality
    if header:
        assert fasta_header == f">{seq_header}"
        assert fastq_header == f"@{seq_header}"
    else:
        assert fasta_header.startswith(">")
        assert fastq_header.startswith("@")

def test_encoding_error_invalid_bases():
    sequence = "ATCGZ"
    invalid_bases = "Z"
    with pytest.raises(EncodingError) as exc_info:
        EncodedNbitSequence.from_sequence(sequence)
    assert exc_info.value.encoding == Encoding.BIT4_FULL_IUPAC
    assert f"['{invalid_bases}']" in str(exc_info.value)
    assert "Unsupported symbols in sequence" in str(exc_info.value)


def test_decode_invalid_tag():
    # Arrange: create an invalid tag (should be '01')
    invalid_bytes = bytes([0b00000000])  # '00' as tag (invalid TAG)
    # Act / Assert
    with pytest.raises(DecodingError):
        decode_Nbit_sequence(invalid_bytes)

def test_decode_nonzero_padding():
    # Arrange: valid tag, but padding bits are not zero
    # Header: '01' + '01' (pad=2), padding bits: '11' (should be '00')
    # Sequence bits: '00' (A)
    bits = '01' + '01' + '11' + '00'
    encoded = bits_to_bytes(bits)
    with pytest.raises(DecodingError):
        decode_Nbit_sequence(encoded)

def test_encode_empty_sequence():
    # Act / Assert
    encoded = EncodedNbitSequence.from_sequence("")
    assert encoded.sequence == ""
    assert encoded.quality is None
    assert encoded.encoding == Encoding.BIT2_ATCG
    assert isinstance(encoded.average_quality, float)
    assert not any([encoded.average_quality <= 0, encoded.average_quality >= 0])

@pytest.mark.parametrize(
    "bases, encoding_type",
    [
        ("ATCG", Encoding.BIT2_ATCG),
        ("ATCGN-.", Encoding.BIT3_Ns_and_GAPs),
        ("ACGT-.NRYSWKMBDHV", Encoding.BIT4_FULL_IUPAC),
    ],
)
def test_encode_decode_long_sequence(bases, encoding_type):
    seq = bases * 1000
    expected_seq = seq.replace(".", "-").upper()
    encoded = EncodedNbitSequence.from_sequence(seq)
    assert encoded.encoding == encoding_type
    assert encoded.sequence == expected_seq
    decoded = EncodedNbitSequence.decode_sequence(encoded.encoded_sequence)
    assert decoded == expected_seq

@pytest.mark.parametrize(
    "bases, encoding_type",
    [
        ("atcgATCG", Encoding.BIT2_ATCG),
        ("ATCGN-.atcgn", Encoding.BIT3_Ns_and_GAPs),
        ("ACGT-.NRYSWKMBDHVacgtnryswkmbdhv", Encoding.BIT4_FULL_IUPAC),
    ],
)
def test_encode_decode_lowercase(bases, encoding_type):
    seq = bases
    expected_seq = seq.replace(".", "-").upper()
    encoded = EncodedNbitSequence.from_sequence(seq)
    assert encoded.encoding == encoding_type
    assert encoded.sequence == expected_seq
    decoded = EncodedNbitSequence.decode_sequence(encoded.encoded_sequence)
    assert decoded == expected_seq

def test_encode_decode_with_header_and_quality():
    seq = "ATCGN.-YWKDHV"
    expected_seq = seq.replace(".", "-")
    quality = "!I#$III!!!!!!"
    header = "testheader"
    encoded = EncodedNbitSequence.from_sequence(seq, quality_score_str=quality, header=header)
    assert encoded.encoding == Encoding.BIT4_FULL_IUPAC
    assert encoded.sequence == expected_seq
    assert encoded.quality == quality
    assert encoded.header == header
    assert encoded.fasta.startswith(f">{header}")
    assert encoded.fastq.startswith(f"@{header}")

def test_encode_quality_score_invalid_characters():
    seq = "ATCGN.-VBDHK"
    quality = " IIIIII"  # Not valid ASCII Phred+33
    with pytest.raises(ValueError):
        EncodedNbitSequence.from_sequence(seq, quality_score_str=quality)

def test_encode_quality_score_length_mismatch_empty_quality():
    seq = "ATCGN.-RW"
    quality = ""
    with pytest.raises(ValueError):
        EncodedNbitSequence.from_sequence(seq, quality_score_str=quality)

if __name__ == "__main__":
    pytest.main()
