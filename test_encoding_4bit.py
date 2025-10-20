#!/usr/bin/env python3

import random
import string

import pytest

from encoding_4bit import ENCODING_TO_BASES, Encoded4bitSequence, Encoding


@pytest.mark.parametrize(
    "dna_sequence, expected_sequence",
    [
        ("atcg", "ATCG"),
        ("AGTC", "AGTC"),
        ("at.cg", "AT-CG"),
        ("AGNT-C", "AGNT-C"),
        ("ARGNT-C", "ARGNT-C"),
        ("AGBNT-C", "AGBNT-C"),
    ],
)
def test_encode_decode_known_examples(dna_sequence, expected_sequence):
    # Act
    encoded_dna = Encoded4bitSequence.from_sequence(dna_sequence)

    # Assert
    assert encoded_dna.sequence == expected_sequence

@pytest.mark.parametrize(
    "dna_sequence",
    [
        "ATCGZ",
    ],
)
def test_encode_invalid_bases(dna_sequence):
    # Act / Assert
    with pytest.raises(ValueError):
        Encoded4bitSequence.from_sequence(dna_sequence)


def test_encode_quality_score_length_mismatch():
    # Act / Assert
    with pytest.raises(ValueError):
        Encoded4bitSequence.from_sequence("ATCG", "!#$%&'()*+,-./01234")


@pytest.mark.parametrize("header", [True, False])
def test_encode_decode_random_fasta(header: bool):
    # Arrange
    length = random.randint(10, 100)
    available_bases = sorted(ENCODING_TO_BASES[Encoding.BIT4_FULL_IUPAC])
    dna_sequence = "".join(random.choice(available_bases) for _ in range(length))
    expected_dna_sequence = dna_sequence.replace(".", "-")
    seq_header = random.choice(string.ascii_letters) if header else None
    quality = None

    # Act
    encoded_dna = Encoded4bitSequence.from_sequence(
        dna_sequence=dna_sequence, header=seq_header
    )

    # Assert
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
    available_bases = sorted(ENCODING_TO_BASES[Encoding.BIT4_FULL_IUPAC])
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
    encoded_dna = Encoded4bitSequence.from_sequence(
        dna_sequence=dna_sequence, quality_score_str=quality, header=seq_header
    )

    # Assert
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


if __name__ == "__main__":
    pytest.main()
