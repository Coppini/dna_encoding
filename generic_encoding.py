#!/usr/bin/env python3

from dataclasses import dataclass
from math import ceil, sqrt
from typing import Iterable, NamedTuple

class EncodingError(ValueError):
    def __init__(self, message="Encoding error: header or bitstring is improperly formatted.", encoding: int | None = None):
        if encoding is not None:
            message = f"BIT{encoding} encoder: {message}"
        super().__init__(message)

class DecodingError(EncodingError):
    def __init__(self, message="Decoding error: bytes, header or bitstring is improperly formatted.", encoding: int | None = None):
        if encoding is not None:
            message = f"BIT{encoding} decoder error: {message}"
        super().__init__(message)

def bits_to_bytes(bitstring: str, pad: bool = False) -> bytes:
    if (remainders := len(bitstring) % 8):
        padding = (8 - remainders)
        if not pad:
            raise EncodingError(
                "bitstring length is not divisible by 8 and can't be properly formatted into bytes"
                f"{len(bitstring)} % 8 = {remainders} (would require padding with {padding} bits)."
            )
        # Pad to a whole number of bytes to avoid OverflowError
        bitstring += "0" * padding
    return int(bitstring, 2).to_bytes(len(bitstring) // 8, "big", signed=False)


def bytes_to_bits(b: bytes) -> str:
    try:
        return "".join(f"{byte:08b}" for byte in b)
    except Exception as e:
        raise DecodingError("Unable to convert bytes to bits") from e


class EncodedQuality(NamedTuple):
    encoded_quality: bytes
    minimum_quality: int

    def decode(self) -> str:
        return decode_quality(encoded_quality=self)

    def decode_into_scores(self) -> str:
        return decode_quality_into_scores(encoded_quality=self)


def encode_quality(quality_score_str: str, ascii_base: int = 33) -> EncodedQuality:
    """
    Encode quality scores string into bytes.

    Args:
        quality_score_str (str): Quality scores string.
        ascii_base (int): ASCII base value (default is 33 for Illumina 1.8+ encoding).

    Returns:
        EncodedQuality: Encoded quality scores and minimum quality, required to unencode it.
    """
    quality_scores = [ord(char) - ascii_base for char in quality_score_str]
    minimum_quality = min(quality_scores)
    maximum_quality = max(quality_scores)

    # Calculate the maximum value that can be represented using the determined bitsize
    bitsize = ceil(sqrt(maximum_quality - minimum_quality + 1))
    max_value = (2**bitsize) - 1
    encoded_quality = bytearray()
    # encoded_quality.append(bitsize)
    count = 0
    prev_score = None
    for score in quality_scores:
        encoded_score = score - minimum_quality
        # encoded = format(score - minimum_quality, f'0{bitsize}b')
        if prev_score is None:
            prev_score = encoded_score
            count = 1
        elif encoded_score != prev_score:
            encoded_quality.append(count)
            encoded_quality.append(prev_score)
            prev_score = encoded_score
            count = 1
        elif count < max_value:
            count += 1
        else:
            encoded_quality.append(count)
            encoded_quality.append(prev_score)
            prev_score = encoded_score
            count = 1

    # Handle the last group of scores
    encoded_quality.append(count)
    encoded_quality.append(prev_score)

    return EncodedQuality(bytes(encoded_quality), minimum_quality)


def decode_quality_into_scores(encoded_quality: EncodedQuality) -> Iterable[int]:
    """
    Decode quality scores bytes into a string.

    Args:
        encoded_quality (EncodedQuality): Encoded quality scores along with minimum quality.

    Returns:
        Iterable[int]: Returns a generator with the quality scores as integers.
    """
    # Extract bitsize from the first byte of encoded bytes
    # bitsize = encoded_quality.encoded_quality[0]

    # Initialize variables for tracking current score and count
    score = None
    count = 0

    # Iterate through the rest of the bytes in encoded bytes
    for i in range(0, len(encoded_quality.encoded_quality), 2):
        count = encoded_quality.encoded_quality[i]
        score = encoded_quality.encoded_quality[i + 1]
        for _ in range(count):
            yield score + encoded_quality.minimum_quality


def decode_quality(encoded_quality: EncodedQuality, ascii_base: int = 33) -> str:
    """
    Decode quality scores bytes into a string.

    Args:
        encoded_quality (EncodedQuality): Encoded quality scores along with minimum quality.
        ascii_base (int): ASCII base value (default is 33 for Illumina 1.8+ encoding).

    Returns:
        str: Decoded quality scores string.
    """
    return "".join(
        chr(quality_score + ascii_base)
        for quality_score in decode_quality_into_scores(encoded_quality)
    )


@dataclass
class EncodedSequence:
    """
    Represents a DNA sequence with its encoding, quality scores, and header information.
    """

    encoded_sequence: bytes  # The encoded sequence
    encoded_quality: EncodedQuality | None = None  # Quality scores as bytes (optional)
    header: str | None = None  # Header information (optional)

    @classmethod
    def from_sequence(
        cls,
        dna_sequence: str,
        quality_score_str: str | None = None,
        header: str | None = None,
        ascii_base: int = 33,
    ):
        """
        Create a EncodedSequence instance from a DNA sequence, quality scores, and header information.

        Args:
            dna_sequence (str): The DNA sequence.
            quality_score_str (Optional[str]): Quality scores for the sequence (optional).
                If provided, must be the same length as the DNA sequence, one quality value for each nucleotide.
            header (Optional[str]): Header information (optional).
            ascii_base (int): ASCII base value (default is 33 for Illumina 1.8+ encoding).

        Returns:
            EncodedSequence: Instance of EncodedSequence.
        """
        if quality_score_str is None:
            encoded_quality = None
        elif len(dna_sequence) != len(quality_score_str):
            raise ValueError(
                "The length of the DNA sequence differs from the length of the provide quality scores: "
                f"{len(dna_sequence)=} != {len(quality_score_str)=}"
            )
        else:
            # Encode quality scores if provided
            encoded_quality = encode_quality(quality_score_str, ascii_base=ascii_base)

        # Encode the DNA sequence
        encoded_sequence = cls.encode_sequence(dna_sequence)

        # Create and return the Encoded2bitSequence instance
        return cls(
            encoded_sequence=encoded_sequence,
            encoded_quality=encoded_quality,
            header=header,
        )

    @staticmethod
    def encode_sequence(sequence: str) -> bytes: ...

    @staticmethod
    def decode_sequence(encoded_sequence: bytes) -> str: ...

    @property
    def sequence(self) -> str:
        """
        Get the decoded DNA sequence.

        Returns:
            str: Decoded DNA sequence.
        """
        return self.decode_sequence(self.encoded_sequence)

    @property
    def quality(self) -> str | None:
        """
        Get the decoded quality scores string.

        Returns:
            Optional[str]: Decoded quality scores string if available, else None.
        """
        if self.encoded_quality:
            return decode_quality(self.encoded_quality)

    @property
    def quality_scores(self) -> list[int] | None:
        """
        Get the quality scores as a list of integers.

        Returns:
            Optional[list[int]]: List of quality scores if available, else None.
        """
        if self.encoded_quality:
            return list(decode_quality_into_scores(self.encoded_quality))

    @property
    def average_quality(self) -> float:
        """
        Calculate the average quality score.

        Returns:
            float: Average quality score.
        """
        if self.encoded_quality:
            quality_scores = list(decode_quality_into_scores(self.encoded_quality))
            return sum(quality_scores) / len(quality_scores)
        return float("NaN")

    def __str__(self) -> str:
        """
        Get a string representation of the Encoded2bitSequence.

        Returns:
            str: String representation of the Encoded2bitSequence.
        """
        if quality := self.quality:
            return f"{self.sequence}\t{quality}"
        return self.sequence

    @property
    def fasta(self) -> str:
        """
        Get the DNA sequence in FASTA format.

        Returns:
            str: DNA sequence in FASTA format.
        """
        return (
            ">"
            + (self.header or str(abs(hash(self.encoded_sequence))))
            + "\n"
            + self.sequence
        )

    @property
    def fastq(self) -> str:
        """
        Get the DNA sequence in FASTQ format.

        Returns:
            str: DNA sequence in FASTQ format.
        """
        if not (quality := self.quality):
            raise ValueError("No quality information")
        return (
            "@"
            + (self.header or str(abs(hash(self.encoded_sequence))))
            + "\n"
            + self.sequence
            + "\n+\n"
            + quality
        )
