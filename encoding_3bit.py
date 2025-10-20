#!/usr/bin/env python3

from dataclasses import dataclass
from math import ceil, sqrt
from typing import Iterable, NamedTuple, Optional

from code_mapping import DECODE_MAPPING, ENCODE_MAPPING, ENCODING_TO_BASES, Encoding

ENCODING = Encoding.BIT3_Ns_and_GAPs
STOP_3BIT = "101"  # reserved sentinel; must be unused by actual bases

class EncodedQuality(NamedTuple):
    encoded_quality: bytes
    minimum_quality: int

    def decode(self) -> str:
        return Encoded3bitSequence.decode_quality(encoded_quality=self)

    def decode_into_scores(self) -> str:
        return list(Encoded3bitSequence.decode_quality_into_scores(encoded_quality=self))


@dataclass
class Encoded3bitSequence:
    """
    Represents a DNA sequence with its encoding, quality scores, and header information.
    """

    encoded_sequence: bytes  # The encoded sequence
    encoded_quality: Optional[EncodedQuality] = (
        None  # Quality scores as bytes (optional)
    )
    header: Optional[str] = None  # Header information (optional)

    @classmethod
    def from_sequence(
        cls,
        dna_sequence: str,
        quality_score_str: Optional[str] = None,
        header: Optional[str] = None,
        ascii_base: int = 33,
    ):
        """
        Create a Encoded3bitSequence instance from a DNA sequence, quality scores, and header information.

        Args:
            dna_sequence (str): The DNA sequence.
            quality_score_str (Optional[str]): Quality scores for the sequence (optional).
                If provided, must be the same length as the DNA sequence, one quality value for each nucleotide.
            header (Optional[str]): Header information (optional).
            ascii_base (int): ASCII base value (default is 33 for Illumina 1.8+ encoding).


        Returns:
            Encoded3bitSequence: Instance of Encoded3bitSequence.
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
            encoded_quality = cls.encode_quality(
                quality_score_str, ascii_base=ascii_base
            )

        # Encode the DNA sequence
        encoded_sequence = cls.encode_sequence(dna_sequence)

        # Create and return the Encoded3bitSequence instance
        return cls(
            encoded_sequence=encoded_sequence,
            encoded_quality=encoded_quality,
            header=header,
        )

    @staticmethod
    def encode_sequence(sequence: str) -> bytes:
        """
        Encode a DNA sequence and return the sequence length, encoded bytes, and encoding type.

        Args:
            sequence (str): The DNA sequence to encode.

        Returns:
            tuple: Sequence length, Encoded bytes, and encoding type.
        """
        sequence = sequence.upper().replace("\n", "").replace("\r", "")
        # Determine length and required padding
        if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
            raise ValueError(f"Unsupported symbols for BIT3: {list(invalid_bases)}")
        
        mapping = ENCODE_MAPPING[ENCODING]
        # Map to 3-bit symbols
        data_bits = "".join(mapping[base] for base in sequence)

        # Add STOP sentinel
        sym_bits = data_bits + STOP_3BIT

        # Choose header (3 bits) telling how many *bit* pads precede first symbol
        # such that total is byte-aligned.
        # total_bits = 3 + lead_pad_bits + 3*len_symbols
        len_symbols = len(sym_bits) // 3
        lead_pad_bits = (- (3 + 3 * len_symbols)) % 8  # in [0..7]
        header = format(lead_pad_bits, "03b")

        bitstream = header + ("0" * lead_pad_bits) + sym_bits
        # Sanity check: must be byte-aligned now
        assert len(bitstream) % 8 == 0, "internal alignment error"

        return int(bitstream, 2).to_bytes(len(bitstream) // 8, "big", signed=False)

    @staticmethod
    def decode_sequence(encoded_bytes: bytes) -> str:
        """
        Decode an encoded DNA sequence.

        Args:
            length (int): Sequence length
            encoded_bytes (bytes): Encoded bytes of the DNA sequence.
            encoding_type (Encoding): Encoding type.

        Returns:
            str: Decoded DNA sequence.
        """
        # Convert the byte array back to a binary string
        bits = "".join(format(b, "08b") for b in encoded_bytes)
        # Header = number of lead pad *bits*
        lead_pad_bits = int(bits[:3], 2)
        i = 3 + lead_pad_bits

        rev = DECODE_MAPPING[ENCODING]  # maps '000'.. to bases
        decoded_sequence = ""
        stop_found = False
        while i + 3 <= len(bits):
            chunk = bits[i:i+3]
            if chunk == STOP_3BIT:
                stop_found = True
                break
            try:
                decoded_sequence += rev[chunk]
            except KeyError:
                raise ValueError(f"Invalid 3-bit symbol {chunk} in stream")
            i += 3
        if not stop_found:
            raise ValueError("STOP sentinel not found (corrupt BIT3 stream)")

        return decoded_sequence


    @staticmethod
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

    @staticmethod
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

    @staticmethod
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
            for quality_score in Encoded3bitSequence.decode_quality_into_scores(encoded_quality)
        )

    @property
    def sequence(self) -> str:
        """
        Get the decoded DNA sequence.

        Returns:
            str: Decoded DNA sequence.
        """
        return self.decode_sequence(self.encoded_sequence)

    @property
    def quality(self) -> Optional[str]:
        """
        Get the decoded quality scores string.

        Returns:
            Optional[str]: Decoded quality scores string if available, else None.
        """
        if self.encoded_quality:
            return self.decode_quality(self.encoded_quality)

    @property
    def quality_scores(self) -> Optional[list[int]]:
        """
        Get the quality scores as a list of integers.

        Returns:
            Optional[list[int]]: List of quality scores if available, else None.
        """
        if self.encoded_quality:
            return list(self.decode_quality_into_scores(self.encoded_quality))

    @property
    def average_quality(self) -> float:
        """
        Calculate the average quality score.

        Returns:
            float: Average quality score.
        """
        if self.encoded_quality:
            quality_scores = list(self.decode_quality_into_scores(self.encoded_quality))
            return sum(quality_scores) / len(quality_scores)
        return float("NaN")

    def __str__(self) -> str:
        """
        Get a string representation of the Encoded3bitSequence.

        Returns:
            str: String representation of the Encoded3bitSequence.
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
