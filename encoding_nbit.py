#!/usr/bin/env python3
from dataclasses import dataclass
from typing import Callable, NamedTuple

from generic_encoding import EncodedSequence, EncodedQuality, encode_quality
from code_mapping import ENCODING_TO_BASES, Encoding

from encoding_2bit import (
    encode_2bit_sequence,
    decode_2bit_sequence,
    TAG_BIT2
)
from encoding_3bit import (
    encode_3bit_sequence,
    decode_3bit_sequence,
    TAG_BIT3
)
from encoding_4bit import (
    encode_4bit_sequence,
    decode_4bit_sequence,
    TAG_BIT4
)

TAG_TO_ENCODING = {
    TAG_BIT2: Encoding.BIT2_ATCG,
    TAG_BIT3: Encoding.BIT3_Ns_and_GAPs,
    TAG_BIT4: Encoding.BIT4_FULL_IUPAC
}

def can_encode_with(encoding: Encoding, sequence: str) -> bool:
    ENCODING_TO_BASES[encoding]
    return set(sequence).issubset(ENCODING_TO_BASES[encoding])

class EncodingAndEncodedBytes(NamedTuple):
    encoding: Encoding
    encoded_bytes: bytes

class EncodingAndSequence(NamedTuple):
    encoding: Encoding
    sequence: str

def choose_minimal_encoding(sequence: str) -> Encoding:
    if can_encode_with(Encoding.BIT2_ATCG, sequence):
        return Encoding.BIT2_ATCG
    if can_encode_with(Encoding.BIT3_Ns_and_GAPs, sequence):
        return Encoding.BIT3_Ns_and_GAPs
    if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[Encoding.BIT4_FULL_IUPAC]):
        raise ValueError(f"Unsupported symbols for BIT4: {list(invalid_bases)}")
    return Encoding.BIT4_FULL_IUPAC

def detect_tag(encoded_bytes: bytes) -> Encoding:
    tag = f"{encoded_bytes[0]:08b}"[:2]
    try:
        return TAG_TO_ENCODING[tag]
    except KeyError:
        raise ValueError(f"Unknown or unsupported tag: {tag}")

def encode_Nbit_sequence(sequence: str, encoding: Encoding | None = None) -> EncodingAndEncodedBytes:
    """Choose the smallest valid encoding that supports `sequence` and encode it."""
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    encoding = encoding or choose_minimal_encoding(sequence)
    if encoding == Encoding.BIT2_ATCG:
        return encoding, encode_2bit_sequence(sequence)
    if encoding == Encoding.BIT3_Ns_and_GAPs:
        return encoding, encode_3bit_sequence(sequence)
    if encoding == Encoding.BIT4_FULL_IUPAC:
        return encoding, encode_4bit_sequence(sequence)
    raise ValueError(f"Unsupported Encoding type: {encoding}")


def decode_Nbit_sequence(encoded_bytes: bytes, encoding: Encoding | None = None) -> EncodingAndSequence:
    """
    Auto-detect the encoding by roundtripping with each decoder.
    We attempt decoders in order 2-bit, 3-bit, then 4-bit; for each,
    we decode -> re-encode with the same codec -> compare bytes.
    The first exact roundtrip match wins.
    """
    encoding = encoding or detect_tag(encoded_bytes)
    if encoding == Encoding.BIT2_ATCG:
        return decode_2bit_sequence(encoded_bytes)
    if encoding == Encoding.BIT3_Ns_and_GAPs:
        return decode_3bit_sequence(encoded_bytes)
    if encoding == Encoding.BIT4_FULL_IUPAC:
        return decode_4bit_sequence(encoded_bytes)
    raise ValueError(f"Invalid encoding: {encoding}")
    


@dataclass
class EncodedNbitSequence(EncodedSequence):
    """
    A sequence container that uses the minimal N-bit encoding automatically
    and auto-detects on decode.
    """
    encoded_sequence: bytes
    encoded_quality: EncodedQuality | None = None
    encoding: Encoding | None = None
    header: str | None = None

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
        encoding, encoded_bytes = encode_Nbit_sequence(dna_sequence)

        # Create and return the Encoded2bitSequence instance
        return cls(
            encoded_sequence=encoded_bytes,
            encoding=encoding,
            encoded_quality=encoded_quality,
            header=header,
        )

    def encode_self(self, sequence: str, encoding: Encoding | None = None) -> bytes:
        encoding, sequence = encode_Nbit_sequence(sequence, encoding=encoding or self.encoding)
        self.encoding = encoding
        return sequence
    
    @staticmethod
    def encode_sequence(sequence: str, encoding: Encoding | None = None) -> bytes:
        encoding, encoded_bytes = encode_Nbit_sequence(sequence, encoding=encoding)
        return encoded_bytes

    def decode_self(self, encoded_sequence: bytes, encoding: Encoding | None = None) -> str:
        encoding, sequence = decode_Nbit_sequence(encoded_sequence, encoding=encoding or self.encoding)
        self.encoding = encoding
        return sequence
    
    @staticmethod
    def decode_sequence(encoded_sequence: bytes, encoding: Encoding | None = None) -> str:
        return decode_Nbit_sequence(encoded_sequence, encoding=encoding)
