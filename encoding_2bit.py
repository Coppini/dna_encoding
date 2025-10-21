#!/usr/bin/env python3

from dataclasses import dataclass

from code_mapping import DECODE_MAPPING, ENCODE_MAPPING, ENCODING_TO_BASES, Encoding
from generic_encoding import (
    EncodedQuality,
    EncodedSequence,
    bits_to_bytes,
    bytes_to_bits,
    EncodingError,
    DecodingError
)

ENCODING = Encoding.BIT2_ATCG
TAG_BIT2 = "01"
HEADER_PAD_LENGTH = 2 
HEADER_LENGTH = len(TAG_BIT2) + HEADER_PAD_LENGTH

##
## From LEFT to RIGHT
## First bit (leftmost): KETO (KETO=1; AMINO=0)
## KETO=GT / AMINO=AC
##
## Second bit: PYRIMIDINE (PYRIMIDINE=1; PURINE=0)
## PYRIMIDINE=CT / PURINE=AG
##


def encode_2bit_sequence(sequence: str) -> bytes:
    """
    Layout (BIT2):
      [2b TAG=01][2b HEADER_PAD][0..6 pad bits of 0][2-bit symbols...]
    HEADER_PAD encodes the number (0..6) of zero bits added to the header
    to ensure proper byte alignment.
    """
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
        raise EncodingError(f"Unsupported symbols in sequence ({sorted(invalid_bases)})", encoding=ENCODING.value)

    mapping = ENCODE_MAPPING[ENCODING]
    data_bits = "".join(mapping[base] for base in sequence)  # 2-bit symbols

    # Compute how many *2-bit pairs* needed to reach next byte boundary
    length_before_padding = HEADER_LENGTH + len(data_bits)
    remainders = (length_before_padding % 8)
    pad_bits = (8 - remainders) if remainders else 0
    header_pad = format(pad_bits // 2, "02b")
    
    header = TAG_BIT2 + header_pad

    bitstring = header + ("0" * pad_bits) + data_bits
    return bits_to_bytes(bitstring)


def decode_2bit_sequence(encoded_bytes: bytes) -> str:
    bits = bytes_to_bits(encoded_bytes)
    if bits[:2] != TAG_BIT2:
        raise DecodingError(
            f"Wrong tag in header (found {bits[:2]}, expected {TAG_BIT2})", encoding=ENCODING.value
        )
    pad_length = int(bits[2:4], 2) * 2
    if bits[HEADER_LENGTH:HEADER_LENGTH+pad_length].strip("0"):
        raise DecodingError(
            ("Non-zero padding bits found in header"
            f" (expected '{pad_length * '0'}', found"
            f" '{bits[HEADER_LENGTH:HEADER_LENGTH+pad_length]}')."), encoding=ENCODING.value
        )
    mapping = DECODE_MAPPING[ENCODING]
    if len(seq_bits := bits[4 + pad_length:]) % 2 != 0:
        raise DecodingError(
            ("bitstring length after header is not divisible by 2"
            f" (found length {len(seq_bits)} % 2 = {len(seq_bits) % 2})."), encoding=ENCODING.value
        )
    return "".join(mapping[bits[i:i+2]] for i in range(4 + pad_length, len(bits), 2))


@dataclass
class Encoded2bitSequence(EncodedSequence):
    """
    Represents a DNA sequence with its encoding, quality scores, and header information.
    """

    encoded_sequence: bytes  # The encoded sequence
    encoded_quality: EncodedQuality | None = None  # Quality scores as bytes (optional)
    header: str | None = None  # Header information (optional)

    @staticmethod
    def encode_sequence(sequence: str) -> bytes:
        return encode_2bit_sequence(sequence)

    @staticmethod
    def decode_sequence(encoded_sequence: bytes) -> str:
        return decode_2bit_sequence(encoded_sequence)
