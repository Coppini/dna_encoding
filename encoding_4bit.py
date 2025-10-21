#!/usr/bin/env python3

from dataclasses import dataclass

from code_mapping import DECODE_MAPPING, ENCODE_MAPPING, ENCODING_TO_BASES, Encoding
from generic_encoding import (
    DecodingError,
    EncodedQuality,
    EncodedSequence,
    EncodingError,
    bits_to_bytes,
    bytes_to_bits,
)

ENCODING = Encoding.BIT4_FULL_IUPAC
TAG_BIT4 = "11"
ODD_TAG = f"{TAG_BIT4}00"
EVEN_TAG = f"{TAG_BIT4}11"
EVEN_TAG_AND_PADDING = f"{EVEN_TAG}0000"

##
## From LEFT to RIGHT
## First bit (leftmost): A included in base (included=1; not-included=0)
## A_INCLUDED=ANRWMDHV / A_NOT_INCLUDED=CGT-.YSKB
##
## Second bit: C included in base (included=1; not-included=0)
## C_INCLUDED=CNYSMBHV / C_NOT_INCLUDED=AGT-.RWKD
##
## Third bit: G included in base (included=1; not-included=0)
## G_INCLUDED=GNRSKBDV / G_NOT_INCLUDED=ACT-.YWMH
##
## Fourth bit: T included in base (included=1; not-included=0)
## T_INCLUDED=TNYWKBDH / T_NOT_INCLUDED=ACG-.RSMV
##


def encode_4bit_sequence(sequence: str) -> bytes:
    """
    Layout (BIT4):
      [2b TAG=11][2-bit or 6-bit padding][4-bit symbols...]
      Padding will be 2-bit if sequence is odd-length, 6-bit if even-length.
    """
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    # Determine length and required padding
    if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
        raise EncodingError(
            f"Unsupported symbols in sequence ({sorted(invalid_bases)})",
            encoding=ENCODING.value,
        )

    mapping = ENCODE_MAPPING[ENCODING]

    data_bits = "".join(mapping[base] for base in sequence)
    header = ODD_TAG if (len(sequence) % 2) else EVEN_TAG_AND_PADDING
    bitstring = header + data_bits

    return bits_to_bytes(bitstring)


def decode_4bit_sequence(encoded_bytes: bytes) -> str:
    bits = bytes_to_bits(encoded_bytes)

    # Checks if odd or even to see how much header to skip
    if bits[:4] == ODD_TAG:
        bits_to_skip = 4  # 2b TAG + 2b PAD
    elif bits[:8] == EVEN_TAG_AND_PADDING:
        bits_to_skip = 8  # 2b TAG + 6b PAD
    else:
        raise DecodingError(
            (
                f"Wrong tag in header (found {bits[:4]} or {bits[:8]},"
                f" expected {ODD_TAG} or {EVEN_TAG_AND_PADDING})"
            ),
            encoding=ENCODING.value,
        )

    rev = DECODE_MAPPING[ENCODING]
    if len(seq_bits := bits[bits_to_skip:]) % 4 > 0:
        raise DecodingError(
            (
                "Bitstring length after header is not divisible by 4"
                f" (found length {len(seq_bits)} % 4 = {len(seq_bits) % 4})."
            ),
            encoding=ENCODING.value,
        )
    return "".join(rev[seq_bits[j : j + 4]] for j in range(0, len(seq_bits), 4))


@dataclass
class Encoded4bitSequence(EncodedSequence):
    """
    Represents a DNA sequence with its encoding, quality scores, and header information.
    """

    encoded_sequence: bytes  # The encoded sequence
    encoded_quality: EncodedQuality | None = None  # Quality scores as bytes (optional)
    header: str | None = None  # Header information (optional)

    @staticmethod
    def encode_sequence(sequence: str) -> bytes:
        return encode_4bit_sequence(sequence)

    @staticmethod
    def decode_sequence(encoded_sequence: bytes) -> str:
        return decode_4bit_sequence(encoded_sequence)
