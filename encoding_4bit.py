#!/usr/bin/env python3

from dataclasses import dataclass

from code_mapping import DECODE_MAPPING, ENCODE_MAPPING, ENCODING_TO_BASES, Encoding
from generic_encoding import EncodedQuality, EncodedSequence, bits_to_bytes, bytes_to_bits

ENCODING = Encoding.BIT4_FULL_IUPAC
TAG_BIT4 = "11"
ODD_TAG = "1100"
EVEN_TAG = "1111"
EVEN_TAG_AND_PADDING = "11110000"

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
    Encode a DNA sequence and return the sequence length, encoded bytes, and encoding type.

    Args:
        sequence (str): The DNA sequence to encode.

    Returns:
        tuple: Sequence length, Encoded bytes, and encoding type.
    """
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    # Determine length and required padding
    if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
        raise ValueError(f"Unsupported symbols for BIT4: {list(invalid_bases)}")

    mapping = ENCODE_MAPPING[ENCODING]

    data_bits = "".join(mapping[base] for base in sequence)
    odd = len(sequence) % 2
    header = ODD_TAG if odd else EVEN_TAG_AND_PADDING
    bitstring = header + data_bits

    # Always byte-aligned now
    assert len(bitstring) % 8 == 0
    return bits_to_bytes(bitstring)


def decode_4bit_sequence(encoded_bytes: bytes) -> str:
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
    bits = bytes_to_bits(encoded_bytes)

    # Checks if odd or even to see how much header to skip
    bits_to_skip = 4 if (bits[:4] == ODD_TAG) else 8

    rev = DECODE_MAPPING[ENCODING]  # maps '000'.. to bases
    seq_bits = bits[bits_to_skip:]
    assert len(seq_bits) % 4 == 0, "corrupt BIT4 alignment"

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
