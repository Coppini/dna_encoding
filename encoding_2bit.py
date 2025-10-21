#!/usr/bin/env python3

from dataclasses import dataclass

from code_mapping import DECODE_MAPPING, ENCODE_MAPPING, ENCODING_TO_BASES, Encoding
from generic_encoding import EncodedQuality, EncodedSequence, bits_to_bytes, bytes_to_bits

ENCODING = Encoding.BIT2_ATCG
TAG_BIT2 = "01"

def encode_2bit_sequence(sequence: str) -> bytes:
    """
    Layout (BIT2):
      [2b TAG=01][2b TRAIL_PAD][2-bit symbols...][0..3 pad bits of 0]
    TRAIL_PAD encodes the number (0..3) of trailing zero bits added for byte alignment.
    """
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    if invalid := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
        raise ValueError(f"Unsupported symbols for BIT2: {sorted(invalid)}")

    enc = ENCODE_MAPPING[ENCODING]
    data_bits = "".join(enc[base] for base in sequence)  # 2-bit symbols

    n = len(sequence)
    lead_pairs = (2 - (n % 4)) % 4                     # <-- KEY CHANGE
    prefix = format(lead_pairs, "02b")

    bitstream = TAG_BIT2 + prefix + ("0" * (2 * lead_pairs)) + data_bits
    # now length is guaranteed multiple of 8
    return bits_to_bytes(bitstream)


def decode_2bit_sequence(encoded_bytes: bytes) -> str:
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
    if bits[:2] != TAG_BIT2:
        raise ValueError(f"BIT2 decoder: wrong tag (expected '{TAG_BIT2}').")

    lead_pairs = int(bits[2:4], 2)
    start = 4 + 2 * lead_pairs
    payload = bits[start:]

    if len(payload) % 2:
        raise ValueError("BIT2 payload not aligned to 2-bit symbols")

    dec = DECODE_MAPPING[ENCODING]
    out = []
    for i in range(0, len(payload), 2):
        sym = payload[i:i+2]
        base = dec.get(sym)
        if base is None:
            raise ValueError(f"BIT2: unknown symbol {sym}")
        out.append(base)
    return "".join(out)



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
