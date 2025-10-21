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

    mapping = ENCODE_MAPPING[ENCODING]
    data_bits = "".join(mapping[base] for base in sequence)  # 2-bit symbols

    # Compute how many *2-bit pairs* needed to reach next byte boundary
    length_before_padding = 4 + len(data_bits) # 2 from TAG + 2 from PAD_LEN info
    pad_bits = (8 - length_before_padding % 8) % 8
    pad_len = pad_bits // 2
    header = TAG_BIT2 + format(pad_len, "02b")

    bitstring = header + ("0" * pad_bits) + data_bits
    assert len(bitstring) % 8 == 0
    return bits_to_bytes(bitstring)


def decode_2bit_sequence(encoded_bytes: bytes) -> str:
    bits = bytes_to_bits(encoded_bytes)
    if bits[:2] != TAG_BIT2:
        raise ValueError("BIT2 decoder: wrong tag")
    pad_len = int(bits[2:4], 2)
    mapping = DECODE_MAPPING[ENCODING]
    decoded_bases = list()
    for i in range(4 + (pad_len*2), len(bits), 2):
        decoded_bases.append(mapping[bits[i:i+2]])
    return "".join(decoded_bases)



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
