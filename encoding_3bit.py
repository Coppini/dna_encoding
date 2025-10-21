#!/usr/bin/env python3

from dataclasses import dataclass

from code_mapping import DECODE_MAPPING, ENCODE_MAPPING, ENCODING_TO_BASES, Encoding
from generic_encoding import EncodedQuality, EncodedSequence, bits_to_bytes, bytes_to_bits

ENCODING = Encoding.BIT3_Ns_and_GAPs
TAG_BIT3 = "10"
STOP_3BIT = "101"  # reserved sentinel; must be unused by actual bases

def encode_3bit_sequence(sequence: str) -> bytes:
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
        raise ValueError(f"Unsupported symbols for BIT3: {list(invalid_bases)}")

    mapping = ENCODE_MAPPING[ENCODING]
    data_bits = "".join(mapping[base] for base in sequence)

    # Choose PRE_PAD in {0,1,2} so first 3-bit symbol starts on 3-bit boundary
    for pre_pad in (0, 1, 2):
        if (5 + pre_pad) % 3 == 0:
            break

    header = TAG_BIT3 + f"{pre_pad:03b}"
    bitstring = header + ("0" * pre_pad) + data_bits + STOP_3BIT

    # Byte-align with trailing zeros (decoder ignores after STOP)
    return bits_to_bytes(bitstring)

def decode_3bit_sequence(encoded_bytes: bytes) -> str:
    bits = bytes_to_bits(encoded_bytes)

    if bits[:2] != TAG_BIT3:
        raise ValueError("BIT3 decoder: wrong tag (expected '10').")

    pre_pad = int(bits[2:5], 2)
    start = 5 + pre_pad  # first 3-bit symbol

    rev = DECODE_MAPPING[ENCODING]
    out = []
    stop_found = False

    # Consume 3-bit symbols until STOP; ignore whatever remains after STOP
    for j in range(start, len(bits), 3):
        chunk = bits[j:j+3]
        if len(chunk) < 3:
            break
        if chunk == STOP_3BIT:
            stop_found = True
            break
        base = rev.get(chunk)
        if base is None:
            raise ValueError(f"Invalid 3-bit symbol {chunk} in stream")
        out.append(base)

    if not stop_found:
        raise ValueError("STOP sentinel not found (corrupt BIT3 stream)")

    return "".join(out)


@dataclass
class Encoded3bitSequence(EncodedSequence):
    """
    Represents a DNA sequence with its encoding, quality scores, and header information.
    """

    encoded_sequence: bytes  # The encoded sequence
    encoded_quality: EncodedQuality | None = None  # Quality scores as bytes (optional)
    header: str | None = None  # Header information (optional)

    @staticmethod
    def encode_sequence(sequence: str) -> bytes:
        return encode_3bit_sequence(sequence)

    @staticmethod
    def decode_sequence(encoded_sequence: bytes) -> str:
        return decode_3bit_sequence(encoded_sequence)
