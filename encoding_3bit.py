#!/usr/bin/env python3

from dataclasses import dataclass

from code_mapping import (
    DECODE_MAPPING,
    ENCODE_MAPPING,
    ENCODING_TO_BASES,
    STOP_3BIT,
    Encoding,
)
from generic_encoding import (
    EncodedQuality,
    EncodedSequence,
    bits_to_bytes,
    bytes_to_bits,
)

ENCODING = Encoding.BIT3_Ns_and_GAPs
TAG_BIT3 = "10"

##
## From LEFT to RIGHT
## First bit (leftmost): CONCRETE_BASE (CONCRETE_BASE=1; NON_STANDARD=0)
## CONCRETE_BASE=ACGT / NON_STANDARD=N-.
##
## Second bit: KETO (KETO=1; AMINO=0)
## KETO=GT / AMINO=AC
##
## Third bit: PYRIMIDINE (PYRIMIDINE=1; PURINE=0)
## PYRIMIDINE=CT / PURINE=AG
##


def encode_3bit_sequence(sequence: str) -> bytes:
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
        raise ValueError(f"Unsupported symbols for BIT3: {list(invalid_bases)}")

    mapping = ENCODE_MAPPING[ENCODING]
    data_bits = "".join(mapping[base] for base in sequence)

    bitstring = TAG_BIT3 + data_bits
    if remainders := (len(bitstring) % 8):
        to_pad = 8 - remainders
        bitstring += "0" * to_pad if to_pad < 3 else (STOP_3BIT + ("0" * (to_pad - 3)))
    assert len(bitstring) % 8 == 0

    # Byte-align with trailing zeros (decoder ignores after STOP)
    return bits_to_bytes(bitstring)


def decode_3bit_sequence(encoded_bytes: bytes) -> str:
    bits = bytes_to_bits(encoded_bytes)

    if bits[:2] != TAG_BIT3:
        raise ValueError("BIT3 decoder: wrong tag (expected '10').")

    rev = DECODE_MAPPING[ENCODING]
    decoded_bases = list()

    # Consume 3-bit symbols until done or STOP; ignore whatever remains after STOP
    for j in range(2, len(bits), 3):
        chunk = bits[j : j + 3]
        try:
            base = rev[chunk]
        except KeyError:
            if (
                chunk == STOP_3BIT
                and not bits[j + 3 :].strip("0")
                or len(chunk) < 3
                and not bits[j:].strip("0")
            ):
                break
            raise ValueError(f"Invalid 3-bit symbol {chunk} in stream")
        decoded_bases.append(base)
    return "".join(decoded_bases)


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
