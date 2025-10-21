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
    DecodingError,
    EncodedQuality,
    EncodedSequence,
    EncodingError,
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
    """
    Layout (BIT3):
      [2b TAG=10][3-bit symbols...][0..2 pad bits of 0 or STOP+0..4 bits of 0]
    If needed, the stream is byte-aligned by adding trailing 1 or 2 zeros (which
    are ignored by the decoder as it does not add to a 3bit code), or a STOP symbol
    and then enough zeros (0..4) to reach the next byte boundary.
    """
    sequence = sequence.upper().replace("\n", "").replace("\r", "")
    if invalid_bases := set(sequence).difference(ENCODING_TO_BASES[ENCODING]):
        raise EncodingError(
            f"Unsupported symbols in sequence ({sorted(invalid_bases)})",
            encoding=ENCODING.value,
        )

    mapping = ENCODE_MAPPING[ENCODING]
    data_bits = "".join(mapping[base] for base in sequence)

    bitstring = TAG_BIT3 + data_bits
    if remainders := (len(bitstring) % 8):
        to_pad = 8 - remainders
        if to_pad < 3:
            # Just add zeros to reach byte alignment, decoder ignores trailing zeros
            # that do not add up to a full 3-bit symbol
            bitstring += "0" * to_pad
        else:
            # If there are 3 or more bits to pad, use STOP symbol + zeros
            pad_0s = to_pad - len(STOP_3BIT)
            bitstring += STOP_3BIT + ("0" * pad_0s)
    return bits_to_bytes(bitstring)


def decode_3bit_sequence(encoded_bytes: bytes) -> str:
    bits = bytes_to_bits(encoded_bytes)

    if bits[:2] != TAG_BIT3:
        raise DecodingError(
            f"Wrong tag in header (found {bits[:2]}, expected {TAG_BIT3})",
            encoding=ENCODING.value,
        )

    rev = DECODE_MAPPING[ENCODING]
    decoded_bases = list()

    # Consume 3-bit symbols until done or STOP; ignore whatever remains after STOP
    for j in range(2, len(bits), 3):
        chunk = bits[j : j + 3]
        try:
            base = rev[chunk]
        except KeyError:
            if chunk == STOP_3BIT:
                if (padding := bits[j + 3 :]).strip("0"):
                    raise DecodingError(
                        (
                            "Non-zero padding bits found after STOP symbol"
                            f" (expected all '0's, found '{padding}')."
                        ),
                        encoding=ENCODING.value,
                    )
            elif len(chunk) < 3:
                if (padding := bits[j:]).strip("0"):
                    raise DecodingError(
                        (
                            "Non-zero padding bits found at end of stream"
                            f" (expected all '0's, found '{padding}')."
                        ),
                        encoding=ENCODING.value,
                    )
            else:
                raise EncodingError(
                    f"Invalid 3-bit symbol {chunk} in stream", encoding=ENCODING.value
                )
            break
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
