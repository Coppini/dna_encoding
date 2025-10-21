#!/usr/bin/env python3

from enum import IntEnum

import bases

STOP_3BIT = "010"  # reserved sentinel; must be unused by actual bases


class Encoding(IntEnum):
    BIT2_ATCG = 2
    BIT3_Ns_and_GAPs = 3
    BIT4_FULL_IUPAC = 4


ENCODING_TO_BASES = {
    Encoding.BIT2_ATCG: bases.BIT2_BASES,
    Encoding.BIT3_Ns_and_GAPs: bases.BIT3_BASES,
    Encoding.BIT4_FULL_IUPAC: bases.BIT4_BASES,
}


def encode_base_2bit_to_int(base: str) -> int:
    b0 = base in bases.KETO_BASES  # keto? (G,T) -> 1
    b1 = base in bases.PYRIMIDINE_BASES  # pyrimidine? (C,T) -> 1
    return (b0 << 1) | b1  # A=00, C=01, G=10, T=11


def encode_base_2bit_to_string(base: str) -> str:
    b1 = int(base in bases.KETO_BASES)  # keto? (G,T) -> 1
    b0 = int(base in bases.PYRIMIDINE_BASES)  # pyrimidine? (C,T) -> 1
    return f"{b1}{b0}"  # A=00, C=01, G=10, T=11


def encode_base_3bit_to_int(base: str) -> int:
    # b2 = 1 for N or gap; b2 = 0 for concrete A/C/G/T
    if base in bases.SINGLE_BASES:
        b2 = 1
        return (b2 << 2) | encode_base_2bit_to_int(base)  # A=100, C=101, G=110, T=111
    b2 = 0
    b1 = b0 = int(base in bases.QUADRUPLE_BASES)  # N -> 1, gap('-' or '.') -> 0
    return (b2 << 2) | (b1 << 1) | b0  # N=011, gap=000


def encode_base_3bit_to_string(base: str) -> str:
    # b2 = 1 for N or gap; b2 = 0 for concrete A/C/G/T
    if base in bases.SINGLE_BASES:
        return f"1{encode_base_2bit_to_string(base)}"  # A=100, C=101, G=110, T=111
    b = int(base in bases.QUADRUPLE_BASES)
    return f"0{b}{b}"  # N -> 1, gap('-' or '.') -> 0


def encode_base_4bit_to_int(base: str) -> int:
    return sum(
        ((single_base in bases.BIT4_BASES[base]) << i)
        for i, single_base in enumerate(sorted("ACGT", reverse=True))
    )


def encode_base_4bit_to_string(base: str) -> str:
    return "".join(
        str(int(single_base in bases.BIT4_BASES[base])) for single_base in "ACGT"
    )


def encode_base(base: str, encoding: Encoding) -> int:
    if encoding == Encoding.BIT2_ATCG:
        return encode_base_2bit_to_int(base)
    if encoding == Encoding.BIT3_Ns_and_GAPs:
        return encode_base_3bit_to_int(base)
    if encoding == Encoding.BIT4_FULL_IUPAC:
        return encode_base_4bit_to_int(base)


def encode_base_to_bit_string(base: str, encoding: Encoding) -> str:
    if encoding == Encoding.BIT2_ATCG:
        return encode_base_2bit_to_string(base)
    if encoding == Encoding.BIT3_Ns_and_GAPs:
        return encode_base_3bit_to_string(base)
    if encoding == Encoding.BIT4_FULL_IUPAC:
        return encode_base_4bit_to_string(base)
    # return format(encode_base(base, encoding), f"0{encoding}b")


# Mapping dictionaries
ENCODE_MAPPING = {
    encoding_type: {
        base: encode_base_to_bit_string(base, encoding_type) for base in bases
    }
    for encoding_type, bases in ENCODING_TO_BASES.items()
}
DECODE_MAPPING = {
    encoding_type: {value: key for key, value in mapping.items() if key != "."}
    for encoding_type, mapping in ENCODE_MAPPING.items()
}
assert not any(STOP_3BIT in v.values() for v in ENCODE_MAPPING.values())
