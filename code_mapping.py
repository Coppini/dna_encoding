#!/usr/bin/env python3

from enum import IntEnum
import bases

class Encoding(IntEnum):
    BIT2_ATCG = 2
    BIT3_Ns_and_GAPs = 3
    BIT4_FULL_IUPAC = 4


ENCODING_TO_BASES = {
    Encoding.BIT2_ATCG: bases.BIT2_BASES,
    Encoding.BIT3_Ns_and_GAPs: bases.BIT3_BASES,
    Encoding.BIT4_FULL_IUPAC: bases.BIT4_BASES,
}

def _encode_base_2bit(base: str) -> int:
    b1 = base in bases.KETO_BASES        # keto? (G,T) -> 1
    b0 = base in bases.PYRIMIDINE_BASES  # pyrimidine? (C,T) -> 1
    return (b1 << 1) | b0                # A=00, C=01, G=10, T=11

def _encode_base_3bit(base: str) -> int:
    # b2 = 1 for N or gap; b2 = 0 for concrete A/C/G/T
    if base in bases.NS_AND_GAPS:
        b2 = 1
        b1 = b0 = int(base in bases.QUADRUPLE_BASES)   # N -> 1, gap('-' or '.') -> 0
        return (b2 << 2) | (b1 << 1) | b0    # N=111, gap=100
    else:
        b2 = 0
        return (b2 << 2) | _encode_base_2bit(base)  # A=000, C=001, G=010, T=011
    
def _encode_base_4bit(base: str) -> int:
    b0 = "A" in bases.BIT4_BASES[base]
    b1 = "C" in bases.BIT4_BASES[base]
    b2 = "G" in bases.BIT4_BASES[base]
    b3 = "T" in bases.BIT4_BASES[base]
    return (b3 << 3) | (b2 << 2) | (b1 << 1) | b0

def encode_base(base: str, encoding: Encoding) -> int:
    if encoding == Encoding.BIT2_ATCG:
        return _encode_base_2bit(base)
    if encoding == Encoding.BIT3_Ns_and_GAPs:
        return _encode_base_3bit(base)
    if encoding == Encoding.BIT4_FULL_IUPAC:
        return _encode_base_4bit(base)
    
def encode_base_to_bit_string(base: str, encoding: Encoding) -> str:
    return format(encode_base(base, encoding), f"0{encoding}b")

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
