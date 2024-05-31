from enum import IntEnum

# Bases
REGULAR_BASES = set("ATCG")
NS_AND_GAPS = set("N-")
DEGENERATE_BASES = set("RYMKSWHBVD")


class Encoding(IntEnum):
    BIT2_ATCG = 2
    BIT3_Ns_and_GAPs = 3
    BIT4_FULL_IUPAC = 4


ENCODING_TO_BASES = {
    Encoding.BIT2_ATCG: REGULAR_BASES,
    Encoding.BIT3_Ns_and_GAPs: REGULAR_BASES.union(NS_AND_GAPS),  # E
    Encoding.BIT4_FULL_IUPAC: REGULAR_BASES.union(NS_AND_GAPS).union(DEGENERATE_BASES),
}

# Mapping dictionaries
ENCODE_MAPPING = {
    encoding_type: {
        base: format(i, f"0{encoding_type}b") for i, base in enumerate(bases)
    }
    for encoding_type, bases in ENCODING_TO_BASES.items()
}
DECODE_MAPPING = {
    encoding_type: {value: key for key, value in mapping.items()}
    for encoding_type, mapping in ENCODE_MAPPING.items()
}
