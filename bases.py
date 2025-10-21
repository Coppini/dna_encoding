#!/usr/bin/env python3

# Bases
## GAP
ZERO_BASES = {"-": "-", ".": "-"}
## REGULAR ATCG SINGLE BASES
SINGLE_BASES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
}

## DEGENERATE BASES WITH TWO POSSIBLE NUCLEOTIDES
DOUBLE_BASES = {
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
}

## DEGENRATE BASES WITH THREE POSSIBLE NUCLEOTIDES
TRIPLE_BASES = {
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
}

## N, DEGENERATE BASE THAT CAN BE ANY NUCLEOTIDE
QUADRUPLE_BASES = {"N": "ACGT"}

BIT2_BASES = SINGLE_BASES

NS_AND_GAPS = {
    key: value
    for key, value in list(ZERO_BASES.items()) + list(QUADRUPLE_BASES.items())
}
BIT3_BASES = {
    key: value for key, value in list(BIT2_BASES.items()) + list(NS_AND_GAPS.items())
}

DEGENERATE_BASES = {
    key: value for key, value in list(DOUBLE_BASES.items()) + list(TRIPLE_BASES.items())
}
BIT4_BASES = {
    key: value
    for key, value in list(BIT3_BASES.items()) + list(DEGENERATE_BASES.items())
}

KETO_BASES = {"G", "T"}
AMINO_BASES = {"A", "C"}
PURINE_BASES = {"A", "G"}
PYRIMIDINE_BASES = {"C", "T"}
