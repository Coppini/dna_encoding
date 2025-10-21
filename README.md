# dna_encoding — Compact DNA/IUPAC encoders (2/3/4‑bit) for Python

This library provides small, composable utilities and dataclasses to **encode and decode DNA sequences**
(including IUPAC degenerate symbols, `N`, and gap characters `-`/`.`) into compact byte streams.
It is designed as an **API/toolbox**: there’s **no `main.py`**, you import what you need.

> **At a glance**
>
> - **2‑bit** encoder for `A/C/G/T`
> - **3‑bit** encoder for `A/C/G/T` + `N` + gaps `-`/`.`
> - **4‑bit** encoder for the **full IUPAC** alphabet (single, double, triple, `N`, gaps)
> - **N‑bit auto** wrapper that **chooses the smallest valid** encoding and decodes automatically
> - Lightweight **quality-score RLE** codec and ready‑to‑use dataclasses with `.fasta` / `.fastq` helpers

---

## Installation

This is a pure‑Python package; just place the module files on your `PYTHONPATH` or vendor them in your project.

```bash
# example
your_project/
  dna_encoding/
    bases.py
    code_mapping.py
    generic_encoding.py
    encoding_2bit.py
    encoding_3bit.py
    encoding_4bit.py
    encoding_nbit.py
```

---

## Supported alphabets

- **BIT2**: `A C G T`  
- **BIT3**: `A C G T N - .` 
- **BIT4 (full IUPAC)**: `A C G T R Y S W K M B D H V N - .`

(`-` and `.` are both interpreted as gaps and will both map back to `-`)

See `bases.py` for authoritative definitions.

---

## Quick start

### Encode/decode a plain ATCG string (2‑bit)

```python
from encoding_2bit import encode_2bit_sequence, decode_2bit_sequence

b = encode_2bit_sequence("ACGTAC")        # -> bytes
s = decode_2bit_sequence(b)               # "ACGTAC"
```

### Encode with `N` and gaps (3‑bit)

```python
from encoding_3bit import encode_3bit_sequence, decode_3bit_sequence

b = encode_3bit_sequence("ACGN-T.")       # -> bytes
s = decode_3bit_sequence(b)               # "ACGN-T."
```

### Full IUPAC (4‑bit)

```python
from encoding_4bit import encode_4bit_sequence, decode_4bit_sequence

b = encode_4bit_sequence("ATGCRYSWKMBDHVN-")  # -> bytes
s = decode_4bit_sequence(b)                   # same string back
```

### Auto‑pick the smallest valid encoding

```python
from encoding_nbit import encode_Nbit_sequence, decode_Nbit_sequence

enc, b = encode_Nbit_sequence("ACGTNNNN")   # enc is an Encoding enum; b are the bytes
s = decode_Nbit_sequence(b)                 # returns the decoded string
```

### Using convenience dataclasses (with optional quality and headers)

```python
from encoding_2bit import Encoded2bitSequence
from encoding_3bit import Encoded3bitSequence
from encoding_4bit import Encoded4bitSequence
from encoding_nbit import EncodedNbitSequence

# Minimal example (sequence only)
obj = Encoded2bitSequence.from_sequence("ACGTACGT")
print(obj.sequence)             # "ACGTACGT"
print(obj.fasta)                # FASTA string

# With FASTQ-quality string and header
q = "IIIIHHHH"                  # ASCII-33 encoded quality (Illumina 1.8+ style)
obj3 = Encoded3bitSequence.from_sequence("ACGTN-.-", quality_score_str=q, header="read/1")
print(obj3.fastq)               # FASTQ string

# Auto-encoding variant stores which encoding was chosen
auto = EncodedNbitSequence.from_sequence("ACGTNNNN", quality_score_str="IIIIIIII", header="auto/1")
print(auto.encoding)            # e.g., Encoding.BIT3_Ns_and_GAPs
```

---

## Module by module

### `bases.py`
- Defines the allowed symbols per encoding tier:
  - `BIT2_BASES`: `A,C,G,T`
  - `BIT3_BASES`: `BIT2 + N + - .`
  - `BIT4_BASES`: `BIT3 + IUPAC degenerates (R,Y,S,W,K,M,B,D,H,V)`
- Helpful sets for biochemical classes: `KETO_BASES (G,T)`, `PURINE_BASES (A,G)`, etc.

### `code_mapping.py`
- `Encoding` enum: `BIT2_ATCG`, `BIT3_Ns_and_GAPs`, `BIT4_FULL_IUPAC`.
- **Symbol↔bit mappings** are generated from `bases.py`:
  - `ENCODE_MAPPING[encoding][base] -> bitstring`
  - `DECODE_MAPPING[encoding][bitstring] -> base` (gaps `.` normalize to `-`)
- 3‑bit encoding reserves `STOP_3BIT="010"` as an in‑stream sentinel (not used by actual bases).
- Utility encoders `encode_base_*bit_to_string/int` (mostly internal).

### `generic_encoding.py`
- Binary helpers: `bits_to_bytes(str) -> bytes`, `bytes_to_bits(bytes) -> str`.
- **Quality codec**:
  - `encode_quality(quality_str, ascii_base=33) -> EncodedQuality`
  - `decode_quality(EncodedQuality, ascii_base=33) -> str`
  - `decode_quality_into_scores(EncodedQuality) -> iterable[int]`
  - Uses a simple **run-length encoding (RLE)** over per-base quality deltas (`min`-normalized).
- **`EncodedSequence` base dataclass** with:
  - `.from_sequence(seq, quality_score_str=None, header=None, ascii_base=33)`
  - `.sequence` (decoded string), `.quality` (decoded quality), `.quality_scores` (ints), `.average_quality`
  - `.fasta` and `.fastq` string builders
  - Subclasses must implement `encode_sequence()` and `decode_sequence()` statics.

### `encoding_2bit.py` — strict `A/C/G/T`
- Public functions: `encode_2bit_sequence(str)->bytes`, `decode_2bit_sequence(bytes)->str`.
- **Bit layout**:  
  `01 | PP | [padding zeros] | 2‑bit symbols ...`  
  - `01` is the tag; `PP` stores how many **2‑bit pairs** of zero‑padding were inserted for byte alignment (0–3).  
  - 2‑bit symbol rule: **Keto/Pyrimidine** (A=00, C=01, G=10, T=11).

- Dataclass: `Encoded2bitSequence(EncodedSequence)` with `.encode_sequence()` / `.decode_sequence()`.

### `encoding_3bit.py` — `ATCG + N + gaps`
- Public: `encode_3bit_sequence(str)->bytes`, `decode_3bit_sequence(bytes)->str`.
- **Tag**: `10`. Stream is a sequence of 3‑bit symbols. For byte alignment, the encoder may append:
  - zeros if <3 bits needed, **or**
  - `STOP_3BIT` (`010`) then zeros if ≥3 bits needed to close the byte.
- Decoder scans 3‑bit chunks after the tag until it sees the sentinel or runs out cleanly.

- Dataclass: `Encoded3bitSequence(EncodedSequence)`.

### `encoding_4bit.py` — full IUPAC (incl. `N` and gaps)
- Public: `encode_4bit_sequence(str)->bytes`, `decode_4bit_sequence(bytes)->str`.
- **Tag/header** always begins with `11`. Two header variants make the stream byte‑aligned regardless of length:
  - **Odd length**: `1100` (skip 4 bits before data)
  - **Even length**: `11110000` (skip 8 bits before data)
- Each base is a 4‑bit mask over presence of `A/C/G/T` (IUPAC degenerates are bit‑wise unions).

- Dataclass: `Encoded4bitSequence(EncodedSequence)`.

### `encoding_nbit.py` — automatic chooser
- `choose_minimal_encoding(sequence) -> Encoding`
- `encode_Nbit_sequence(sequence, encoding: Encoding|None=None) -> (Encoding, bytes)`  
  - Uppercases and strips newlines; chooses the smallest encoding that can represent all symbols.
- `decode_Nbit_sequence(encoded_bytes, encoding: Encoding|None=None) -> str`  
  - Detects the tag (`01`,`10`,`11`) and dispatches the correct decoder.
- `EncodedNbitSequence(EncodedSequence)` dataclass adds:
  - `encoding: Encoding|None`
  - `.from_sequence(...)` which also stores the chosen `encoding`
  - `.encode_self(sequence, encoding:Encoding|None=None)` updates `self.encoding`
  - `.decode_self(encoded, encoding:Encoding|None=None)` updates `self.encoding`

---

## Error handling & invariants

- All encoders **uppercase** and **strip newlines**.
- A `ValueError` is raised if the sequence contains symbols not supported by the selected encoder.
- Encoders guarantee **byte‑aligned** output; decoders assume and validate their header/tag format.
- `EncodedSequence.from_sequence` raises `ValueError` if `len(seq) != len(quality_str)` when quality provided.

---

## Recipes

### Choose encoding explicitly

```python
from encoding_nbit import EncodedNbitSequence
from code_mapping import Encoding

rec = EncodedNbitSequence.from_sequence("ACGTNNNN")
# Re-encode with 4-bit even if 3-bit would work:
forced = EncodedNbitSequence.encode_sequence("ACGTNNNN", encoding=Encoding.BIT4_FULL_IUPAC)
```

---

## FAQs

- **Why does 3‑bit have a `STOP` symbol?**  
  Because 3‑bit symbols don’t evenly pack into bytes; the sentinel lets us end cleanly while keeping the stream self‑delimiting for padding.

- **What happens to `.` vs `-`?**  
  They both decode as the gap `-` to keep a single canonical representation.

- **Can I mix encodings?**
  Each byte stream is a single encoding. Use `encoding_nbit` to detect and decode automatically.

---

## Reference

- **Tag bytes / headers**
  - BIT2: leading `01` + 2‑bit pad length (`00`..`11`), then optional zero padding, then data
  - BIT3: leading `10`, data symbols, optional STOP (`010`) then zeros
  - BIT4: headers start with `11` and indicate odd (`1100`) vs even (`11110000`) length

- **Bit meaning per symbol**
  - 2‑bit: `(Keto?, Pyrimidine?)` → A=00, C=01, G=10, T=11
  - 3‑bit: `(Concrete?, Keto?, Pyrimidine?)` with `STOP=010`
  - 4‑bit: `(A,C,G,T)` presence mask

---

## Contributing & tests

The code is intentionally small and readable.

Add unit tests that round‑trip random strings for each alphabet
(and mixed IUPAC) and verify byte‑alignment and header/tag properties.

