# FEC_RS

A small, dependency-free C implementation of **packet-level Forward Error Correction (FEC)** using a systematic Reed-Solomon erasure code over GF(2⁸).

`K` data packets are encoded into `N` packets: the original `K` plus `N - K` parity packets. The receiver can rebuild the original data from **any** `K` of the `N` packets, so up to `N - K` lost packets can be recovered.

> This is an **erasure** code: it recovers packets whose positions are known to be missing (for example, dropped UDP packets). It does not find or correct corrupted bytes at unknown positions.

## How it works

For the full math (finite fields, log tables, why any `K` packets suffice, and a byte-by-byte walkthrough of the demo), see [THEORY.md](THEORY.md).

1. **GF(2⁸) arithmetic.** At startup, log and anti-log tables are built from the primitive polynomial `x⁸ + x⁴ + x³ + x² + 1` (`0x11D`). Addition is XOR, and multiplication and division are table lookups.
2. **Encoding matrix.** An `N × K` Vandermonde matrix `V` (row `r` = `[1, r, r², …]`) is built, then multiplied by the inverse of its top `K × K` block. This produces a *systematic* matrix `E`: its top `K` rows form the identity, so the first `K` encoded packets are the original data unchanged, and the bottom `N - K` rows produce the parity.
3. **Encode.** `encoded (N × PKT_LEN) = E (N × K) · data (K × PKT_LEN)`.
4. **Lose packets.** Drop up to `N - K` rows at known indices and keep `K` of the rows that arrive.
5. **Decode.** Take the `K` rows of `E` that match the received packets, invert that `K × K` matrix with Gauss-Jordan elimination, and multiply it by the received packets to get the original data back. Any `K` rows of `E` are invertible, so this works for every loss pattern within the limit.

## Configuration

The parameters are compile-time constants in [`rs.h`](rs.h):

| Macro     | Default | Meaning                                       |
|-----------|---------|-----------------------------------------------|
| `K`       | 4       | Number of original data packets                |
| `N`       | 6       | Total packets after encoding (`N - K` parity) |
| `PKT_LEN` | 8       | Bytes per packet                              |

Constraint: `1 ≤ K ≤ N ≤ 256`, which is checked at compile time. The demo reads `K` rows of `PKT_LEN` hex bytes from `data.txt`, so update that file if you change `K` or `PKT_LEN`.

The field polynomial is `POLYNOMIAL` in [`logarithm.h`](logarithm.h). It must be primitive; a non-primitive polynomial is detected at startup and rejected.

## Build and run

Requires a C compiler (`gcc` by default) and `make`.

```bash
git clone https://github.com/kchen-dale/FEC_RS.git
cd FEC_RS
make run        # builds ./program and runs the demo
```

The demo (`test.c`) encodes the sample data in `data.txt`, drops packets `0` and `3`, decodes the remaining packets, and prints every intermediate matrix. It ends with:

```
RS recovered the lost packets successfully!
```

and exits with status `0`, or `1` on failure. Run it from the repository root so it can find `data.txt`.

Other targets:

```bash
make            # build ./program only
make galois     # build ./galois: bitwise (shift-and-XOR) GF(2^8) multiply, for comparison
make clean
```

## API

Declared in [`rs.h`](rs.h) and [`util.h`](util.h):

```c
initializeTables();                                    // call once before any GF operation

Matrix *V   = get_vandermonde();                       // N x K
Matrix *E   = get_encoding_matrix(V);                  // N x K, systematic
Matrix *enc = rs_encoding(E, data);                    // data: K x PKT_LEN -> N x PKT_LEN

int lost[] = {0, 3};
Matrix *rx  = get_encoded_matrix_after_loss(enc, lost, 2);  // K x PKT_LEN (received rows)
Matrix *D   = get_decoding_matrix(E, lost, 2);              // K x K
Matrix *out = rs_decoding(D, rx);                           // K x PKT_LEN == data

free_matrix(/* each matrix above */);
```

`get_decoding_matrix` and `get_encoded_matrix_after_loss` return `NULL` if more than `N - K` packets are reported lost. Every returned `Matrix *` is heap-allocated and must be released with `free_matrix`.

## File structure

| File | Purpose |
|------|---------|
| [`logarithm.c`](logarithm.c) / [`logarithm.h`](logarithm.h) | Build the GF(2⁸) log/exp tables from the primitive polynomial |
| [`util.c`](util.c) / [`util.h`](util.h) | GF arithmetic, the `Matrix` type, multiplication, Gauss-Jordan inversion, file input |
| [`rs.c`](rs.c) / [`rs.h`](rs.h) | RS encoding/decoding matrices and the encode/decode API; `K`, `N`, `PKT_LEN` |
| [`test.c`](test.c) | End-to-end demo: encode, drop, decode, verify |
| [`THEORY.md`](THEORY.md) | The math behind the code, with a worked example |
| [`data.txt`](data.txt) | Sample input: `K` packets of `PKT_LEN` hex bytes |
| [`galois.c`](galois.c) | Standalone bitwise GF(2⁸) multiply, for comparison with the table method |
| [`Makefile`](Makefile) | Build rules |

## Limitations

- `K`, `N` and `PKT_LEN` are fixed at compile time.
- Matrix multiplication is a straightforward triple loop with no SIMD or table-driven bulk operations, so it is written for clarity rather than throughput.
- Errors such as a singular matrix or failed allocation print a message and call `exit()` instead of returning an error code.
