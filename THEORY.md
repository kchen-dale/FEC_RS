# Theory: Reed–Solomon Erasure Codes for Packet-Level FEC

This document explains the math behind this repository: finite-field arithmetic, log tables, how the encoding matrix is built, and why decoding always works. Each step links to the code that implements it. Section 5 walks through the demo in [`test.c`](test.c) number by number.

**Contents**

1. [Finite fields](#1-finite-fields)
2. [Arithmetic in GF(2^m)](#2-arithmetic-in-gf2m)
3. [Primitive polynomials and log tables](#3-primitive-polynomials-and-log-tables)
4. [Reed–Solomon erasure coding](#4-reedsolomon-erasure-coding)
5. [Worked example: n = 6, k = 4](#5-worked-example-n--6-k--4)
6. [Mapping to the code](#6-mapping-to-the-code)
7. [References](#7-references)

**Notation**

- $k$ = number of data packets (`K`), $n$ = total packets after encoding (`N`), $L$ = bytes per packet (`PKT_LEN`).
- Indices start at 0, as in the code. Data packets are $0,\dots,k-1$, and parity packets are $k,\dots,n-1$.
- Bytes are written in hex (`0x57`). Byte matrices are printed without the `0x`, exactly as the program prints them.
- In $\mathrm{GF}(2^m)$, the symbol $+$ means field addition, which is bitwise XOR ($\oplus$).

---

## 1. Finite fields

### 1.1 What is a field?

A **field** is a set with addition and multiplication that satisfies, for all $a, b, c$:

| Property | Addition | Multiplication |
|---|---|---|
| Closure | $a+b \in F$ | $a \cdot b \in F$ |
| Associativity | $(a+b)+c = a+(b+c)$ | $(ab)c = a(bc)$ |
| Commutativity | $a+b = b+a$ | $ab = ba$ |
| Identity | $a + 0 = a$ | $a \cdot 1 = a$ (with $0 \ne 1$) |
| Inverse | $a + (-a) = 0$ | $a \cdot a^{-1} = 1$ for $a \ne 0$ |
| Distributivity | $a(b+c) = ab + ac$ | |

Informally, a field is a number system where you can add, subtract, multiply and divide (except by 0), and the usual rules of algebra hold. $\mathbb{Q}$ and $\mathbb{R}$ are fields. $\mathbb{Z}$ is not, because 2 has no integer inverse.

For coding, the important point is that **linear algebra works over any field**. Gaussian elimination, determinants, matrix inverses and the fact that a polynomial of degree $d$ has at most $d$ roots all carry over. RS codes need exactly these tools, applied to a finite set of symbols such as the 256 values of a byte.

### 1.2 Galois fields

A **Galois field** is a finite field. Two classical facts describe all of them:

- A field with $q$ elements exists **if and only if** $q = p^m$ for a prime $p$ (the *characteristic*) and an integer $m \ge 1$.
- For each such $q$, the field is unique up to relabelling. It is written $\mathrm{GF}(p^m)$.

So fields with 8 or 256 elements exist, but fields with 6 or 10 elements do not. In particular, the bytes $`\{0,\dots,255\}`$ with arithmetic **mod 256 are not a field**, because 2 has no inverse mod 256. $\mathrm{GF}(2^8)$ needs a different construction.

**Characteristic 2.** In $\mathrm{GF}(2^m)$, $1+1=0$, so $a+a=0$ and $-a=a$. **Addition and subtraction are the same operation (XOR).** That is why no minus signs appear below.

### 1.3 Constructing GF(2^m)

The elements are the polynomials of degree $< m$ with coefficients in $`\{0,1\}`$:

```math
a(x) = a_{m-1}x^{m-1} + \dots + a_1 x + a_0, \qquad a_i \in \{0,1\}.
```

- **Addition:** add coefficients mod 2.
- **Multiplication:** multiply the polynomials, then take the remainder modulo a fixed **irreducible** polynomial $P(x)$ of degree $m$. Irreducible means it cannot be factored over $\mathrm{GF}(2)$. Irreducibility plays the role primality plays for integers mod $p$: it guarantees that every non-zero element has an inverse.

Each element has three equivalent forms. Bit $i$ is the coefficient of $x^i$:

| Polynomial | Binary | Integer |
|---|---|---|
| $x^2+1$ | `101` | 5 |

The computer stores a $\mathrm{GF}(2^8)$ element as one byte. The polynomial form is only for reasoning about it.

**Example: GF(2³) with P(x) = x³ + x + 1.** Let $\alpha$ be the class of $x$, so $\alpha^3 = \alpha + 1$. Section 3 shows that every non-zero element is a power of $\alpha$. The last column below uses that result.

| Integer | Binary | Polynomial | Power of α |
|:-:|:-:|:-:|:-:|
| 0 | 000 | $0$ | — (no logarithm) |
| 1 | 001 | $1$ | $\alpha^0$ |
| 2 | 010 | $x$ | $\alpha^1$ |
| 3 | 011 | $x+1$ | $\alpha^3$ |
| 4 | 100 | $x^2$ | $\alpha^2$ |
| 5 | 101 | $x^2+1$ | $\alpha^6$ |
| 6 | 110 | $x^2+x$ | $\alpha^4$ |
| 7 | 111 | $x^2+x+1$ | $\alpha^5$ |

The integer order is **not** the power order: $\alpha^3 = 3$ but $\alpha^4 = 6$. Zero is never a power of $\alpha$.

---

## 2. Arithmetic in GF(2^m)

### 2.1 Addition and subtraction

Addition is XOR, with no carries and no reduction step:

```math
(x^2+x) + (x+1) = x^2+1 \quad\Longleftrightarrow\quad 110 \oplus 011 = 101 \quad (6 \oplus 3 = 5).
```

### 2.2 Multiplication

Multiply, then reduce mod $P(x)$. In $\mathrm{GF}(2^3)$ with $P(x) = x^3+x+1$:

```math
(x^2+x)(x+1) = x^3 + x^2 + x^2 + x = x^3 + x.
```

Since $x^3 \equiv x+1$:

```math
x^3 + x \equiv (x+1) + x = 1.
```

So $6 \cdot 3 = 1$ in this field: `6` and `3` are inverses of each other. With powers, $\alpha^4 \cdot \alpha^3 = \alpha^7 = 1$.

### 2.3 Inversion

The inverse of $a \ne 0$ is the unique $b$ with $a \cdot b \equiv 1 \pmod{P}$. There are three ways to find it:

1. **Extended Euclidean algorithm** on $a(x)$ and $P(x)$. It returns $u, v$ with $`u\,a + v\,P = 1`$, so $a^{-1} = u \bmod P$.
2. **Exponentiation.** The non-zero elements form a group of order $2^m-1$, so $a^{-1} = a^{2^m-2}$.
3. **Log table** (Section 3): $\log a^{-1} = (2^m - 1 - \log a) \bmod (2^m-1)$. This is what the code does.

*Example:* In $\mathrm{GF}(2^3)$, $x^2 = \alpha^2$, so $(x^2)^{-1} = \alpha^{7-2} = \alpha^5 = x^2+x+1$. Check: $x^2(x^2+x+1) = x^4+x^3+x^2 = (x^2+x)+(x+1)+x^2 = 1$.

### 2.4 Worked example in GF(2⁸): 0x57 · 0x83

The code uses

```math
P(x) = x^8 + x^4 + x^3 + x^2 + 1 \qquad (\texttt{0x11D}),
```

so the reduction rule is $x^8 \equiv x^4 + x^3 + x^2 + 1$.

**Step 1: write the operands as polynomials.**

```math
\texttt{0x57} = 0101\,0111_2 = x^6+x^4+x^2+x+1, \qquad \texttt{0x83} = 1000\,0011_2 = x^7+x+1.
```

**Step 2: multiply.**

```math
\begin{aligned}
(x^6+x^4+x^2+x+1)\cdot x^7 &= x^{13}+x^{11}+x^9+x^8+x^7\\
(x^6+x^4+x^2+x+1)\cdot x   &= x^7+x^5+x^3+x^2+x\\
(x^6+x^4+x^2+x+1)\cdot 1   &= x^6+x^4+x^2+x+1
\end{aligned}
```

When the three rows are added, the pairs $x^7$, $x^2$ and $x$ cancel:

```math
R(x) = x^{13}+x^{11}+x^9+x^8+x^6+x^5+x^4+x^3+1.
```

**Step 3: reduce, removing the highest term each time.**

- $x^{13} = x^5 \cdot x^8 \equiv x^9+x^8+x^7+x^5$. The pairs $x^9$, $x^8$ and $x^5$ cancel, leaving $R \equiv x^{11}+x^7+x^6+x^4+x^3+1$.
- $x^{11} = x^3 \cdot x^8 \equiv x^7+x^6+x^5+x^3$. The pairs $x^7$, $x^6$ and $x^3$ cancel, leaving $R \equiv x^5+x^4+1$.

**Step 4: convert back.** $`x^5+x^4+1 = 0011\,0001_2 =`$ **`0x31`**.

> The result depends on $P(x)$. AES uses `0x11B`, and there `0x57 · 0x83 = 0xC1`. The two fields are isomorphic but label their bytes differently, so encoder and decoder must use the same polynomial.

### 2.5 Shift-and-XOR multiplication

The hand computation above turns into a short loop. Since $`a \cdot b = \sum_i b_i\,(x^i a)`$, scan the bits of $b$ and keep multiplying $a$ by $x$ (a left shift). If the shift would set bit 8, reduce immediately by XORing in the low byte of $P$ (`0x1D`):

```c
uint8_t result = 0;
while (b) {
    if (b & 1) result ^= a;
    a = (a & 0x80) ? (a << 1) ^ 0x1D : a << 1;
    b >>= 1;
}
```

This is `gf_multiply` in [`galois.c`](galois.c). `make galois && ./galois` prints `0x31`. It takes up to 8 iterations per product. The log-table method below needs only three lookups.

---

## 3. Primitive polynomials and log tables

### 3.1 Primitive elements and primitive polynomials

The non-zero elements of $\mathrm{GF}(2^m)$ form a multiplicative group of order $2^m-1$, and this group is always **cyclic**. Some element $\alpha$, called a **primitive element**, has powers $\alpha^0, \alpha^1, \dots, \alpha^{2^m-2}$ that cover every non-zero element.

An irreducible $P(x)$ is **primitive** if its root $\alpha = x$ is a primitive element. Equivalently, the smallest $e > 0$ with $x^e \equiv 1 \pmod P$ is $e = 2^m - 1$.

Every primitive polynomial is irreducible, but not the other way round. For example, $x^4+x^3+x^2+x+1$ is irreducible, but $x^5 \equiv 1$ modulo it. Its root has order 5 instead of 15, and its powers miss ten of the non-zero elements.

Choosing a primitive $P$ means the integer `2` (the polynomial $x$) generates the whole group, so the log table can be built by repeatedly multiplying by 2. Degree 8 has exactly 16 primitive polynomials. Written as the low byte, they are:

```
0x1D 0x2B 0x2D 0x4D 0x5F 0x63 0x65 0x69 0x71 0x87 0x8D 0xA9 0xC3 0xCF 0xE7 0xF5
```

This repo uses `0x1D` (`POLYNOMIAL` in [`logarithm.h`](logarithm.h)). The AES polynomial `0x11B` is irreducible but **not** primitive: the order of `2` is only 51. That is why AES log tables use `3` as the generator.

### 3.2 The cyclic structure of GF(2³)

With $\alpha^3 = \alpha + 1$:

```math
\begin{aligned}
\alpha^4 &= \alpha\cdot\alpha^3 = \alpha^2+\alpha\\
\alpha^5 &= \alpha\cdot\alpha^4 = \alpha^3+\alpha^2 = \alpha^2+\alpha+1\\
\alpha^6 &= \alpha\cdot\alpha^5 = \alpha^3+\alpha^2+\alpha = (\alpha+1)+\alpha^2+\alpha = \alpha^2+1\\
\alpha^7 &= \alpha\cdot\alpha^6 = \alpha^3+\alpha = 1
\end{aligned}
```

The powers visit all seven non-zero elements and return to 1. As arrays:

| $k$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| `exp[k]` | 1 | 2 | 4 | 3 | 6 | 7 | 5 |

| $a$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| `log[a]` | 0 | 1 | 3 | 2 | 6 | 4 | 5 |

`log[0]` is undefined, because no power of $\alpha$ is 0.

### 3.3 Arithmetic with tables

Writing $a = \alpha^{\log a}$ turns multiplication into addition of exponents mod $2^m-1$. For $a, b \ne 0$:

```math
\begin{aligned}
a \cdot b &= \exp[(\log a + \log b) \bmod (2^m-1)]\\
a / b &= \exp[(\log a - \log b) \bmod (2^m-1)]\\
a^{-1} &= \exp[(2^m - 1 - \log a) \bmod (2^m-1)]\\
a^e &= \exp[(e \cdot \log a) \bmod (2^m-1)]
\end{aligned}
```

Zero is a special case: $0 \cdot b = 0$, $0/b = 0$, $a/0$ is an error, $0^0 = 1$, and $0^e = 0$ for $e > 0$.

- **GF(2³):** $6 \cdot 7 = \exp[(4+5) \bmod 7] = \exp[2] = 4$, i.e. $\alpha^4 \cdot \alpha^5 = \alpha^2$. Direct check: $(\alpha^2+\alpha)(\alpha^2+\alpha+1) = \alpha^4+\alpha^3+\alpha^2+\alpha^3+\alpha^2+\alpha = \alpha^4+\alpha = \alpha^2$.
- **GF(2⁸):** $\log[\texttt{0x57}] = 189$ and $\log[\texttt{0x83}] = 247$, so $(189+247) \bmod 255 = 181$ and $\exp[181] =$ `0x31`. This matches Section 2.4.

### 3.4 Implementation notes ([`logarithm.c`](logarithm.c), [`util.c`](util.c))

- **Building the tables.** Start at $g = 1$. For $i = 0..254$: set `log[g] = i`, then `g <<= 1`, and if bit 8 is set, `g ^= 0x11D`. If any $g$ repeats, the polynomial is not primitive and the program aborts.
- **No modulo in `multiply`.** `EXP_TABLE` holds the cycle twice (`exp[i] == exp[i+255]`, 510 entries). Because `log a + log b ≤ 508`, multiplication is just `EXP_TABLE[LOG_TABLE[a] + LOG_TABLE[b]]`. `divide` adds 255 when the difference is negative.

---

## 4. Reed–Solomon erasure coding

### 4.1 The problem

The sender has $k$ data packets of $L$ bytes each and transmits $n > k$ packets. The receiver knows **which** packets are missing, for example from sequence numbers. Losses at known positions are called **erasures**. We want two properties:

- **Systematic.** The first $k$ transmitted packets are the data itself, so when nothing is lost there is nothing to decode.
- **MDS (optimal).** The data can be rebuilt from **any** $k$ of the $n$ packets, so any $n-k$ losses are tolerated. No code can do better.

**Byte positions are independent.** Every operation is linear and works column by column. Byte $j$ of all $n$ packets forms one codeword over $\mathrm{GF}(2^8)$. So it is enough to describe the code for a message vector $\mathbf{m} \in \mathrm{GF}(2^8)^k$ and a codeword $\mathbf{c} \in \mathrm{GF}(2^8)^n$. Packets in one block must all have the same length $L$; shorter ones are padded.

### 4.2 Polynomials and the Vandermonde matrix

Choose $n$ **distinct** evaluation points $x_0, \dots, x_{n-1}$. A coefficient vector $\mathbf{a} = (a_0, \dots, a_{k-1})$ defines a polynomial $f(x) = a_0 + a_1x + \dots + a_{k-1}x^{k-1}$. Evaluating $f$ at all the points is a product with the $n \times k$ **Vandermonde matrix**:

```math
\begin{bmatrix} f(x_0) \\ f(x_1) \\ \vdots \\ f(x_{n-1}) \end{bmatrix}
=
\underbrace{\begin{bmatrix}
1 & x_0 & x_0^2 & \cdots & x_0^{k-1} \\
1 & x_1 & x_1^2 & \cdots & x_1^{k-1} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & x_{n-1} & x_{n-1}^2 & \cdots & x_{n-1}^{k-1}
\end{bmatrix}}_{V}
\begin{bmatrix} a_0 \\ a_1 \\ \vdots \\ a_{k-1} \end{bmatrix}
```

By convention $0^0 = 1$, so the point $x_i = 0$ is allowed and gives the row $(1, 0, \dots, 0)$.

**Lemma (any $k$ rows of $V$ are invertible).** For any set $S$ of $k$ rows,

```math
\det V_S = \prod_{i \lt j,\ i,j \in S} (x_j - x_i) \ne 0,
```

because the points are distinct. Another way to see it: if $V_S\mathbf{a} = 0$, then a polynomial of degree $< k$ has $k$ distinct roots, so it must be zero.

$\mathrm{GF}(2^8)$ has only 256 elements, so **$n \le 256$**. The code uses the points $x_i = i$ for $i = 0, \dots, n-1$ ([`vandermonde`](util.c)) and checks `1 <= K <= N <= 256` at compile time ([`rs.h`](rs.h)).

> Classical RS codes (CDs, QR codes, DVB) use the points $x_i = \alpha^i$ and a generator polynomial. Any distinct points give a (generalized) RS code with the same MDS property. The points $0, 1, \dots, n-1$ are simply the easiest to compute.

### 4.3 The systematic encoding matrix

Encoding with $\mathbf{c} = V\mathbf{m}$ is MDS but not systematic. To make it systematic, split $V$ into its top $k$ rows and the remaining $n-k$ rows, and multiply by the inverse of the top block:

```math
V = \begin{bmatrix} V_{\text{up}} \\ V_{\text{low}} \end{bmatrix},
\qquad
G = V\,V_{\text{up}}^{-1} = \begin{bmatrix} I_k \\ Q \end{bmatrix},
\qquad
Q = V_{\text{low}}\,V_{\text{up}}^{-1}.
```

This is [`get_encoding_matrix`](rs.c). $Q$ is the $(n-k) \times k$ **parity matrix**.

**What G means: interpolation.** Let $\mathbf{a} = V_{\text{up}}^{-1}\mathbf{m}$. Then $V_{\text{up}}\mathbf{a} = \mathbf{m}$ says $f(x_i) = m_i$ for $i < k$: $f$ is the unique polynomial of degree $< k$ through the $k$ data values. The codeword is

```math
\mathbf{c} = G\mathbf{m} = V\mathbf{a} = \big(f(x_0), \dots, f(x_{n-1})\big),
```

so **each parity symbol is the interpolating polynomial evaluated at an extra point**. Decoding reverses this: any $k$ values of $f$ determine $f$, and with it the data.

**Theorem (MDS).** For any set $S$ of $k$ rows, $G_S$ is invertible.

*Proof.* $`G_S = V_S\,V_{\text{up}}^{-1}`$, and both factors are invertible by the lemma. ∎

So the minimum distance is $d = n - k + 1$, the Singleton bound, and **any $n-k$ erasures can be recovered**.

> **A common pitfall.** It is tempting to stack the identity directly on the lower Vandermonde rows, $`\begin{bmatrix} I_k \\ V_{\text{low}} \end{bmatrix}`$. That matrix is **not** MDS in general. With $x_i = i$, $n = 10$ and $k = 5$, rows $`\{3,6,7,8,9\}`$ are linearly dependent, so losing packets 0, 1, 2, 4 and 5 would be unrecoverable. This mistake appeared in Plank's well-known 1997 tutorial and was corrected in 2005 [2, 3]. Multiplying by $V_{\text{up}}^{-1}$, as this repo does, is the correct construction.

### 4.4 Encoding

Stack the data packets as the rows of a $k \times L$ matrix $D$. The $n \times L$ matrix of transmitted packets is

```math
C = G\,D = \begin{bmatrix} D \\ QD \end{bmatrix}.
```

This is [`rs_encoding`](rs.c). Each parity byte is $`C_{i,j} = \bigoplus_t Q_{i-k,\,t} \cdot D_{t,j}`$. Encoding costs $`(n-k)\,k\,L`$ field multiplications, and $G$ depends only on $(n, k)$, so it is computed once.

### 4.5 Decoding

Let $S$ be the indices of any $k$ received packets, and let $R_S$ be those packets. Since $C = GD$,

```math
R_S = G_S\,D \quad\Longrightarrow\quad \boxed{D = G_S^{-1}\,R_S}
```

The steps, and the code for each:

1. Pick $k$ received indices $S$. The code takes the first $k$ survivors ([`get_matrix_with_drop`](rs.c)).
2. Keep the rows of $G$ indexed by $S$ to form $G_S$.
3. Invert $G_S$ with Gauss–Jordan elimination over $\mathrm{GF}(2^8)$ ([`invert`](util.c)). Augment $[G_S \mid I]$, reduce the left half to $I$, and read the inverse from the right half. This is $O(k^3)$ and depends only on the loss pattern ([`get_decoding_matrix`](rs.c)).
4. Compute $D = G_S^{-1} R_S$ ([`rs_decoding`](rs.c)).

**Only the lost rows need work.** Suppose data packet $t$ arrived at position $p$ of $R_S$. Then row $p$ of $G_S$ is $\mathbf{e}_t^T$, so $\mathbf{e}_p^T G_S = \mathbf{e}_t^T$, and row $t$ of $G_S^{-1}$ is $\mathbf{e}_p^T$: it just copies the received packet. Real work is needed only for the lost data packets, at $`(\#\text{lost data}) \cdot k \cdot L`$ multiplications. If only parity packets were lost, nothing needs to be done. This repo multiplies by the full $G_S^{-1}$ for simplicity.

### 4.6 Erasures vs. errors

This code corrects **erasures** (missing packets at known positions). An RS code with $n-k$ parity symbols can also correct $t$ **errors** at unknown positions plus $e$ erasures whenever $2t + e \le n-k$. That requires syndrome decoding (Berlekamp–Massey with Forney's algorithm), which is not implemented here. On packet networks, checksums (CRC, UDP) drop corrupted packets, which turns errors into erasures, so erasure decoding is enough.

---

## 5. Worked example: n = 6, k = 4

Every number below is printed by `make run` (with $P =$ `0x11D` and points $x_i = i$).

**Vandermonde matrix.** Row $i$ is $(1, i, i^2, i^3)$. For example, $3^2 = (x+1)^2 = x^2+1 =$ `05` and $3^3 = (x^2+1)(x+1) =$ `0f`.

```
V                    G = V · V_up⁻¹
01 00 00 00          01 00 00 00
01 01 01 01          00 01 00 00
01 02 04 08          00 00 01 00
01 03 05 0f          00 00 00 01
01 04 10 40          1b 1c 12 14   ← Q
01 05 11 55          1c 1b 14 12   ← Q
```

**Encoding.** The data is the ASCII string `"Fortinetprotectsthedigitalworld"` padded with one `00`:

```
D (data)                     C = G·D (transmitted)
46 6f 72 74 69 6e 65 74      46 6f 72 74 69 6e 65 74
70 72 6f 74 65 63 74 73      70 72 6f 74 65 63 74 73
74 68 65 64 69 67 69 74      74 68 65 64 69 67 69 74
61 6c 77 6f 72 6c 64 00      61 6c 77 6f 72 6c 64 00
                             c1 7c 24 88 38 48 68 fe   ← parity
                             3d 37 1b b2 46 51 31 ce   ← parity
```

Checking one byte by hand (row 4 of $G$ times column 0 of $D$):

```math
\texttt{1b}\cdot\texttt{46} \oplus \texttt{1c}\cdot\texttt{70} \oplus \texttt{12}\cdot\texttt{74} \oplus \texttt{14}\cdot\texttt{61} = \texttt{d4} \oplus \texttt{29} \oplus \texttt{fb} \oplus \texttt{c7} = \texttt{c1}.
```

**Loss and decoding.** Packets 0 and 3 are lost, so $`S = \{1, 2, 4, 5\}`$:

```
R_S (received)               G_S            G_S⁻¹ (decoding matrix)
70 72 6f 74 65 63 74 73      00 01 00 00    f5 69 f4 69   ← rebuilds packet 0
74 68 65 64 69 67 69 74      00 00 01 00    01 00 00 00   ← copies packet 1
c1 7c 24 88 38 48 68 fe      1b 1c 12 14    00 01 00 00   ← copies packet 2
3d 37 1b b2 46 51 31 ce      1c 1b 14 12    29 f5 53 8e   ← rebuilds packet 3
```

As Section 4.5 predicts, rows 1 and 2 of $G_S^{-1}$ are unit vectors that copy the received data packets. Only the rows for the lost packets mix all four received packets. For example, packet 0 is rebuilt as $\texttt{f5}\cdot R_0 \oplus \texttt{69}\cdot R_1 \oplus \texttt{f4}\cdot R_2 \oplus \texttt{69}\cdot R_3$.

$G_S^{-1} R_S$ reproduces $D$ exactly, and the program prints `RS recovered the lost packets successfully!`. By the MDS theorem, the same works for all $\binom{6}{4} = 15$ ways of receiving four of the six packets.

---

## 6. Mapping to the code

| Code | Theory |
|---|---|
| [`logarithm.c`](logarithm.c): `initializeTables` | Log/exp tables for `0x11D`, primitivity check (§3.4) |
| [`util.c`](util.c): `multiply`, `divide`, `galois_exp` | Table-based field arithmetic (§3.3) |
| [`galois.c`](galois.c): `gf_multiply` | Shift-and-XOR multiplication (§2.5) |
| [`util.c`](util.c): `vandermonde`, `get_vandermonde` | $V$ with points $x_i = i$ (§4.2) |
| [`rs.c`](rs.c): `get_encoding_matrix` | $`G = V\,V_{\text{up}}^{-1}`$ (§4.3) |
| [`rs.c`](rs.c): `rs_encoding` | $C = GD$ (§4.4) |
| [`rs.c`](rs.c): `get_encoded_matrix_after_loss` | $R_S$, the first $k$ received rows (§4.5) |
| [`rs.c`](rs.c): `get_decoding_matrix` | $G_S^{-1}$ (§4.5) |
| [`util.c`](util.c): `invert`, `gaussian_elimination` | Gauss–Jordan elimination over $\mathrm{GF}(2^8)$ |
| [`rs.c`](rs.c): `rs_decoding` | $D = G_S^{-1} R_S$ (§4.5) |
| [`rs.h`](rs.h): `K`, `N`, `PKT_LEN` | $k$, $n$, $L$ with $1 \le k \le n \le 256$ |

---

## 7. References

1. I. S. Reed and G. Solomon, "Polynomial codes over certain finite fields," *J. SIAM*, 8(2):300–304, 1960.
2. J. S. Plank, "A tutorial on Reed–Solomon coding for fault-tolerance in RAID-like systems," *Software—Practice & Experience*, 27(9):995–1012, 1997.
3. J. S. Plank and Y. Ding, "Note: Correction to the 1997 tutorial on Reed–Solomon coding," *Software—Practice & Experience*, 35(2):189–194, 2005.
4. L. Rizzo, "Effective erasure codes for reliable computer communication protocols," *ACM SIGCOMM CCR*, 27(2):24–36, 1997.
5. S. Lin and D. J. Costello, *Error Control Coding*, 2nd ed., Pearson Prentice Hall, 2004.
