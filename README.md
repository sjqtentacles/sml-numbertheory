# sml-numbertheory

[![CI](https://github.com/sjqtentacles/sml-numbertheory/actions/workflows/ci.yml/badge.svg)](https://github.com/sjqtentacles/sml-numbertheory/actions/workflows/ci.yml)

A number-theory toolkit for Standard ML: prime sieving, integer
factorization, modular inverses, the Chinese Remainder Theorem, Euler's
totient, the Jacobi/Legendre symbols, modular square roots
(Tonelli–Shanks), primitive roots, and the Möbius function — pure,
deterministic, and byte-identical under [MLton](http://mlton.org/) and
[Poly/ML](https://www.polyml.org/).

It is layered on the vendored
[`sml-bigint`](https://github.com/sjqtentacles/sml-bigint), so every result is
exact and arbitrary-precision. The three heavy primitives — the GCD, modular
exponentiation, and Miller-Rabin primality — are taken **directly** from
`sml-bigint` (`BigInt.gcd`, `BigInt.modpow`, `BigInt.isProbablePrime`); this
library never re-implements them. What it adds is the higher-level toolkit
built on top.

Everything is pure Standard ML over the Basis library: no FFI, threads, clock,
or RNG, so a given input always produces the same output.

## API

```sml
signature NUMBER_THEORY =
sig
  val sieve      : int -> int list                                  (* primes <= n *)
  val factor     : BigInt.int -> BigInt.int list                    (* sorted, with multiplicity *)
  val modInverse : BigInt.int * BigInt.int -> BigInt.int option     (* extended Euclid *)
  val crt        : (BigInt.int * BigInt.int) list -> BigInt.int option
  val eulerPhi   : BigInt.int -> BigInt.int                         (* Euler's totient *)
  val jacobi     : BigInt.int * BigInt.int -> int                   (* ~1 | 0 | 1 *)
  val legendre   : BigInt.int * BigInt.int -> int                   (* ~1 | 0 | 1 *)
  val tonelliShanks : BigInt.int * BigInt.int -> BigInt.int option  (* sqrt mod prime *)
  val primitiveRoot : BigInt.int -> BigInt.int option               (* smallest generator *)
  val moebius       : BigInt.int -> int                             (* ~1 | 0 | 1 *)
end
```

| Function | What it does | Notes |
| --- | --- | --- |
| `sieve n` | All primes `2 <= p <= n`, ascending | sieve of Eratosthenes; `[]` for `n < 2` |
| `factor n` | Prime factors of `n >= 1`, non-decreasing, product `= n` | trial division by small primes, then Pollard's rho; primality via Miller-Rabin. `Domain` for `n <= 0` |
| `modInverse (a, m)` | `SOME x`, `0 <= x < \|m\|`, `a*x = 1 (mod m)`, else `NONE` | extended Euclid; `Domain` when `m = 0` |
| `crt [(a1,m1), ...]` | Unique `x` in `[0, M)` with `x = ai (mod mi)`, else `NONE` | general (non-coprime) CRT; `[]` → `SOME 0`; each `mi >= 1` (`Domain` otherwise) |
| `eulerPhi n` | Euler's totient `phi(n)`, `n >= 1` | from the factorization; `Domain` for `n <= 0` |
| `jacobi (a, n)` | Jacobi symbol `(a/n)` for odd `n >= 1` | `Domain` when `n` is even or `< 1` |
| `legendre (a, p)` | Legendre symbol `(a/p)` for an odd prime `p` | coincides with Jacobi on a prime modulus |
| `tonelliShanks (n, p)` | `SOME r`, `0 <= r < p`, `r*r = n (mod p)`, else `NONE` | modular square root mod the odd prime `p`; `NONE` for a non-residue; returns the smaller of the two roots (deterministic); `n = 0` → `SOME 0` |
| `primitiveRoot p` | `SOME g`, the smallest primitive root mod `p`, else `NONE` | a generator of the units mod `p`; intended for prime `p` (where one always exists); `Domain` for `p < 1` |
| `moebius n` | Möbius function `μ(n)`, `n >= 1`: `1`/`~1`/`0` | `~1`/`1` by parity of the (distinct) prime count when square-free, `0` otherwise; `μ(1)=1`; `Domain` for `n <= 0` |

## How it works

- **`sieve`** marks composites in a boolean array; only odd multiples from
  `i*i` upward are struck, so it stays linear in the array size.
- **`factor`** first peels off every prime factor below 1000 by trial
  division (this alone clears the common small-factor cases), then splits the
  remaining cofactor with **Pollard's rho** (`f(x) = x^2 + c (mod n)`, Floyd
  cycle detection, retried with a fresh `c` on collapse). Each candidate is
  classified with `BigInt.isProbablePrime` (deterministic Miller-Rabin), so
  the returned factors are genuinely prime and their product is the input.
- **`modInverse`** runs the **extended Euclidean algorithm** to recover the
  Bézout coefficient; the inverse exists iff `gcd(a, m) = 1`.
- **`crt`** folds the congruences pairwise with the general CRT merge, which
  handles non-coprime moduli and reports `NONE` on an inconsistent system.
- **`eulerPhi`** multiplies `(p-1) * p^(k-1)` over the prime powers `p^k` in
  the factorization.
- **`jacobi`** is the classic reduction loop (strip factors of two using the
  `n mod 8` rule, then quadratic reciprocity); `legendre` delegates to it.
- **`tonelliShanks`** screens with Euler's criterion (`NONE` for a
  non-residue), writes `p-1 = q*2^s`, picks the smallest quadratic
  non-residue as a 2-power witness, and runs the **Tonelli–Shanks** descent;
  the `p ≡ 3 (mod 4)` shortcut `n^((p+1)/4)` falls out as the `s = 1` case.
  Of the two roots `r` and `p-r` the smaller is returned, so the answer is
  deterministic.
- **`primitiveRoot`** factors `φ(p)` once and accepts the least `g` (coprime
  to `p`) for which `g^(φ(p)/q) ≠ 1` for every prime `q | φ(p)` — i.e. `g` has
  full order `φ(p)`.
- **`moebius`** reads the (sorted, with-multiplicity) factorization: `0` as
  soon as a prime repeats, otherwise `(-1)` raised to the number of distinct
  primes.

## Example

```sml
structure NT = NumberTheory
structure B  = BigInt

val n  = B.mul (valOf (B.fromString "2147483647"),         (* 2^31 - 1 *)
                valOf (B.fromString "2305843009213693951")) (* 2^61 - 1 *)
val fs = NT.factor n          (* [2147483647, 2305843009213693951] *)

val inv = NT.modInverse (B.fromInt 3, B.fromInt 11)         (* SOME 4 *)
val x   = NT.crt [ (B.fromInt 2, B.fromInt 3)
                 , (B.fromInt 3, B.fromInt 5)
                 , (B.fromInt 2, B.fromInt 7) ]             (* SOME 23 *)
val phi = NT.eulerPhi (B.fromInt 36)                        (* 12 *)

val r   = NT.tonelliShanks (B.fromInt 10, B.fromInt 13)     (* SOME 6 *)
val g   = NT.primitiveRoot (B.fromInt 23)                   (* SOME 5 *)
val mu  = NT.moebius (B.fromInt 30)                         (* ~1 *)
```

`examples/demo.sml` runs a fuller tour; `make example` builds and runs it:

```
sml-numbertheory demo
=====================

primes up to 30 : 2 3 5 7 11 13 17 19 23 29
pi(1000)        = 168

n = (2^31-1)(2^61-1) = 4951760154835678088235319297
factor n        = 2147483647 * 2305843009213693951

3^-1 mod 11     = 4

crt {2 mod 3, 3 mod 5, 2 mod 7} = 23   (mod 105)

phi(36)         = 12
phi(2^31-1)     = 2147483646

(a/11), a=1..10 : 1 ~1 1 1 1 ~1 ~1 ~1 1 ~1

sqrt(10) mod 13 = 6   (6*6 = 36 = 10 mod 13)
primitiveRoot 23 = 5
mu(n), n=1..10  : 1 ~1 ~1 0 ~1 1 ~1 0 0 1
```

## Build & test

Requires [MLton](http://mlton.org/) and/or [Poly/ML](https://polyml.org/).

```sh
make test        # build + run the suite under MLton
make test-poly   # run the suite under Poly/ML (use-and-run)
make all-tests   # both
make example     # build + run examples/demo.sml
make clean
```

## Installing with smlpkg

```sh
smlpkg add github.com/sjqtentacles/sml-numbertheory
smlpkg sync
```

Reference `src/numbertheory.mlb` from your own `.mlb` (it pulls in the
vendored `sml-bigint`), or, under Poly/ML, `use` the sources in dependency
order (`bigint.sig`, `bigint.sml`, `numbertheory.sig`, `numbertheory.sml`) —
see the `test-poly` target in the `Makefile`.

## Layout

Layout B (dependent): own sources live in `src/`; `sml-bigint` is vendored
under `lib/` and loaded first.

```
sml.pkg                                       smlpkg manifest
Makefile                                      MLton + Poly/ML targets
.github/workflows/ci.yml                      CI: MLton + Poly/ML
src/
  numbertheory.sig / numbertheory.sml         the toolkit
  numbertheory.mlb                            public basis (pulls in sml-bigint)
lib/github.com/sjqtentacles/sml-bigint/       vendored dependency
examples/
  demo.sml                                    sieve / factor / crt / phi / Legendre / sqrt / root / mu tour
  sources.mlb
test/
  harness.sml                                 shared assertion harness
  test_numbertheory.sml                       canonical vectors (123 checks)
  entry.sml / main.sml
```

## Tests

123 deterministic checks: the sieve against `pi(100)=25`, `pi(1000)=168`,
`pi(10000)=1229`; factorizations of known composites including the classic
Pollard-rho target `8051 = 83*97`, `10062757 = 1009*9973`, Project Euler's
`600851475143`, and large semiprimes built from the Mersenne primes
`2^19-1`, `2^31-1`, `2^61-1` (so trial division alone cannot crack them, with
each result re-checked as `product = input` and every factor prime); modular
inverses verified against `a*x = 1 (mod m)`; CRT solutions checked against
brute force (including non-coprime and inconsistent systems); Euler's totient
on primes, prime powers, and products; and Jacobi/Legendre vectors — the
worked example `(1001/9907) = ~1`, the `(-1/p)` rule, composite-modulus Jacobi
values, and Legendre cross-checked against Euler's criterion `a^((p-1)/2)`;
Tonelli–Shanks modular square roots (including `sqrt(10) mod 13 = 6`, a
`p ≡ 1 (mod 4)` modulus, the zero root, and non-residues returning `NONE`),
each `SOME r` re-checked as `r*r = n (mod p)` with `r` the canonical smaller
root; the smallest primitive roots `3→2, 5→2, 7→3, 11→2, 13→2, 23→5`; and the
Möbius function on the canonical `μ(1..)` values. Run `make all-tests` to
verify identical output under both compilers.

## License

MIT — see [LICENSE](LICENSE).
