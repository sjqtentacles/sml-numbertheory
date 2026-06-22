(* numbertheory.sig

   A small number-theory toolkit in pure Standard ML, layered on top of the
   vendored [sml-bigint].  This module deliberately does NOT re-implement the
   primitives BigInt already provides; instead it CALLS them:

     - [BigInt.gcd]              the binary (Stein) GCD,
     - [BigInt.modpow]          modular exponentiation, and
     - [BigInt.isProbablePrime] deterministic Miller-Rabin,

   and builds the higher-level operations -- the sieve of Eratosthenes,
   integer factorization (trial division + Pollard's rho), modular inverse
   (extended Euclid), the Chinese Remainder Theorem, Euler's totient, and the
   Jacobi/Legendre symbols -- on top of them.

   Every value of type [BigInt.int] below is the arbitrary-precision integer;
   the host machine integer ([Int.int]) appears only in [sieve], whose input
   and output are small by nature.  All functions are pure and deterministic:
   the same arguments always produce the same result, byte-identically under
   MLton and Poly/ML. *)

signature NUMBER_THEORY =
sig
  (* [sieve n] is the ascending list of all primes p with 2 <= p <= n,
     computed with the sieve of Eratosthenes.  Empty for n < 2. *)
  val sieve : int -> int list

  (* [factor n] is the multiset of prime factors of n (n >= 1) in
     non-decreasing order, so that the product of the list equals n and every
     element is prime.  [factor 1] is the empty list.  Small factors are
     stripped by trial division and the remaining cofactor is split with
     Pollard's rho, primality being decided by [BigInt.isProbablePrime].
     Raises [Domain] for n <= 0. *)
  val factor : BigInt.int -> BigInt.int list

  (* [modInverse (a, m)] is SOME x with 0 <= x < |m| and a*x = 1 (mod m),
     computed by the extended Euclidean algorithm, or NONE when gcd(a,m) <> 1
     (no inverse exists).  Raises [Domain] when m = 0. *)
  val modInverse : BigInt.int * BigInt.int -> BigInt.int option

  (* [crt [(a1,m1), ...]] solves the system x = ai (mod mi) by the Chinese
     Remainder Theorem.  On success it returns SOME x, the unique solution in
     [0, M) where M = lcm-style product of the moduli; pairwise-coprimality is
     not required, but NONE is returned when the congruences are inconsistent.
     The empty system returns SOME 0.  Each modulus must be positive
     (raises [Domain] otherwise). *)
  val crt : (BigInt.int * BigInt.int) list -> BigInt.int option

  (* [eulerPhi n] is Euler's totient phi(n) for n >= 1: the count of integers
     in [1, n] coprime to n.  [eulerPhi 1] is 1.  Computed from the prime
     factorization.  Raises [Domain] for n <= 0. *)
  val eulerPhi : BigInt.int -> BigInt.int

  (* [jacobi (a, n)] is the Jacobi symbol (a/n) for odd n >= 1, taking the
     values ~1, 0 or 1.  For n = 1 the symbol is 1; it is 0 exactly when
     gcd(a,n) <> 1.  Raises [Domain] when n is even or n < 1. *)
  val jacobi : BigInt.int * BigInt.int -> int

  (* [legendre (a, p)] is the Legendre symbol (a/p) for an odd prime p: 0 when
     p | a, 1 when a is a non-zero quadratic residue mod p, and ~1 otherwise.
     It coincides with the Jacobi symbol on a prime modulus and is computed as
     such; the caller is responsible for p being an odd prime (the value is
     only meaningful then).  Raises [Domain] when p is even or p < 1. *)
  val legendre : BigInt.int * BigInt.int -> int

  (* [tonelliShanks (n, p)] is SOME r with 0 <= r < p and r*r = n (mod p) -- a
     modular square root of n mod the odd prime p -- or NONE when n is a
     quadratic non-residue mod p.  When n = 0 (mod p) the root is 0.  Of the
     two roots r and p-r the smaller is returned, so the result is
     deterministic.  The caller is responsible for p being an odd prime (the
     value is only meaningful then); it is computed by the Tonelli-Shanks
     algorithm. *)
  val tonelliShanks : BigInt.int * BigInt.int -> BigInt.int option

  (* [primitiveRoot p] is SOME g, the smallest primitive root modulo p -- a
     generator of the multiplicative group of units mod p -- or NONE when no
     primitive root exists.  It is intended for a prime modulus p, where a
     primitive root always exists; more generally it returns the least g
     coprime to p whose multiplicative order equals eulerPhi p, scanning
     g = 1, 2, ... and so reporting NONE only after exhausting [1, p).  Raises
     [Domain] for p < 1. *)
  val primitiveRoot : BigInt.int -> BigInt.int option

  (* [moebius n] is the Moebius function mu(n) for n >= 1: 1 when n is
     square-free with an even number of prime factors, ~1 when square-free with
     an odd number, and 0 when n is divisible by a square > 1.  [moebius 1] is
     1.  Computed from the prime factorization.  Raises [Domain] for n <= 0. *)
  val moebius : BigInt.int -> int
end
