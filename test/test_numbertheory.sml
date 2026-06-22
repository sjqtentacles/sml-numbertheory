(* test_numbertheory.sml

   Deterministic test suite for sml-numbertheory.  Reference values are the
   canonical ones from the literature: prime-counting values pi(100)=25 and
   pi(1000)=168, factorizations of known composites (including the classic
   Pollard-rho target 8051 = 83*97 and large semiprimes built from Mersenne
   primes so trial division alone cannot crack them), modular inverses checked
   against a*x = 1 (mod m), CRT solutions checked against brute force, Euler's
   totient on primes / prime powers / products, and the worked Jacobi example
   (1001/9907) = ~1 alongside Legendre symbols cross-checked with Euler's
   criterion.  Output is plain text and identical across MLton and Poly/ML. *)

structure NumberTheoryTests =
struct

  structure NT = NumberTheory
  structure B  = BigInt

  (* ---- helpers ---- *)

  fun b n  = B.fromInt n
  fun bs s = valOf (B.fromString s)
  val str  = B.toString
  val zero = b 0
  val one  = b 1
  fun beq (x, y) = B.compare (x, y) = EQUAL
  fun bmod (a, m) = #2 (B.divMod (a, m))

  (* product of a BigInt list (empty product = 1) *)
  fun prod xs = List.foldl (fn (x, acc) => B.mul (x, acc)) one xs

  (* non-decreasing? *)
  fun nondec [] = true
    | nondec [_] = true
    | nondec (x :: y :: r) = B.compare (x, y) <> GREATER andalso nondec (y :: r)

  fun allPrime xs = List.all (fn p => B.isProbablePrime (p, b 12)) xs

  (* a factorization is correct iff the product is the input, every factor is
     prime, and the list is sorted *)
  fun checkFactor name n =
    let val fs = NT.factor n
    in Harness.check (name ^ " product==n") (beq (prod fs, n))
     ; Harness.check (name ^ " all-prime")  (allPrime fs)
     ; Harness.check (name ^ " sorted")     (nondec fs)
    end

  fun checkFactorList name (expected, n) =
    Harness.checkStringList name (List.map str expected, List.map str (NT.factor n))

  (* modInverse: assert a*x = 1 (mod m) and 0 <= x < |m| *)
  fun checkInvMul name (a, m) =
    case NT.modInverse (a, m) of
        SOME x =>
          Harness.check name
            (beq (bmod (B.mul (a, x), m), bmod (one, m))
             andalso B.compare (x, zero) <> LESS
             andalso B.compare (x, B.abs m) = LESS)
      | NONE => Harness.check name false

  fun checkInvNone name (a, m) =
    Harness.check name (case NT.modInverse (a, m) of NONE => true | _ => false)

  fun checkInvVal name (expected, (a, m)) =
    Harness.checkString name
      (str expected, case NT.modInverse (a, m) of SOME x => str x | NONE => "NONE")

  (* crt over small int systems, rendered as strings *)
  fun crtStr sys =
    case NT.crt (List.map (fn (a, m) => (b a, b m)) sys) of
        SOME x => B.toString x
      | NONE => "NONE"

  fun checkCrt name (expected, sys) =
    Harness.checkString name (expected, crtStr sys)

  (* brute-force smallest non-negative solution in [0, product) *)
  fun bruteCrt sys =
    let
      val p = List.foldl (fn ((_, m), acc) => acc * m) 1 sys
      fun ok x = List.all (fn (a, m) => x mod m = a mod m) sys
      fun find x = if x >= p then NONE else if ok x then SOME x else find (x + 1)
    in find 0 end

  fun checkCrtBrute name sys =
    Harness.checkString name
      ((case bruteCrt sys of SOME x => Int.toString x | NONE => "NONE"), crtStr sys)

  fun checkPhiI name (expected, n) =
    Harness.checkString name (str (b expected), str (NT.eulerPhi (b n)))

  fun checkPhiB name (expected, n) =
    Harness.checkString name (str expected, str (NT.eulerPhi n))

  fun checkJac name (expected, (a, n)) =
    Harness.checkInt name (expected, NT.jacobi (b a, b n))

  fun checkLeg name (expected, (a, p)) =
    Harness.checkInt name (expected, NT.legendre (b a, b p))

  (* independent Legendre via Euler's criterion a^((p-1)/2) mod p *)
  fun legByEuler (a, p) =
    let
      val pp = b p
      val aa = bmod (b a, pp)
    in
      if beq (aa, zero) then 0
      else
        let val r = B.modpow (aa, #1 (B.divMod (B.sub (pp, one), b 2)), pp)
        in if beq (r, one) then 1 else ~1 end
    end

  fun checkLegEuler name (a, p) =
    Harness.checkInt name (legByEuler (a, p), NT.legendre (b a, b p))

  (* ---- sieve ---- *)

  fun sieveTests () =
    ( Harness.section "sieve (Eratosthenes)"
    ; Harness.checkIntList "sieve 1 = []"  ([], NT.sieve 1)
    ; Harness.checkIntList "sieve 0 = []"  ([], NT.sieve 0)
    ; Harness.checkIntList "sieve ~5 = []" ([], NT.sieve ~5)
    ; Harness.checkIntList "sieve 2 = [2]" ([2], NT.sieve 2)
    ; Harness.checkIntList "sieve 10"      ([2,3,5,7], NT.sieve 10)
    ; Harness.checkIntList "sieve 13 (inclusive endpoint)"
        ([2,3,5,7,11,13], NT.sieve 13)
    ; Harness.checkIntList "sieve 30"
        ([2,3,5,7,11,13,17,19,23,29], NT.sieve 30)
    ; Harness.checkInt "pi(100) = 25"  (25,  List.length (NT.sieve 100))
    ; Harness.checkInt "pi(1000) = 168" (168, List.length (NT.sieve 1000))
    ; Harness.checkInt "pi(10000) = 1229" (1229, List.length (NT.sieve 10000)) )

  (* ---- factor ---- *)

  fun factorTests () =
    let
      (* large semiprimes with no small factors -> forces Pollard's rho *)
      val m19 = bs "524287"               (* 2^19 - 1, Mersenne prime *)
      val m31 = bs "2147483647"           (* 2^31 - 1, Mersenne prime *)
      val m61 = bs "2305843009213693951"  (* 2^61 - 1, Mersenne prime *)
    in
      Harness.section "factor (trial division + Pollard rho)"
    ; Harness.checkStringList "factor 1 = []" ([], List.map str (NT.factor one))
    ; checkFactorList "factor 2 = [2]"      ([b 2], b 2)
    ; checkFactorList "factor 97 (prime)"   ([b 97], b 97)
    ; checkFactorList "factor 12 = 2,2,3"   ([b 2, b 2, b 3], b 12)
    ; checkFactorList "factor 360"          ([b 2,b 2,b 2,b 3,b 3,b 5], b 360)
    ; checkFactorList "factor 8051 = 83,97" ([b 83, b 97], b 8051)
    ; checkFactorList "factor 10062757 = 1009,9973"
        ([b 1009, b 9973], b 10062757)
    ; checkFactor "factor 2^31-1 (prime)" m31
    ; checkFactorList "factor (2^19-1)(2^31-1)"
        ([m19, m31], B.mul (m19, m31))
    ; checkFactorList "factor (2^31-1)(2^61-1) [rho]"
        ([m31, m61], B.mul (m31, m61))
    ; checkFactor "factor 600851475143" (bs "600851475143")   (* Project Euler 3 *)
    ; Harness.checkRaises "factor 0 raises" (fn () => NT.factor zero)
    ; Harness.checkRaises "factor ~6 raises" (fn () => NT.factor (b ~6))
    end

  (* ---- modInverse ---- *)

  fun modInverseTests () =
    ( Harness.section "modInverse (extended Euclid)"
    ; checkInvVal "3^-1 mod 11 = 4"   (b 4,  (b 3, b 11))
    ; checkInvVal "10^-1 mod 17 = 12" (b 12, (b 10, b 17))
    ; checkInvVal "7^-1 mod 26 = 15"  (b 15, (b 7, b 26))
    ; checkInvVal "(~3)^-1 mod 11 = 7" (b 7, (b ~3, b 11))
    ; checkInvVal "1^-1 mod 1 = 0"    (b 0,  (b 1, b 1))
    ; checkInvNone "2^-1 mod 4 = NONE" (b 2, b 4)
    ; checkInvNone "6^-1 mod 9 = NONE" (b 6, b 9)
    ; checkInvMul "123456789^-1 mod 1000000007"
        (b 123456789, b 1000000007)
    ; checkInvMul "big a^-1 mod prime"
        (bs "987654321987654321", bs "2305843009213693951")
    ; Harness.checkRaises "modInverse (_,0) raises"
        (fn () => NT.modInverse (b 3, zero)) )

  (* ---- crt ---- *)

  fun crtTests () =
    ( Harness.section "crt (Chinese Remainder Theorem)"
    ; checkCrt "empty system = 0"        ("0", [])
    ; checkCrt "single (7 mod 5) = 2"    ("2", [(7, 5)])
    ; checkCrt "classic 2,3,2 / 3,5,7"   ("23", [(2,3),(3,5),(2,7)])
    ; checkCrt "0,3,4 / 3,4,5"           ("39", [(0,3),(3,4),(4,5)])
    ; checkCrt "non-coprime consistent"  ("8",  [(2,6),(8,12)])
    ; checkCrt "non-coprime inconsistent = NONE" ("NONE", [(1,2),(2,4)])
    ; checkCrtBrute "brute 1,2 / 3,5"            [(1,3),(2,5)]
    ; checkCrtBrute "brute 2,3,5 / 5,7,11"       [(2,5),(3,7),(5,11)]
    ; checkCrtBrute "brute 4 congruences"        [(1,2),(2,3),(3,5),(4,7)]
    ; Harness.checkRaises "crt with modulus 0 raises"
        (fn () => NT.crt [(b 1, zero)]) )

  (* ---- eulerPhi ---- *)

  fun eulerPhiTests () =
    let val m31 = bs "2147483647"
    in
      Harness.section "eulerPhi (Euler's totient)"
    ; checkPhiI "phi(1) = 1"   (1, 1)
    ; checkPhiI "phi(7) = 6"   (6, 7)
    ; checkPhiI "phi(13) = 12" (12, 13)
    ; checkPhiI "phi(9) = 6 (3^2)"   (6, 9)
    ; checkPhiI "phi(8) = 4 (2^3)"   (4, 8)
    ; checkPhiI "phi(16) = 8"        (8, 16)
    ; checkPhiI "phi(25) = 20 (5^2)" (20, 25)
    ; checkPhiI "phi(27) = 18 (3^3)" (18, 27)
    ; checkPhiI "phi(10) = 4"  (4, 10)
    ; checkPhiI "phi(12) = 4"  (4, 12)
    ; checkPhiI "phi(15) = 8"  (8, 15)
    ; checkPhiI "phi(36) = 12" (12, 36)
    ; checkPhiI "phi(8051) = 7872 (83*97)" (7872, 8051)
    ; checkPhiB "phi(2^31-1) = 2^31-2" (B.sub (m31, one), m31)
    ; Harness.checkRaises "phi 0 raises" (fn () => NT.eulerPhi zero)
    end

  (* ---- jacobi / legendre ---- *)

  fun jacobiTests () =
    ( Harness.section "jacobi / legendre symbols"
    (* Legendre mod 7: QRs are {1,2,4} *)
    ; checkLeg "(1/7) = 1"  (1,  (1, 7))
    ; checkLeg "(2/7) = 1"  (1,  (2, 7))
    ; checkLeg "(3/7) = ~1" (~1, (3, 7))
    ; checkLeg "(4/7) = 1"  (1,  (4, 7))
    ; checkLeg "(5/7) = ~1" (~1, (5, 7))
    ; checkLeg "(6/7) = ~1" (~1, (6, 7))
    ; checkLeg "(7/7) = 0"  (0,  (7, 7))
    (* Legendre mod 11: QRs are {1,3,4,5,9} *)
    ; checkLeg "(2/11) = ~1" (~1, (2, 11))
    ; checkLeg "(3/11) = 1"  (1,  (3, 11))
    ; checkLeg "(5/11) = 1"  (1,  (5, 11))
    ; checkLeg "(7/11) = ~1" (~1, (7, 11))
    (* (-1/p) = (-1)^((p-1)/2) *)
    ; checkJac "(~1/5) = 1"   (1,  (~1, 5))
    ; checkJac "(~1/7) = ~1"  (~1, (~1, 7))
    ; checkJac "(~1/11) = ~1" (~1, (~1, 11))
    ; checkJac "(~1/13) = 1"  (1,  (~1, 13))
    (* Jacobi with composite lower modulus *)
    ; checkJac "(1/1) = 1"   (1, (1, 1))
    ; checkJac "(0/1) = 1"   (1, (0, 1))
    ; checkJac "(5/1) = 1"   (1, (5, 1))
    ; checkJac "(0/15) = 0"  (0, (0, 15))
    ; checkJac "(2/15) = 1"  (1, (2, 15))
    ; checkJac "(6/15) = 0 (gcd 3)" (0, (6, 15))
    ; checkJac "(19/45) = 1" (1, (19, 45))
    (* the classic worked example from textbooks *)
    ; checkJac "(1001/9907) = ~1" (~1, (1001, 9907))
    (* cross-check Legendre against Euler's criterion *)
    ; checkLegEuler "(a/13) via Euler, a=2" (2, 13)
    ; checkLegEuler "(a/13) via Euler, a=7" (7, 13)
    ; checkLegEuler "(a/97) via Euler, a=11" (11, 97)
    ; checkLegEuler "(a/97) via Euler, a=50" (50, 97)
    ; checkLegEuler "(a/101) via Euler, a=2" (2, 101)
    ; Harness.checkRaises "jacobi even modulus raises"
        (fn () => NT.jacobi (b 3, b 8)) )

  fun run () =
    ( sieveTests ()
    ; factorTests ()
    ; modInverseTests ()
    ; crtTests ()
    ; eulerPhiTests ()
    ; jacobiTests () )

end
