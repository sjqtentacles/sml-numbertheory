(* numbertheory.sml

   A number-theory toolkit on top of the vendored [sml-bigint].  The three
   heavy primitives -- the GCD, modular exponentiation, and Miller-Rabin
   primality -- are taken straight from [BigInt] ([gcd], [modpow],
   [isProbablePrime]); this module never re-implements them.  What it adds:

     - [sieve]      sieve of Eratosthenes over machine ints,
     - [factor]     trial division by small primes, then Pollard's rho,
     - [modInverse] extended Euclidean algorithm,
     - [crt]        Chinese Remainder Theorem (general, non-coprime),
     - [eulerPhi]   totient from the prime factorization, and
     - [jacobi]/[legendre]  the quadratic-residue symbols.

   Everything is pure and deterministic. *)

structure NumberTheory :> NUMBER_THEORY =
struct

  structure B = BigInt

  val zero = B.fromInt 0
  val one  = B.fromInt 1
  val two  = B.fromInt 2

  (* floored remainder/quotient; for a positive modulus the remainder is the
     least non-negative residue *)
  fun bmod (a, m) = #2 (B.divMod (a, m))
  fun bdiv (a, m) = #1 (B.divMod (a, m))
  fun beq (a, b) = B.compare (a, b) = EQUAL
  fun isZero x = beq (x, zero)
  fun isEven x = isZero (bmod (x, two))

  (* ---- sieve of Eratosthenes ---- *)

  fun sieve n =
    if n < 2 then []
    else
      let
        val a = Array.array (n + 1, true)
        val () = Array.update (a, 0, false)
        val () = Array.update (a, 1, false)
        fun mark i =
          if i * i > n then ()
          else
            ( if Array.sub (a, i)
              then
                let
                  fun strike j =
                    if j > n then ()
                    else (Array.update (a, j, false); strike (j + i))
                in strike (i * i) end
              else ()
            ; mark (i + 1) )
        val () = mark 2
        fun collect (i, acc) =
          if i < 2 then acc
          else collect (i - 1, if Array.sub (a, i) then i :: acc else acc)
      in collect (n, []) end

  (* ---- primality: call BigInt's deterministic Miller-Rabin ---- *)

  val mrRounds = B.fromInt 12
  fun isPrime n = B.isProbablePrime (n, mrRounds)

  (* ---- Pollard's rho: a non-trivial divisor of a composite, odd, small-
     factor-free n.  Floyd cycle detection on f(x) = x^2 + c (mod n); retry
     with a fresh c when the run collapses onto n. ---- *)

  fun pollardRho n =
    let
      fun g (x, c) = bmod (B.add (B.mul (x, x), c), n)
      fun tryC c =
        let
          fun loop (x, y, d) =
            if beq (d, one) then
              let
                val x' = g (x, c)
                val y' = g (g (y, c), c)
                val d' = B.gcd (B.abs (B.sub (x', y')), n)
              in loop (x', y', d') end
            else d
          val d = loop (two, two, one)
        in
          if beq (d, n) then tryC (B.add (c, one)) else d
        end
    in tryC one end

  (* insertion sort of a BigInt list (ascending) *)
  fun insert (x, []) = [x]
    | insert (x, y :: ys) =
        if B.compare (x, y) <> GREATER then x :: y :: ys else y :: insert (x, ys)
  fun sortFactors xs = List.foldl insert [] xs

  (* ---- factorization ---- *)

  fun factor n =
    if B.compare (n, one) = LESS then raise Domain
    else if beq (n, one) then []
    else
      let
        (* strip every prime factor below the small-prime bound *)
        fun trial (m, [], acc) = (m, acc)
          | trial (m, p :: ps, acc) =
              let val pb = B.fromInt p
              in
                if B.compare (B.mul (pb, pb), m) = GREATER then (m, acc)
                else
                  let
                    fun strip (m, acc) =
                      if isZero (bmod (m, pb)) then strip (bdiv (m, pb), pb :: acc)
                      else (m, acc)
                    val (m', acc') = strip (m, acc)
                  in trial (m', ps, acc') end
              end
        val (rest, smallFactors) = trial (n, sieve 1000, [])
        (* split whatever remains with rho, testing primality at each step *)
        fun rho (m, acc) =
          if beq (m, one) then acc
          else if isPrime m then m :: acc
          else
            let val d = pollardRho m
            in rho (bdiv (m, d), rho (d, acc)) end
      in
        sortFactors (rho (rest, smallFactors))
      end

  (* ---- extended Euclid: returns (g, x, y) with a*x + b*y = g ---- *)

  fun egcd (a, b) =
    if isZero b then (a, one, zero)
    else
      let
        val (q, r) = B.divMod (a, b)
        val (g, x, y) = egcd (b, r)
      in (g, y, B.sub (x, B.mul (q, y))) end

  fun modInverse (a, m) =
    if isZero m then raise Domain
    else
      let
        val mm = B.abs m
        val (g, x, _) = egcd (bmod (a, mm), mm)
      in
        if beq (g, one) then SOME (bmod (x, mm)) else NONE
      end

  (* ---- Chinese Remainder Theorem (general; non-coprime moduli allowed) ---- *)

  fun crt pairs =
    let
      val () =
        List.app
          (fn (_, m) => if B.compare (m, one) = LESS then raise Domain else ())
          pairs
      (* combine x = r (mod M) with x = a (mod m) *)
      fun merge ((r, M), (a, m)) =
        let
          val g = B.gcd (M, m)
          val diff = B.sub (a, r)
        in
          if not (isZero (bmod (diff, g))) then NONE   (* inconsistent *)
          else
            let
              val lcm = B.mul (bdiv (M, g), m)
              val mg  = bdiv (m, g)
              val coef = bmod (bdiv (M, g), mg)
            in
              case modInverse (coef, mg) of
                  NONE => NONE
                | SOME inv =>
                    let
                      val t = bmod (B.mul (bdiv (diff, g), inv), mg)
                      val x = bmod (B.add (r, B.mul (M, t)), lcm)
                    in SOME (x, lcm) end
            end
        end
      fun loop (acc, []) = SOME acc
        | loop (acc, (a, m) :: rest) =
            (case merge (acc, (bmod (a, m), m)) of
                 NONE => NONE
               | SOME acc' => loop (acc', rest))
    in
      case pairs of
          [] => SOME zero
        | (a0, m0) :: rest =>
            (case loop ((bmod (a0, m0), m0), rest) of
                 NONE => NONE
               | SOME (x, _) => SOME x)
    end

  (* ---- Euler's totient from the factorization ---- *)

  fun eulerPhi n =
    if B.compare (n, one) = LESS then raise Domain
    else if beq (n, one) then one
    else
      let
        fun go ([], _, phi) = phi
          | go (p :: ps, prev, phi) =
              if (case prev of SOME q => beq (p, q) | NONE => false)
              then go (ps, SOME p, B.mul (phi, p))             (* repeated prime *)
              else go (ps, SOME p, B.mul (phi, B.sub (p, one)))(* new prime *)
      in go (factor n, NONE, one) end

  (* ---- Jacobi symbol (and Legendre, which delegates) ---- *)

  val four  = B.fromInt 4
  val eight = B.fromInt 8

  fun jacobi (a, n) =
    if isEven n orelse B.compare (n, one) = LESS then raise Domain
    else
      let
        fun loop (a, n, t) =
          let val a = bmod (a, n)
          in
            if isZero a then (if beq (n, one) then t else 0)
            else
              let
                fun rem2 (a, t) =
                  if isEven a then
                    let
                      val t' =
                        case B.toInt (bmod (n, eight)) of
                            SOME 3 => ~t
                          | SOME 5 => ~t
                          | _ => t
                    in rem2 (bdiv (a, two), t') end
                  else (a, t)
                val (a, t) = rem2 (a, t)
                val t =
                  if B.toInt (bmod (a, four)) = SOME 3
                     andalso B.toInt (bmod (n, four)) = SOME 3
                  then ~t else t
              in loop (n, a, t) end
          end
      in loop (a, n, 1) end

  fun legendre (a, p) = jacobi (a, p)

end
