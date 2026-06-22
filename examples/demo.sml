(* demo.sml

   A small command-line tour of sml-numbertheory: count primes with the sieve,
   factor an RSA-style semiprime built from two large primes (Pollard's rho),
   compute a modular inverse, solve a system of congruences with the CRT,
   evaluate Euler's totient, and print a row of Legendre symbols.  Everything
   is exact and arbitrary-precision via the vendored sml-bigint. *)

structure NT = NumberTheory
structure B  = BigInt

fun line s = print (s ^ "\n")
fun bs s = valOf (B.fromString s)

val () = line "sml-numbertheory demo"
val () = line "====================="
val () = line ""

(* --- sieve --- *)
val () = line ("primes up to 30 : " ^
               String.concatWith " " (List.map Int.toString (NT.sieve 30)))
val () = line ("pi(1000)        = " ^ Int.toString (List.length (NT.sieve 1000)))
val () = line ""

(* --- factor a semiprime made of two Mersenne primes --- *)
val p = bs "2147483647"          (* 2^31 - 1 *)
val q = bs "2305843009213693951" (* 2^61 - 1 *)
val n = B.mul (p, q)
val () = line ("n = (2^31-1)(2^61-1) = " ^ B.toString n)
val () = line ("factor n        = " ^
               String.concatWith " * " (List.map B.toString (NT.factor n)))
val () = line ""

(* --- modular inverse --- *)
val () = line ("3^-1 mod 11     = " ^
               (case NT.modInverse (B.fromInt 3, B.fromInt 11) of
                    SOME x => B.toString x | NONE => "none"))
val () = line ""

(* --- CRT: x = 2 (mod 3), 3 (mod 5), 2 (mod 7) --- *)
val sys = [ (B.fromInt 2, B.fromInt 3)
          , (B.fromInt 3, B.fromInt 5)
          , (B.fromInt 2, B.fromInt 7) ]
val () = line ("crt {2 mod 3, 3 mod 5, 2 mod 7} = " ^
               (case NT.crt sys of SOME x => B.toString x | NONE => "none")
               ^ "   (mod 105)")
val () = line ""

(* --- Euler's totient --- *)
val () = line ("phi(36)         = " ^ B.toString (NT.eulerPhi (B.fromInt 36)))
val () = line ("phi(2^31-1)     = " ^ B.toString (NT.eulerPhi p))
val () = line ""

(* --- Legendre symbols (a / 11) for a = 1..10 --- *)
fun legRow () =
  let
    fun one a = Int.toString (NT.legendre (B.fromInt a, B.fromInt 11))
  in
    String.concatWith " " (List.map one [1,2,3,4,5,6,7,8,9,10])
  end
val () = line ("(a/11), a=1..10 : " ^ legRow ())
val () = line ""

(* --- modular square root via Tonelli-Shanks --- *)
val () = line ("sqrt(10) mod 13 = " ^
               (case NT.tonelliShanks (B.fromInt 10, B.fromInt 13) of
                    SOME r => B.toString r | NONE => "none")
               ^ "   (6*6 = 36 = 10 mod 13)")

(* --- smallest primitive root mod 23 --- *)
val () = line ("primitiveRoot 23 = " ^
               (case NT.primitiveRoot (B.fromInt 23) of
                    SOME g => B.toString g | NONE => "none"))

(* --- Moebius function mu(n) for n = 1..10 --- *)
fun muRow () =
  let fun one n = Int.toString (NT.moebius (B.fromInt n))
  in String.concatWith " " (List.map one [1,2,3,4,5,6,7,8,9,10]) end
val () = line ("mu(n), n=1..10  : " ^ muRow ())
