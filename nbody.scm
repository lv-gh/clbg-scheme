#!r6rs

;; The Computer Language Benchmarks Game
;; http://benchmarksgame.alioth.debian.org/
;;
;; Unoptimized although Scheme'ish and r6rs compatible implementation.
;;
;; Author: Laimonas Vėbra


(import (rnrs base)
        (rnrs lists)
        (rnrs control)
        (rnrs programs) 
        (rnrs io simple)
        (rnrs syntax-case)
        (rnrs records syntactic))

(define *dt* 0.01)
(define *year-days* 365.24)
(define *sun-mass* (* 4 (expt 3.141592653589793 2)))

(define (round* number count)
  ;; Round number to a fraction scale (count of decimal places).
  (let ((scale (expt 10 count)))
    (/ (round (* scale number)) scale)))

(define (print* energy)
  (display (round* energy 9)) (newline))

(define-syntax Λ
  ;; Syntactic sugar; lambda alias.
  (syntax-rules ()
    ([_ <formals> <body> ...] (lambda <formals> <body> ...))))

(define (a•x⃗ a x⃗)
  ;; Multiplies vector by a number
  (vector-map * x⃗ (make-vector (vector-length x⃗) a)))

(define-record-type body
  (fields
   (mutable x⃗ x⃗ x⃗!)       ; position vector (x y z)
   (mutable v⃗ v⃗ v⃗!)       ; velocity vecor (vx vy vz)
   (immutable mass mass)) ; mass

  ;; Custom constructor (mass and velocity values are normalized)
  (protocol (Λ (p)
               (Λ (x⃗ v⃗ mass)
                  (let ([v⃗ (a•x⃗ *year-days* v⃗)]
                        [mass (* mass *sun-mass*)])
                    (p x⃗ v⃗ mass))))))

;; Velocity and position change (advance)
(define (v⃗+ body v⃗ᵥ)
  (v⃗! body (vector-map + (v⃗ body) v⃗ᵥ)))
(define (x⃗+ body x⃗ᵥ)
  (x⃗! body (vector-map + (x⃗ body) x⃗ᵥ)))
 

(define (x⃗•x⃗ x⃗)
  ;; Dot product of a vector
  (apply + (vector->list (vector-map * x⃗ x⃗))))
 
(define (d⃗ bᵢ bⱼ)
  ;; Distance vector '(dx dy dz) between two bodies
  (vector-map - (x⃗ bᵢ) (x⃗ bⱼ)))


(define *sun*
  (make-body '#(0.0 0.0 0.0) ;; x⃗ (x y z)
             '#(0.0 0.0 0.0) ;; v⃗ (vx vy vz)
             1.0))           ;; mass

(define *jupiter*
  (make-body '#(4.84143144246472090 -1.16032004402742839 -1.03622044471123109e-1)
             '#(1.66007664274403694e-3 7.69901118419740425e-3 -6.90460016972063023e-5)
             9.54791938424326609e-4))

(define *saturn*
  (make-body '#(8.34336671824457987 4.12479856412430479 -4.03523417114321381e-1)
             '#(-2.76742510726862411e-3 4.99852801234917238e-3 2.30417297573763929e-5)
             2.85885980666130812e-4))

(define *uranus*
  (make-body '#(1.28943695621391310e1 -1.51111514016986312e1 -2.23307578892655734e-1)
             '#(2.96460137564761618e-03 2.37847173959480950e-03 -2.96589568540237556e-05)
             4.36624404335156298e-05))

(define *neptune*
  (make-body '#(1.53796971148509165e+01 -2.59193146099879641e+01 1.79258772950371181e-01)
             '#(2.68067772490389322e-03 1.62824170038242295e-03 -9.51592254519715870e-05)
             5.15138902046611451e-05))

(define *system* (list *sun* *jupiter* *saturn* *uranus* *neptune*))
;;
;; Distance:
;;   d⃗ᵢⱼ = '(dx dy dz) = '(xᵢ-xⱼ yᵢ-yⱼ zᵢ-zⱼ)
;;   r = sqrt(d⃗•d⃗) = sqrt(dx² + dy² + dz²)
;;
;; Momentum:
;;   p⃗ = mᵢ•v⃗ᵢ, v⃗ᵢ = '(vxᵢ vyᵢ vzᵢ)  
;;
;; Kinetic energy:
;;   Eₖ = 1/2 • mᵢ • (v⃗ᵢ•v⃗ᵢ)
;;
;; Gravitational potencial energy between two bodies:
;;   Uᵢⱼ = -1 • mᵢmⱼ/rᵢⱼ
;;
;; Velocity change (advance step)
;;   dv⃗ᵢ = -1 • a⃗ᵢ [F/mᵢ = (1/mᵢ • d⃗ᵢⱼ (dx dy dz)/r • mᵢmⱼ/r²] • dt
;;   dv⃗ⱼ = -1 • dv⃗ᵢ • mᵢ/mⱼ
;;

(define (offset-momentum)
  (let* ([p⃗ (Λ (body) (a•x⃗ (mass body) (v⃗ body)))]
         [Σp⃗ (Λ (bodies) 
                 (fold-left (Λ (a bᵢ)
                               (vector-map + a (p⃗ bᵢ))) '#(0.0 0.0 0.0) bodies))])
    (v⃗! *sun*  ;; update sun velocity
         (a•x⃗ (/ -1.0 (mass *sun*)) (Σp⃗ *system*)))))


(define (energy)
  (letrec* ([Eₖ (Λ (bᵢ) (* 0.5 (mass bᵢ) (x⃗•x⃗ (v⃗ bᵢ))))]
            [U (Λ (bᵢ bⱼ) (/ (* (mass bᵢ) (mass bⱼ) -1.0)
                             (sqrt (x⃗•x⃗ (d⃗ bᵢ bⱼ)))))]
            [ΣE (Λ (bodies)  ;; kinetic and gravitational potencial energy sum
                   (let ([bᵢ (car bodies)]
                         [b_ (cdr bodies)])
                     (+ (fold-left (Λ (a bⱼ) (+ a (U bᵢ bⱼ))) (Eₖ bᵢ) b_)
                        (if (null? b_) 0.0 (ΣE b_)))))])
           (ΣE *system*)))

(define (advance)
  (letrec ([v⃗← (Λ (bᵢ bⱼ) ;; Velocity advance/update
                   (let* ([d⃗ᵢⱼ (d⃗ bᵢ bⱼ)]
                          [dv⃗ᵢ (a•x⃗ (/ (* -1.0 (mass bⱼ) *dt*)
                                         (expt (x⃗•x⃗ d⃗ᵢⱼ) 3/2)) d⃗ᵢⱼ)])
                     (v⃗+ bᵢ dv⃗ᵢ)
                     (v⃗+ bⱼ (a•x⃗ (/ (mass bᵢ) (mass bⱼ) -1.0) dv⃗ᵢ))))]
           [x⃗← (Λ (bᵢ)    ;; Position advance/update
                   (x⃗+ bᵢ (a•x⃗ *dt* (v⃗ bᵢ))))]  
           [loop (Λ (bodies)
                    (let ([bᵢ (car bodies)]
                          [b_ (cdr bodies)])
                      (for-each (Λ (bⱼ) (v⃗← bᵢ bⱼ)) b_)
                      (x⃗← bᵢ)
                      (if (not (null? b_)) (loop b_))))])
    (loop *system*)))


(let* ([arg (cdr (command-line))]
       [num (if (null? arg) 500000
                (string->number (car arg)))])
  (offset-momentum)
  (print* (energy))
  (do ([i 0 (+ i 1)])
    ((= i num)) (advance))
  (print* (energy)))


