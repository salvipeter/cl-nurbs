;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Affine operators

(defun safe-/ (a b)
  "Returns 0.0 when dividing by 0.0."
  (if (= b 0.0)
      0.0
      (/ a b)))

(defun 2d->3d (lst)
  (if (= (length lst) 2)
      (append lst '(0.0))
      lst))

(defun vlength (u)
  (sqrt (apply #'+ (mapcar #'(lambda (x) (* x x)) u))))

(defun v+ (&rest u)
  (apply #'mapcar #'+ u))

(defun v- (&rest u)
  (apply #'mapcar #'- u))

(defun v* (u &rest scale)
  (mapcar #'(lambda (x) (apply #'* x scale)) u))

(defun vnormalize (u)
  (v* u (safe-/ 1.0 (vlength u))))

(defun scalar-product (u v)
  (apply #'+ (mapcar #'* u v)))

(defun cross-product (u v)
  "Cross product for 3D vectors."
  (list (- (* (second u) (third v)) (* (third u) (second v)))
	(- (* (third u) (first v)) (* (first u) (third v)))
	(- (* (first u) (second v)) (* (second u) (first v)))))

(defun point-distance (a b)
  (vlength (v- b a)))

(defun affine-combine (u i v)
  "Linear interpolation that returns U when I = 0 and V when I = 1."
  (mapcar #'(lambda (x y) (+ (* x (- 1 i)) (* y i))) u v))

(defun interpolate (x i y)
  "Linear interpolation that returns X when I = 0 and Y when I = 1."
  (+ (* x (- 1 i)) (* y i)))
