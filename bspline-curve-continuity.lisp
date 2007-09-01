;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Continuity

(defun zero-point (dimension)
  (if (= dimension 2)
      (list 0 0)
      (list 0 0 0)))

(defun point-delta (curve i k max)
  "Kth delta with P_I as the lowest point in the form
[coefficient of P_MAX, remainder of delta]."
  (with-accessors ((p degree)
		   (knots knot-vector)
		   (points control-points))
      curve
      (if (= k 0)
	  (if (= i max)
	      (list 1 (zero-point (bsc-dimension curve)))
	      (list 0 (elt points i)))
	  (let ((p1 (point-delta curve (1+ i) (1- k) max))
		(p2 (point-delta curve i (1- k) max))
		(factor (safe-/ (- p k -1)
				(- (elt knots (+ i p 1))
				   (elt knots (+ i k))))))
	    (list (* (- (first p1) (first p2)) factor)
		  (v* (v- (second p1) (second p2)) factor))))))

(defun bsc-continuous-point (curve k derivative &key from-end)
  "Finds the Kth control point for CURVE so that its Kth derivative
in its first point (or its last point, if FROM-END is T) will be equal to
DERIVATIVE."
  (let* ((n (1- (length (control-points curve))))
	 (delta (point-delta curve (if from-end (- n k) 0)
			     k (if from-end (- n k) k))))
    (v* (v- derivative (second delta)) (/ (first delta)))))

(defun bsc-continuous-curve (curve base-curve continuity &key from-end)
  "Returns a modified version of CURVE, where CONTINUITY with
BASE-CURVE is achieved.
If FROM-END is NIL, the end of BASE-CURVE meets the start of CURVE.
If FROM-END is T, the end of CURVE meets the start of BASE-CURVE."
  (let ((n (1- (length (control-points curve))))
	(new-curve (copy-bspline-curve curve)))
    (do ((k 0 (1+ k)))
	((> k continuity) new-curve)
      (setf (elt (control-points new-curve) (if from-end (- n k) k))
	    (let ((deriv (bsc-evaluate base-curve
				       (if from-end
					   (bsc-lower-parameter base-curve)
					   (bsc-upper-parameter base-curve))
				       :derivative k)))
	      (bsc-continuous-point new-curve k deriv :from-end from-end))))))
