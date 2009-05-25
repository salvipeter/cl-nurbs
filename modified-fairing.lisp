;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

(defun average-double-array (matrix iterations)
  "Averages the contents destructively of MATRIX preserving border elements."
  (flet ((da-aref (da i j) (elt (elt da i) j)))
    (iter (with n = (length matrix))
	  (with m = (length (elt matrix 0)))
	  (with tmp = (make-array (list n m)))
	  (repeat iterations)
	  (iter (for i from 1 below (1- n))
		(iter (for j from 1 below (1- m))
		      (setf (aref tmp i j)
			    (/ (+ (da-aref matrix (1- i) j)
				  (da-aref matrix (1+ i) j)
				  (da-aref matrix i (1- j))
				  (da-aref matrix i (1+ j)))
			       4))))
	  (iter (for i from 1 below (1- n))
		(iter (for j from 1 below (1- m))
		      (setf (elt (elt matrix i) j) (aref tmp i j)))))))

;;; Modified version of FAIR-IN-ONE-DIRECTION, where the target curvatures
;;; are computed globally for the whole surface.
(defun fair-in-one-direction (surface resolution iteration distance
			      &key (u-direction t))
  (flet ((ufirst (lst) (if u-direction (first lst) (second lst)))
	 (usecond (lst) (if u-direction (second lst) (first lst)))
	 (uswap (lst) (if u-direction (reverse lst) lst))
	 (make-double-array (sizes)
	   (let ((result (make-array (first sizes))))
	     (dotimes (i (first sizes))
	       (setf (elt result i) (make-array (second sizes))))
	     result)))
    (let* ((parameters (make-double-array (uswap resolution)))
	   (curvatures (make-double-array (uswap resolution)))
	   (curves (make-array (usecond resolution)))
	   (max (1- (usecond resolution))))
      (iter (with lower = (bss-lower-parameter surface))
	    (with upper = (bss-upper-parameter surface))
	    (for i from 0 to max)
	    (for v = (usecond (affine-combine lower (/ i max) upper)))
	    (for curve = (bss-get-surface-curve surface v :u-curve u-direction))
	    (setf (elt curves i) curve)
	    (iter (for u in-vector
		       (arc-length-sampling curve
					    (bsc-lower-parameter curve)
					    (bsc-upper-parameter curve)
					    (ufirst resolution)))
		  (for j upfrom 0)
		  (setf (elt (elt parameters i) j) u
			(elt (elt curvatures i) j) (bsc-curvature curve u))))
      ;;; TODO: set the side curvatures here
      (average-double-array curvatures iteration)
      (iter (for i from 0 to max)
	    (for left = (integrate (elt curves i) (elt parameters i)
				   (elt curvatures i) distance
				   :from-right nil))
	    (for right = (integrate (elt curves i) (elt parameters i)
				    (elt curvatures i) distance
				    :from-right t))
	    (collect (blend-points left (nreverse right)))))))

;;; Test

;; (defparameter *bottom*
;;   (first (read-rbn "/home/salvi/project/cl-nurbs/models/bottom.rbn")))
;; (defparameter *mesh*
;;   (bss-faired-mesh *bottom* '(100 100) 100 0.1d0))
;; (write-points2-pts *mesh* "/tmp/faired-mesh.pts")

;;; Az eredmeny jobb, mint az eredeti algoritmus.
;;; Azonban kerdeses, hogy jo-e, hogy az ivhossz szerint mintavetelezett
;;; gorbepontokat matrixkent kezeljuk (olyan nagyon rossz mondjuk nem lehet).
