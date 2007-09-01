;;; -*- mode: lisp; syntax: common-lisp -*-
;;;
;;; Downhill Simplex method
;;;
;;; Peter Salvi, 2007
;;;
;;; Based on the one in `Numerical Recipes in C'

(defpackage :downhill-simplex
  (:use :common-lisp)
  (:export :minimize))

(in-package :downhill-simplex)

(defun make-starting-simplex (point)
  "Returns a list of N+1 N-dimensional points, first of which is a copy of
POINT itself, then its translations by all base vectors."
  (let ((n (length point)))
    (append (list (copy-list point))
	    (loop for i from 0 below n collect
		  (let ((new-point (copy-list point)))
		    (setf (nth i new-point) (+ (nth i new-point) 1))
		    new-point)))))

(defun high-next-low (lst)
  "Returns a list with the indices of the highest, next-to-highest
and lowest values, respectively."
  (let ((sorted (sort (copy-list lst) #'>)))
    (list (position (first sorted) lst)
	  (position (second sorted) lst)
	  (position (car (last sorted)) lst))))

(defun scale-vertex (simplex vertex factor)
  "Scale VERTEX in SIMPLEX by FACTOR. Returns a point."
  (let ((high-point (nth vertex simplex))
	(point-sum (apply #'mapcar #'+ simplex))
	(alpha (/ (- 1 factor) (length (first simplex)))))
    (let ((sum-high (mapcar #'- point-sum high-point)))
      (mapcar #'(lambda (x y) (+ (* x alpha) (* y factor)))
	      sum-high high-point))))

(defun scale-simplex (simplex vertex factor)
  "Scale the whole SIMPLEX towards VERTEX by FACTOR. Returns a simplex."
  (let ((low-point (nth vertex simplex)))
    (mapcar #'(lambda (point)
		(mapcar #'(lambda (x y) (* (+ x y) factor)) point low-point))
	    simplex)))

(defun minimize (function start iteration)
  "Execute ITERATION steps of downhill simplex method, minimizing FUNCTION.
Starting simplex consists of the START point and its translations in every
direction."
  (let* ((simplex (make-starting-simplex start))
	 (values (mapcar function simplex)))
    (dotimes (i iteration)
      (let* ((hnl (high-next-low values))
	     (reflection (scale-vertex simplex (first hnl) -1.0))
	     (reflection-value (funcall function reflection)))
	(cond ((<= reflection-value (nth (third hnl) values))
	       (let* ((extension (scale-vertex simplex (first hnl) -2.0))
		      (extension-value (funcall function extension)))
		 (if (<= extension-value reflection-value)
		     (setf (nth (first hnl) simplex) extension
			   (nth (first hnl) values) extension-value)
		     (setf (nth (first hnl) simplex) reflection
			   (nth (first hnl) values) reflection-value))))
	      ((>= reflection-value (nth (second hnl) values))
	       (let* ((contraction (scale-vertex simplex (first hnl) 0.5))
		      (contraction-value (funcall function contraction)))
		 (if (>= contraction-value
			 (min reflection-value (nth (first hnl) values)))
		     (setf simplex (scale-simplex simplex (third hnl) 0.5)
			   values (mapcar function simplex))
		     (setf (nth (first hnl) simplex) contraction
			   (nth (first hnl) values) contraction-value))))
	      (t (setf (nth (first hnl) simplex) reflection
		       (nth (first hnl) values) reflection-value)))))
    (nth (third (high-next-low values)) simplex)))
