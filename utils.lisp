;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

(defmacro define-constant (name value &optional doc)
  `(defconstant ,name (if (boundp ',name) (symbol-value ',name) ,value)
     ,@(when doc (list doc))))


;;; 2D array utilities

(defun array-get-row (array n)
  (iter (for i from 0 below (array-dimension array 1))
	(collect (aref array n i))))

(defun array-get-column (array n)
  (iter (for i from 0 below (array-dimension array 0))
	(collect (aref array i n))))

(defun array-set-row (array n row)
  (iter (for i from 0 below (array-dimension array 1))
	(setf (aref array n i) (elt row i))))

(defun array-set-column (array n column)
  (iter (for i from 0 below (array-dimension array 0))
	(setf (aref array i n) (elt column i))))

(defun subarray-2d (array min1 min2 max1 max2)
  "Returns a copy of the subarray of the given 2-dimensional array.
The array should contain lists. NIL arguments mean extreme values."
  (let ((min1 (or min1 0))
	(min2 (or min2 0))
	(max1 (or max1 (array-dimension array 0)))
	(max2 (or max2 (array-dimension array 1))))
    (let ((result (make-array (list (- max1 min1) (- max2 min2)))))
      (iter (for i from min1 below max1)
	    (iter (for j from min2 below max2)
		  (setf (aref result (- i min1) (- j min2))
			(copy-list (aref array i j)))))
      result)))
