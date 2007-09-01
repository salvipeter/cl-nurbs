;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Continuity

(defun bss-continuous-surface (surface base-surface continuity
			   &key (u-direction t) from-end)
  "Returns a modified version of SURFACE, where CONTINUITY with
BASE-SURFACE is achieved in the u direction (or v, if U-DIRECTION is NIL).
If FROM-END is NIL, the end of BASE-SURFACE meets the start of SURFACE.
If FROM-END is T, the end of SURFACE meets the start of BASE-SURFACE."
  (let ((new-surface (copy-bspline-surface surface)))
    (with-accessors ((degrees degrees)
		     (knots knot-vectors)
		     (net control-net))
	new-surface
      (iter (with dim = (array-dimensions net))
	    (for i from 0 below (if u-direction (second dim) (first dim)))
	    (for base = (bss-construction-curve base-surface i
						:u-direction u-direction))
	    (for original = (bss-construction-curve surface i
						    :u-direction u-direction))
	    (for new = (continuous-curve original base continuity
					 :from-end from-end))
	    (if u-direction
		(array-set-column net i (control-points new))
		(array-set-row net i (control-points new)))))
    new-surface))
