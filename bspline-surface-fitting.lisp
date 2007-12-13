;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Fitting

(defun uniform-parameter-points-2d (points &optional
				    (start-u 0.0) (end-u 1.0)
				    (start-v 0.0) (end-v 1.0))
  "Returns a list of the form (U1 V1 X1 Y1 Z1 U2 V2 X2...),
where U1..UN/V1..VN are equidistant points of the [(U1,V1), (U2,V2)] interval."
  (let ((n (array-dimension points 0))
	(m (array-dimension points 1))
	(len-u (- end-u start-u))
	(len-v (- end-v start-v)))
    (iter (for i from 0 below n)
	  (with ppts)
	  (iter (for j from 0 below m)
		(setf ppts (append ppts
				   (list (+ (/ (* len-u i) (1- n)) start-u)
					 (+ (/ (* len-v j) (1- m)) start-v))
				   (copy-list (aref points i j)))))
	  (finally (return ppts)))))

#+fff
(defun bspline-surface-from-sf (sf)
  "Extracts the fitted b-spline surface from a SF object."
  (let* ((nr-knots-u (sf-get-nr-knots-u sf))
	 (nr-knots-v (sf-get-nr-knots-v sf))
	 (nr-cpts-u (sf-get-nr-ctrl-points-u sf))
	 (nr-cpts-v (sf-get-nr-ctrl-points-v sf))
	 (np (* nr-cpts-u nr-cpts-v 3)))
    (make-bspline-surface (list (sf-get-degree-u sf) (sf-get-degree-v sf))
			  (list
			   (with-foreign-object (knots-u :double nr-knots-u)
			     (sf-get-knot-vector-u sf knots-u)
			     (double-array->vector knots-u nr-knots-u))
			   (with-foreign-object (knots-v :double nr-knots-v)
			     (sf-get-knot-vector-v sf knots-v)
			     (double-array->vector knots-v nr-knots-v)))
			  (coerce
			   (dimensionate
			    (dimensionate
			     (with-foreign-object (points :double np)
			       (sf-get-control-points sf points)
			       (double-array->vector points np))
			     3)
			    nr-cpts-v)
			   'list))))

#+fff
(defun bss-fit (surface points tolerance)
    "Approximates POINTS with a surface, while retaining the
features of SURFACE.
POINTS should be a two-dimensional array."
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let ((resolution (array-dimensions points))
	  (parameter-points (uniform-parameter-points-2d points))
	  (tolerance (* (vlength (apply #'v- (bss-bounding-box surface)))
			tolerance))
	  (sf (sf-create))
	  result)
      (with-double-arrays ((ppoint-array parameter-points)
			   (knot-array-u (first knots))
			   (knot-array-v (second knots)))
	(sf-add-point-group sf (coerce tolerance 'double-float)
			    (* (first resolution) (second resolution))
			    ppoint-array)
	(sf-set-degree-u sf (first degrees))
	(sf-set-degree-v sf (second degrees))
	(sf-set-knot-vector-u sf (length (first knots)) knot-array-u)
	(sf-set-knot-vector-v sf (length (second knots)) knot-array-v)
	(sf-set-smoothness-functional sf :smf-crv)
	(sf-set-optimize-parameters sf nil)
	(unwind-protect
	     (setf result (if (eql (sf-fit sf) :success)
			      (bspline-surface-from-sf sf)
			      nil))
	  (sf-destroy sf)))
      result)))

#-fff
(defun fit (surface points tolerance)
    "Approximates POINTS with a surface, while retaining the
features of SURFACE.
POINTS should be a two-dimensional array."
  'TODO)
