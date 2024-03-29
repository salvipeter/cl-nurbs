;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Fairing

(defun bss-knot-removal-reinsertion (surface uv)
  "Removes and reinserts a knot by moving the specified control point."
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (when (not (= (first degrees) (second degrees) 3))
      (error "this function can only be used with 3rd-degree surfaces"))
    (let* ((u (first uv))
	   (v (second uv))
	   (i (+ u 2))
	   (j (+ v 2))
	   (ku (first knots))
	   (kv (second knots))
	   (lu (affine-combine (aref net (- u 2) v)
			       (safe-/ (- (elt ku (1+ i)) (elt ku (- i 3)))
				       (- (elt ku i) (elt ku (- i 3))))
			       (aref net (- u 1) v)))
	   (ru (affine-combine (aref net (+ u 2) v)
			       (safe-/ (- (elt ku (+ i 3)) (elt ku (1- i)))
				       (- (elt ku (+ i 3)) (elt ku i)))
			       (aref net (+ u 1) v)))
	   (krr-u (affine-combine ru
				  (safe-/ (- (elt ku (+ i 2)) (elt ku i))
					  (- (elt ku (+ i 2))
					     (elt ku (- i 2))))
				  lu))
	   (lv (affine-combine (aref net u (- v 2))
			       (safe-/ (- (elt kv (1+ j)) (elt kv (- j 3)))
				       (- (elt kv j) (elt kv (- j 3))))
			       (aref net u (- v 1))))
	   (rv (affine-combine (aref net u (+ v 2))
			       (safe-/ (- (elt kv (+ j 3)) (elt kv (1- j)))
				       (- (elt kv (+ j 3)) (elt kv j)))
			       (aref net u (+ v 1))))
	   (krr-v (affine-combine rv
				  (safe-/ (- (elt kv (+ j 2)) (elt kv j))
					  (- (elt kv (+ j 2))
					     (elt kv (- j 2))))
				  lv)))
      (affine-combine krr-u 0.5 krr-v))))

(defun get-curve-at-cp (surface i &key (u-curve t))
  "Returns the u (v) curve assembled from the Ith control points
in the v (u) direction. For internal use only, reuses data of the surface."
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let* ((1st (if u-curve 0 1)) (2nd (if u-curve 1 0))
 	   (length (nth 1st (array-dimensions net))))
      (make-bspline-curve (nth 1st degrees)
			  (nth 1st knots)
			  (iter (for j below length)
				(collect (aref net
					       (nth 1st (list j i))
					       (nth 2nd (list j i)))))))))

(defun bss-construction-curve (surface i &key (u-direction t))
  "Returns the Ith construction curve of SURFACE."
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (make-bspline-curve (if u-direction (first degrees) (second degrees))
			(if u-direction (first knots) (second knots))
			(if u-direction
			    (array-get-column net i)
			    (array-get-row net i)))))

(defun bss-get-surface-curve (surface p &key (u-curve t))
  "Returns the u (v) curve at the specified v (u) parameter."
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let* ((1st (if u-curve 0 1))
	   (points (iter (for i below (nth 1st (array-dimensions net)))
			 (collect
			  (bsc-evaluate (get-curve-at-cp surface i :u-curve
							 (not u-curve))
					p)))))
      (make-bspline-curve (nth 1st degrees)
			  (copy-seq (nth 1st knots))
			  points))))

(defun fair-in-one-direction (surface resolution iteration distance
			      &key (u-direction t))
  (iter (with 1st = (if u-direction 0 1))
	(with 2nd = (if u-direction 1 0))
	(with lower = (bss-lower-parameter surface))
	(with upper = (bss-upper-parameter surface))
	(with max = (1- (nth 2nd resolution)))
	(for i to max)
	(for p = (nth 2nd (affine-combine lower (/ i max) upper)))
	(for curve = (bss-get-surface-curve surface p :u-curve u-direction))
	(collect (bsc-faired-polygon curve
				     (nth 1st resolution) iteration distance
				     :from (nth 1st lower)
				     :to (nth 1st upper)))))

(defun bss-faired-mesh (surface resolution iteration distance)
  "Results in RESOLUTION*RESOLUTION points, or RESOLUTION_1*RESOLUTION_2
points if RESOLUTION is a list."
  (let ((distance (* distance
		     (vlength (apply #'v- (bss-bounding-box surface))))))
    (unless (listp resolution)
      (setf resolution (list resolution resolution)))
    (let ((mesh-u (fair-in-one-direction surface resolution iteration
					 distance :u-direction t))
	  (mesh-v (fair-in-one-direction surface resolution iteration
					 distance :u-direction nil))
	  (result (make-array resolution)))
      (dotimes (i (first resolution))
	(dotimes (j (second resolution))
	  (setf (aref result i j) (affine-combine (elt (nth j mesh-u) i)
						  0.5
						  (elt (nth i mesh-v) j)))))
      result)))

;; (defun bss-fair (surface &key (measure #'fairness) (resolution 100)
;; 		 (target-iteration 100) (simplex-iteration 5)
;; 		 (fairing-iteration 5) (lock-endpoints t)
;; 		 (from (bss-lower-parameter surface))
;; 		 (to (bss-upper-parameter surface)))
;;   "Iterative fairing algorithm using the downhill simplex method.
;; MEASURE is a function giving a relative measure of the surface's fairness,
;; RESOLUTION is the sampling rate of the target curvature,
;; TARGET-ITERATION is the # of iterations in creating the target curvature,
;; SIMPLEX-ITERATION is the # of iterations the simplex method runs for a point,
;; FAIRING-ITERATION is the # of iterations every control point is faired,
;; and if LOCK-ENDPOINTS is T, the outer control points won't be moved.
;; TODO: it should be enough to check the fairness in the vicinity of the point.
;; TODO: the order of control point fairing should be based on the fairness."
;;   (let* ((parameters (arc-length-sampling curve from to resolution))
;; 	 (curvatures (target-curvature curve parameters target-iteration))
;; 	 (new-surface (copy-bspline-surface surface))
;; 	 (net (control-net new-surface))
;; 	 (size (array-dimensions net)))
;;     (dotimes (k fairing-iteration)
;;       (let ((low (if lock-endpoints 1 0))
;; 	    (high-u (if lock-endpoints (1- (first size)) (first size)))
;; 	    (high-v (if lock-endpoints (1- (second size)) (second size))))
;; 	(iter (for i from low below high-u)
;; 	      (iter (for j from low below high-v)
;; 		    (flet ((fairness-fn (point)
;; 			     (let ((old (elt points i j)))
;; 			       (setf (elt points i j) point)
;; 			       (prog1 (funcall measure new-surface
;; 					       parameters curvatures)
;; 				 (setf (elt points i j) old)))))
;; 		      (setf (elt points i j)
;; 			    (downhill-simplex:minimize
;; 			     #'fairness-fn (elt points i j)
;; 			     simplex-iteration)))))))
;;     new-surface))
