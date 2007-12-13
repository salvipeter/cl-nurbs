;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Fairing

(defun bsc-knot-removal-reinsertion (curve u)
  "Removes and reinserts a knot by moving the specified control point."
  (with-accessors ((degree degree)
		   (knots knot-vector)
		   (points control-points))
      curve
    (unless (= degree 3)
      (error "this function only works for 3rd degree bspline curves"))
    (let* ((i (+ u 2))
	   (left (affine-combine (elt points (- u 1))
				 (safe-/ (- (elt knots i) (elt knots (+ i 1)))
					 (- (elt knots i) (elt knots (- i 3))))
				 (elt points (- u 2))))
	   (right (affine-combine (elt points (+ u 1))
				  (safe-/ (- (elt knots (- i 1)) (elt knots i))
					  (- (elt knots (+ i 3))
					     (elt knots i)))
				  (elt points (+ u 2)))))
      (affine-combine left
		      (safe-/ (- (elt knots i) (elt knots (- i 2)))
			      (- (elt knots (+ i 2)) (elt knots (- i 2))))
		      right))))

(defun arc-length-sampling (curve from to resolution &optional (iteration 20))
  "Returns a vector of RESOLUTION parameter values that define
equivalent-length segments of CURVE in the [FROM, TO] interval."
  (let ((step (/ (- to from) (1- resolution)))
	(step-length (/ (bsc-estimate-arc-length curve from to)
			(1- resolution)))
	(result (make-array resolution)))
    (setf (elt result 0) from)
    (iter (for i from 1 below (1- resolution))
	  (for low = (elt result (1- i)))
	  (for upp = (min to (+ low (* step 2.0))))
	  (setf (elt result i)
		(iter (repeat iteration)
		      (for mid = (/ (+ low upp) 2.0))
		      (for segment = (bsc-estimate-arc-length curve from mid))
		      (for correct = (* step-length i))
		      (cond ((< segment correct) (setf low mid))
			    ((> segment correct) (setf upp mid))
			    (t nil))
		      (finally (return mid)))))
    (setf (elt result (1- resolution)) to)
    result))

(defun segment-lengths (curve parameters &optional resolution)
  "Returns a vector of the segment lengths for testing purposes."
  (map 'vector
       (if resolution
	   #'(lambda (x y) (bsc-iterative-arc-length curve x y resolution))
	   #'(lambda (x y) (bsc-estimate-arc-length curve x y)))
       (subseq parameters 0 (1- (length parameters))) (subseq parameters 1)))

(defun target-curvature (curve parameters iteration &key start-value end-value)
  "Averages the curvature of CURVE at PARAMETERS ITERATION times.
The first and last values can be fixed through START-VALUE and END-VALUE."
  (let ((curvatures (map 'list #'(lambda (x) (bsc-curvature curve x))
			 parameters)))
    (when start-value (setf (first curvatures) start-value))
    (when end-value (setf (car (last curvatures)) end-value))
    (iter (repeat iteration)
	  (setf curvatures (append (list (first curvatures))
				   (mapcar #'(lambda (x y) (/ (+ x y) 2.0))
					   (butlast curvatures 2)
					   (rest (rest curvatures)))
				   (last curvatures)))
	  (finally (return (coerce curvatures 'vector))))))

(defun integrate (curve parameters curvatures distance &key from-right)
  "Euler-integration at PARAMETERS of CURVE,using CURVATURES
as the second derivatives and DISTANCE as the deviation constraint.
TODO: Slow because of unncecessary evaluations."
  (let* ((resolution (length parameters))
	 (points (bsc-evaluate-on-parameters curve parameters))
	 (outward (bsc-out-direction-on-parameters curve parameters))
	 (step (/ (bsc-estimate-arc-length curve
					   (elt parameters 0)
					   (elt parameters (1- resolution)))
		  (1- resolution)))
	 (1/step (/ 1.0 step))
	 (result (make-array resolution)))
    (flet ((flip (x) (if from-right (- resolution 1 x) x))
	   (sqr (x) (* x x)))
      (setf (elt result 0) (copy-list (elt points (flip 0))))
      (iter (with last = (v* (v- (elt points (flip 1)) (elt points (flip 0)))
			     1/step))
	    (for i from 0 below (1- resolution))
	    (for deviation = (v- (elt points (flip i)) (elt result i)))
	    (for offset = (v+ (v* last step)
			      (v* deviation
				  (min (sqr (/ (vlength deviation) distance))
				       1.0))))
	    (for last-offset = (v* (elt outward (flip i))
				   (elt curvatures (flip i)) step))
	    (setf (elt result (1+ i)) (v+ (elt result i) offset)
		  last (vnormalize (v+ last last-offset)))
	    (finally (return result))))))

;; (defun blend-function (x)
;;   "Blend function in the interval [0, 1]: 3x^2-2x^3."
;;   (- (* 3.0 (* x x)) (* 2.0 (expt x 3))))

(defun blend-function (x)
  "Blend function in the interval [0, 1]: 6x^5-15x^4+10x^3."
  (+ (* 6.0 (expt x 5)) (- (* 15.0 (expt x 4))) (* 10.0 (expt x 3))))

(defun blend-points (left right)
  "Blends the two point vectors."
  (let* ((resolution (length left))
	 (result (make-array resolution)))
    (iter (for i from 0 below resolution)
	  (for alpha = (blend-function (/ i (1- resolution))))
	  (setf (elt result i)
		(affine-combine (elt left i) alpha (elt right i)))
	  (finally (return result)))))

(defun bsc-faired-polygon (curve resolution iteration distance
			   &key from to start-curvature end-curvature)
  "Fairs the [FROM, TO] interval of CURVE. The number of
resulting points is based on RESOLUTION; the target curvature is averaged
ITERATION times and the deviation of the points is constrained to be less
than DISTANCE.
Results in RESOLUTION points."
  (let* ((from (or from (bsc-lower-parameter curve)))
	 (to (or to (bsc-upper-parameter curve)))
	 (parameters (arc-length-sampling curve from to resolution))
	 (curvatures (target-curvature curve parameters iteration
				       :start-value start-curvature
				       :end-value end-curvature))
	 (left (integrate curve parameters curvatures distance :from-right nil))
	 (right (integrate curve parameters curvatures distance :from-right t)))
    (blend-points left (nreverse right))))

(defun fairness (curve parameters curvatures)
  "A fairness measure: the squared sum of the differences
of the real curvature from the target curvature."
  (let ((squared-differences
	 (map 'list #'(lambda (x y)
			(let ((difference (- (bsc-curvature curve x) y)))
			  (* difference difference)))
	      parameters curvatures)))
    (apply #'+ squared-differences)))

(defun bsc-fair (curve &key (measure #'fairness) (resolution 100)
		 (target-iteration 100) start-curvature end-curvature
		 (simplex-iteration 5) (fairing-iteration 5) (lock-endpoints t)
		 (from (bsc-lower-parameter curve))
		 (to (bsc-upper-parameter curve)))
  "Iterative fairing algorithm using the downhill simplex method.
MEASURE is a function giving a relative measure of the curve's fairness,
RESOLUTION is the sampling rate of the target curvature,
TARGET-ITERATION is the # of iterations in creating the target curvature,
SIMPLEX-ITERATION is the # of iterations the simplex method runs for a point,
FAIRING-ITERATION is the # of iterations every control point is faired,
and if LOCK-ENDPOINTS is T, the endpoints of the curve won't be moved.
TODO: it should be enough to check the fairness in the vicinity of the point.
TODO: the order of control point fairing should be based on the fairness."
  (let* ((parameters (arc-length-sampling curve from to resolution))
	 (curvatures (target-curvature curve parameters target-iteration
				       :start-value start-curvature
				       :end-value end-curvature))
	 (new-curve (copy-bspline-curve curve))
	 (points (control-points new-curve)))
    (dotimes (i fairing-iteration)
      (let ((low (if lock-endpoints 1 0))
	    (high (if lock-endpoints (1- (length points)) (length points))))
	(iter (for j from low below high)
	      (flet ((fairness-fn (point)
		       (let ((old (elt points j)))
			 (setf (elt points j) point)
			 (prog1
			     (funcall measure new-curve parameters curvatures)
			   (setf (elt points j) old)))))
		(setf (elt points j)
		      (downhill-simplex:minimize
		       #'fairness-fn (elt points j) simplex-iteration))))))
    new-curve))
