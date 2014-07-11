;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; B-spline curve

(defclass bspline-curve ()
  ((degree
    :initarg :degree
    :reader degree)
   (knot-vector
    :initarg :knot-vector
    :accessor knot-vector
    :documentation "A sequence of floats.")
   (control-points
    :initarg :control-points
    :accessor control-points
    :documentation "A sequence of points.")))

(defun make-bspline-curve (degree knot-vector control-points)
  (assert (and (integerp degree) (> degree 0)) (degree)
	  "Degree should be a positive integer, not ~d." degree)
  (let ((kl (length knot-vector))
	(cl (length control-points)))
    (assert (= kl (+ cl degree 1)) (knot-vector control-points)
	    "~d /= ~d + ~d + 1~%~
             Knot vector length should be equal~%~
             to the # of control points + degree + 1" kl cl degree)
    (let ((control-array (if (arrayp control-points)
			     control-points
			     (make-array cl :initial-contents control-points))))
      (make-instance 'bspline-curve
		     :degree degree
		     :knot-vector knot-vector
		     :control-points control-array))))

(defun copy-bspline-curve (bspline-curve)
  (with-accessors ((degree degree)
		   (knots knot-vector)
		   (points control-points))
      bspline-curve
    (let ((new-points (make-array (length points))))
      (dotimes (i (length points))
	(setf (elt new-points i) (copy-list (elt points i))))
      (make-bspline-curve degree (copy-seq knots) new-points))))

(defun bsc-dimension (curve)
  (length (elt (control-points curve) 0)))

(defun bsc-lower-parameter (curve)
  (elt (knot-vector curve) (degree curve)))

(defun bsc-upper-parameter (curve)
  (let ((knots (knot-vector curve)))
    (elt knots (- (length knots) (1+ (degree curve))))))

(defun bsc-bounding-box (curve)
  (let ((lst (apply #'mapcar #'list (coerce (control-points curve) 'list))))
    (list (mapcar #'(lambda (x) (apply #'min x)) lst)
	  (mapcar #'(lambda (x) (apply #'max x)) lst))))

(defun bsc-bounding-box-axis (curve)
  (vlength (apply #'v- (bsc-bounding-box curve))))


;; Evaluation

(let ((left 1))
  (defun lower-bound (seq x lower upper)
    "Finds the largest I that SEQ[I] <= X.
If X = UPPER it returns the largest I that SEQ[I] < X.
Uses previous result as a first guess."
    (unless (<= lower x upper)
      (error "~f is out of bounds: [~f; ~f]" x lower upper))
    (cond ((= x upper) (setf left (position x seq :from-end t :test #'>)))
	  ((and (< left (1- (length seq)))
		(>= x (elt seq left))
		(< x (elt seq (1+ left))))
	   left)
	  (t (setf left (position x seq :from-end t :test #'>=))))))

(defun find-bounding-interval (curve u)
  (lower-bound (knot-vector curve) u
	       (bsc-lower-parameter curve)
	       (bsc-upper-parameter curve)))

(defun eval-blossom (points knots u degree k derivative)
  "Calculates the Kth level blossom for the parameters U[0..DEGREE].
POINTS is a DEGREE+1, KNOTS is a 2*DEGREE long sequence.
Calculate the derivative triangles up to level DERIVATIVE."
  (if (< k 0)
      points
      (let ((blossom (eval-blossom points knots u degree (1- k) derivative)))
	(if (>= k derivative)
	    (iter (for j from 0 below (- degree k))
		  (collect (affine-combine (elt blossom j)
					   (safe-/ (- (elt u k)
						      (elt knots (+ k j)))
						   (- (elt knots (+ degree j))
						      (elt knots (+ k j))))
					   (elt blossom (1+ j)))))
	    (iter (for j from 0 below (- degree k))
		  (collect (v* (v- (elt blossom (1+ j)) (elt blossom j))
			       (safe-/ (- degree k)
				       (- (elt knots (+ degree j))
					  (elt knots (+ k j)))))))))))

(defun bsc-evaluate (curve u &key (derivative 0))
  "Evaluates CURVE at PARAMETER and calculates the given DERIVATIVE, if given.
Uses deBoor blossoming.
TODO: Should return the list of derivatives from 0 to DERIVATIVE."
  (let ((degree (degree curve))
	(r (find-bounding-interval curve u)))
    (elt (eval-blossom (subseq (control-points curve) (- r degree) (1+ r))
		       (subseq (knot-vector curve)
			       (- r (1- degree)) (+ r degree 1))
		       (make-array degree :initial-element u)
		       degree (1- degree) derivative) 0)))

(defun bsc-evaluate-on-parameters (curve parameters)
  "Convenience function for evaluation on a sequence of parameters."
  (map 'vector #'(lambda (x) (bsc-evaluate curve x)) parameters))


;; High-level functions

(defun bsc-2d-normal (curve u)
  "Normal vector of a 2-dimensional B-spline curve."
  (let ((d (bsc-evaluate curve u :derivative 1)))
    (vnormalize (list (second d) (- (first d))))))

(defun bsc-out-direction (curve u)
  "Outward direction on the osculating plane at U for a 3-dimensional spline."
  (let ((d1 (bsc-evaluate curve u :derivative 1))
	(d2 (bsc-evaluate curve u :derivative 2)))
    (vnormalize (cross-product (cross-product d1 d2) d1))))

(defun bsc-out-direction-on-parameters (curve parameters)
  "Convenience function for calculating outward direction on a
sequence of parameters."
  (if (= (bsc-dimension curve) 2)
      (map 'vector #'(lambda (x) (bsc-2d-normal curve x)) parameters)
      (map 'vector #'(lambda (x) (bsc-out-direction curve x)) parameters)))

(defun bsc-curvature (curve u)
  (let ((d1 (bsc-evaluate curve u :derivative 1))
	(d2 (bsc-evaluate curve u :derivative 2)))
    (if (= (bsc-dimension curve) 2)
	(safe-/ (scalar-product d1 (list (- (second d2)) (first d2)))
		(expt (vlength d1) 3))
	(safe-/ (vlength (cross-product d1 d2)) (expt (vlength d1) 3)))))

(defun bsc-torsion (curve u)
  (let* ((d1 (bsc-evaluate curve u :derivative 1))
         (d2 (bsc-evaluate curve u :derivative 2))
         (d3 (bsc-evaluate curve u :derivative 3))
         (cross (cross-product d1 d2)))
    (safe-/ (scalar-product cross d3) (vlength2 cross))))

(defun bsc-insert-knot (curve u &optional (repetition 1))
  "Inserts the knot U into the knot vector of CURVE REPETITION times.

Translation of the algorithm in the NURBS Book, pp. 151."
  (let* ((p (degree curve))
	 (knots (knot-vector curve))
	 (points (control-points curve))
	 (s (count u knots))
	 (k (position u knots :test #'>= :from-end t)))
    (assert (<= (+ repetition s) p))
    (let ((new-knots (concatenate 'vector
				  (subseq knots 0 (1+ k))
				  (iter (repeat repetition) (collect u))
				  (subseq knots (1+ k))))
	  (temp-array (make-array (1+ p)))
	  (new-points (make-array (+ (length points) repetition))))
      (iter (for i from 0 to (- k p))
	    (setf (elt new-points i) (elt points i)))
      (iter (for i from (- k s) below (length points))
	    (setf (elt new-points (+ i repetition)) (elt points i)))
      (iter (for i from 0 to (- p s))
	    (setf (elt temp-array i) (elt points (+ (- k p) i))))
      (iter (for j from 1 to repetition)
	    (for L = (+ (- k p) j))
	    (iter (for i from 0 to (- p j s))
		  (let ((alpha (/ (- u (elt knots (+ L i)))
				  (- (elt knots (+ i k 1))
				     (elt knots (+ L i))))))
		    (setf (elt temp-array i)
			  (affine-combine (elt temp-array i)
					  alpha
					  (elt temp-array (1+ i))))))
	    (setf (elt new-points L)
		  (elt temp-array 0)
		  (elt new-points (- (+ k repetition) j s))
		  (elt temp-array (- p j s))))
      (let ((L (+ (- k p) repetition)))
	(iter (for i from (1+ L) below (- k s))
	      (setf (elt new-points i) (elt temp-array (- i L)))))
      (make-bspline-curve p new-knots new-points))))

(defun bsc-split-curve (curve u)
  (let ((degree (degree curve)))
    (iter (repeat (- degree (count u (knot-vector curve))))
	  (setf curve (bsc-insert-knot curve u)))
    (with-accessors ((knots knot-vector)
		     (points control-points))
	curve
      (let ((k1 (1+ (position u knots :from-end t)))
	    (k2 (position u knots)))
	(list (make-bspline-curve degree
				  (concatenate 'vector
					       (subseq knots 0 k1) (list u))
				  (subseq points 0 (- k1 degree)))
	      (make-bspline-curve degree
				  (concatenate 'vector
					       (list u) (subseq knots k2))
				  (subseq points (1- k2))))))))

(defun bsc-subcurve (curve min max)
  (cond ((and min max)
	 (second (bsc-split-curve (first (bsc-split-curve curve max)) min)))
	(min
	 (second (bsc-split-curve curve min)))
	(max
	 (first (bsc-split-curve curve max)))
	(t curve)))

(define-constant +gaussian-quadrature+
  '((-0.861136312 0.347854845) (-0.339981044 0.652145155)
    (0.339981044 0.652145155) (0.861136312 0.347854845))
  "Four-point Gaussian quadrature.")

;; (define-constant +gaussian-quadrature+
;;   '((-0.906180 0.236927) (-0.538469 0.478629) (0 0.568889)
;;     (0.538469 0.478629) (0.906180 0.236927))
;;   "Five-point Gaussian quadrature.")

(defun bsc-estimate-arc-length (curve &optional
				(from (bsc-lower-parameter curve))
				(to (bsc-upper-parameter curve)))
  "Estimates the arc length of CURVE in the
[FROM, TO] parameter interval, using Gaussian quadratures."
  (if (>= from to)
      0.0
      (let ((next (min to (elt (knot-vector curve)
			       (1+ (find-bounding-interval curve from))))))
	(+ (iter (for gauss in +gaussian-quadrature+)
		 (for u = (/ (+ (* (- next from) (first gauss)) from next) 2.0))
		 (for dn = (vlength (bsc-evaluate curve u :derivative 1)))
		 (sum (* dn (second gauss) (- next from) 0.5)))
	   (bsc-estimate-arc-length curve next to)))))

(defun bsc-iterative-arc-length (curve begin end &optional (resolution 100))
  "Slow (but arbitrarily accurate) arc length calculation."
  (iter (with step = (/ (- end begin) resolution))
	(for u from begin below end by step)
	(sum (vlength (v- (bsc-evaluate curve (min end (+ u step)))
			  (bsc-evaluate curve u))))))
