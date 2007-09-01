;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; B-spline surface

(defclass bspline-surface ()
  ((degrees
    :initarg :degrees
    :reader degrees
    :documentation "A list of two integers.")
   (knot-vectors
    :initarg :knot-vectors
    :accessor knot-vectors
    :documentation "A list of two sequence of floats.")
   (control-net
    :initarg :control-net
    :accessor control-net
    :documentation "A two-dimensional array of points.")))

(defun make-bspline-surface (degrees knot-vectors control-net)
  "DEGREES can be given as a single integer if both degrees are the same.
CONTROL-NET can be given as a nested list of points."
  (unless (listp degrees)
    (setf degrees (list degrees degrees)))
  (unless (arrayp control-net)
    (setf control-net (make-array (list (length control-net)
					(length (first control-net)))
				  :initial-contents control-net)))
  (let ((kl1 (length (first knot-vectors)))
	(kl2 (length (second knot-vectors)))
	(cl1 (array-dimension control-net 0))
	(cl2 (array-dimension control-net 1))
	(deg1 (first degrees))
	(deg2 (second degrees)))
    (assert (and (integerp deg1) (> deg1 0) (integerp deg2) (> deg2 0))
	    (degrees)
	    "Degrees should be positive integers, not ~a." degrees)
    (assert (and (= kl1 (+ cl1 deg1 1)) (= kl2 (+ cl2 deg2 1)))
	    (knot-vectors control-net)
	    "~d /= ~d + ~d + 1 or ~d /= ~d + ~d + 1~%~
             Knot vector length should be equal~%~
             to the # of control points + degree + 1" kl1 cl1 deg1 kl2 cl2 deg2)
    (make-instance 'bspline-surface
		   :degrees degrees
		   :knot-vectors knot-vectors
		   :control-net control-net)))

(defun copy-bspline-surface (surface)
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let ((new-net (make-array (array-dimensions net))))
      (dotimes (i (array-dimension net 0))
	(dotimes (j (array-dimension net 1))
	  (setf (aref new-net i j) (copy-list (aref net i j)))))
      (make-bspline-surface (copy-list degrees)
			    (list (copy-seq (first knots))
				  (copy-seq (second knots)))
			    new-net))))

(defun bsc-extrude (curve length)
  "Creates a BSPLINE-SURFACE from the extrusion of the two-dimensional
CURVE by LENGTH units along the z axis."
  (with-accessors ((degree degree)
		   (knots knot-vector)
		   (points control-points))
      curve
    (let ((net (map 'vector #'(lambda (x) (list (append x (list 0.0))
						(append x (list (/ length
								   2.0)))
						(append x (list length))))
		    points)))
      (make-bspline-surface (list degree 1)
			    (list (copy-seq knots)
				  (vector 0.0 0.0 0.5 1.0 1.0))
			    (make-array (list (length points) 3)
					:initial-contents net)))))

(defun bss-lower-parameter (surface)
  (with-accessors ((degrees degrees)
		   (knots knot-vectors))
      surface
    (list (elt (first knots) (first degrees))
	  (elt (second knots) (second degrees)))))

(defun bss-upper-parameter (surface)
  (with-accessors ((degrees degrees)
		   (knots knot-vectors))
      surface
    (let ((u-knots (first knots)) (v-knots (second knots)))
      (list (elt u-knots (- (length u-knots) (1+ (first degrees))))
	    (elt v-knots (- (length v-knots) (1+ (second degrees))))))))

(defun bss-bounding-box (surface)
  (with-accessors ((net control-net))
      surface
    (let* ((pts (iter (for j from 0 below (array-dimension net 1))
		      (append (iter (for i from 0 below (array-dimension net 0))
				    (collect (aref net i j))))))
	   (lst (apply #'mapcar #'list pts)))
      (list (mapcar #'(lambda (x) (apply #'min x)) lst)
	    (mapcar #'(lambda (x) (apply #'max x)) lst)))))


;; Evaluation

(defun find-bounding-intervals (surface uv)
  (let ((lower (bss-lower-parameter surface))
	(upper (bss-upper-parameter surface)))
    (mapcar #'lower-bound (knot-vectors surface) uv lower upper)))

(defun bss-evaluate (surface uv &key (derivative '(0 0)))
  "Evaluates SURFACE at PARAMETER and calculates the given DERIVATIVE, if given.
TODO: Slow (but easy) implementation that uses BSPLINE-CURVE's EVALUATE."
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let* ((bounds (find-bounding-intervals surface uv))
	   (u-bounds (first bounds)) (v-bounds (second bounds))
	   (u-degree (first degrees)) (v-degree (second degrees))
	   (u-points (make-array (array-dimension net 0)))
	   (v-points (make-array (array-dimension net 1))))
      (iter (for i from (- u-bounds u-degree) to u-bounds)
	    (for curve =
		 (iter (for j from (- v-bounds v-degree) to v-bounds)
		       (setf (elt v-points j) (aref net i j))
		       (finally (return (make-bspline-curve v-degree
							    (second knots)
							    v-points)))))
	    (setf (elt u-points i)
		  (bss-evaluate curve (second uv)
				:derivative (second derivative))))
      (bss-evaluate (make-bspline-curve u-degree (first knots) u-points)
		    (first uv) :derivative (first derivative)))))


;; High-level functions

(defun bss-surface-normal (surface uv)
  (let ((der-u (bss-evaluate surface uv :derivative '(1 0)))
	(der-v (bss-evaluate surface uv :derivative '(0 1))))
    (vnormalize (cross-product der-u der-v))))

(defun bss-principal-curvatures (surface uv)
  "List of the principal curvatures (KMAX, KMIN)."
  (let* ((der-00 (bss-evaluate surface uv :derivative '(1 0)))
	 (der-01 (bss-evaluate surface uv :derivative '(0 1)))
	 (der-10 (bss-evaluate surface uv :derivative '(2 0)))
	 (der-11 (bss-evaluate surface uv :derivative '(1 1)))
	 (der-12 (bss-evaluate surface uv :derivative '(0 2)))
	 (normal (vnormalize (cross-product der-00 der-01)))
	 (E (scalar-product der-00 der-00))
	 (F (scalar-product der-00 der-01))
	 (G (scalar-product der-01 der-01))
	 (L (scalar-product normal der-10))
	 (M (scalar-product normal der-11))
	 (N (scalar-product normal der-12))
	 (m (+ (* N E) (* -2 F M) (* L G)))
	 (g (* 2 (- (* L N) (* M M))))
	 (a (* 2 (- (* E G) (* F F))))
	 (d (sqrt (- (* m m) (* a g)))))
    (list (/ (+ m d) a) (/ (- m d) a))))

(defun bss-gaussian-curvature (surface uv)
  (let* ((der-00 (bss-evaluate surface uv :derivative '(1 0)))
	 (der-01 (bss-evaluate surface uv :derivative '(0 1)))
	 (der-10 (bss-evaluate surface uv :derivative '(2 0)))
	 (der-11 (bss-evaluate surface uv :derivative '(1 1)))
	 (der-12 (bss-evaluate surface uv :derivative '(0 2)))
	 (normal (vnormalize (cross-product der-00 der-01)))
	 (E (scalar-product der-00 der-00))
	 (F (scalar-product der-00 der-01))
	 (G (scalar-product der-01 der-01))
	 (L (scalar-product normal der-10))
	 (M (scalar-product normal der-11))
	 (N (scalar-product normal der-12)))
    (safe-/ (- (* L N) (* M M)) (- (* E G) (* F F)))))

(defun bss-mean-curvature (surface uv)
  (let* ((der-00 (bss-evaluate surface uv :derivative '(1 0)))
	 (der-01 (bss-evaluate surface uv :derivative '(0 1)))
	 (der-10 (bss-evaluate surface uv :derivative '(2 0)))
	 (der-11 (bss-evaluate surface uv :derivative '(1 1)))
	 (der-12 (bss-evaluate surface uv :derivative '(0 2)))
	 (normal (vnormalize (cross-product der-00 der-01)))
	 (E (scalar-product der-00 der-00))
	 (F (scalar-product der-00 der-01))
	 (G (scalar-product der-01 der-01))
	 (L (scalar-product normal der-10))
	 (M (scalar-product normal der-11))
	 (N (scalar-product normal der-12)))
    (safe-/ (+ (* N E) (* -2 M F) (* L G)) (- (* E G) (* F F)))))

;; (defun insert-knot-u (bspline-surface u)
;; "Inserts the not U into the knot vector of SURFACE."
;;   (with-accessors ((degrees degrees)
;; 		   (knots knot-vectors)
;; 		   (net control-net))
;;       bspline-surface
;;     (let ((new-knots (sort (concatenate 'vector (first knots) (list u)) #'<))
;; 	  (new-net (make-array (list (1+ (array-dimension net 0))
;; 				     (array-dimension net 1))))
;; 	  (u-degree (first degrees))
;; 	  (u-knots (first knots)))
;;       (let ((k (1- (position u new-knots :from-end t)))
;; 	    (alpha (make-array u-degree)))
;; 	(loop for i from 0 below u-degree
;; 	      for j = (+ k (- u-degree) 1 i)
;; 	      do (setf (elt alpha i)
;; 		       (/ (- u (elt u-knots j))
;; 			  (- (elt u-knots (+ j u-degree)) (elt u-knots j)))))
;; 	(loop for j from 0 below (array-dimension net 1) do
;; 	      (loop for i from 0 to (- k u-degree) do
;; 		    (setf (aref new-net i j) (copy-list (aref net i j))))
;; 	      (loop for i from (+ k (- u-degree) 1) to k do
;; 		    (setf (aref new-net i j)
;; 			  (affine-combine (aref net (1- i) j)
;; 					  (elt alpha (- i k (- u-degree) 1))
;; 					  (aref net i j))))
;; 	      (loop for i from (1+ k) below (array-dimension new-net 0) do
;; 		    (setf (aref new-net i j)
;; 			  (copy-list (aref net (1- i) j))))))
;;       (make-bspline-surface (copy-list degrees)
;; 			    (list new-knots (copy-seq (second knots)))
;; 			    new-net))))

;; (defun insert-knot-v (bspline-surface v)
;; "Inserts the not V into the knot vector of SURFACE."
;;   (with-accessors ((degrees degrees)
;; 		   (knots knot-vectors)
;; 		   (net control-net))
;;       bspline-surface
;;     (let ((new-knots (sort (concatenate 'vector (second knots) (list v)) #'<))
;; 	  (new-net (make-array (list (array-dimension net 0)
;; 				     (1+ (array-dimension net 1)))))
;; 	  (v-degree (second degrees))
;; 	  (v-knots (second knots)))
;;       (let ((k (1- (position v new-knots :from-end t)))
;; 	    (alpha (make-array v-degree)))
;; 	(loop for i from 0 below v-degree
;; 	      for j = (+ k (- v-degree) 1 i)
;; 	      do (setf (elt alpha i)
;; 		       (/ (- v (elt v-knots j))
;; 			  (- (elt v-knots (+ j v-degree)) (elt v-knots j)))))
;; 	(loop for i from 0 below (array-dimension net 0) do
;; 	      (loop for j from 0 to (- k v-degree) do
;; 		    (setf (aref new-net i j) (copy-list (aref net i j))))
;; 	      (loop for j from (+ k (- v-degree) 1) to k do
;; 		    (setf (aref new-net i j)
;; 			  (affine-combine (aref net i (1- j))
;; 					  (elt alpha (- j k (- v-degree) 1))
;; 					  (aref net i j))))
;; 	      (loop for j from (1+ k) below (array-dimension new-net 1) do
;; 		    (setf (aref new-net i j)
;; 			  (copy-list (aref net i (1- j)))))))
;;       (make-bspline-surface (copy-list degrees)
;; 			    (list (copy-seq (first knots)) new-knots)
;; 			    new-net))))

;; (defmethod insert-knot ((surface bspline-surface) uv &key (u-direction t))
;;   (if u-direction
;;       (insert-knot-u surface uv)
;;       (insert-knot-v surface uv)))

;; (defun split-surface-u (bspline-surface u)
;;   (let ((surface bspline-surface)
;; 	(degrees (degrees bspline-surface)))
;;     (loop repeat (- (first degrees)
;; 		    (count u (first (knot-vectors bspline-surface)))) do
;; 	  (setf surface (insert-knot surface u :u-direction t)))
;;     (with-accessors ((knots knot-vectors)
;; 		     (net control-net))
;; 	surface
;;       (let* ((main-knots (first knots))
;; 	     (k1 (1+ (position u main-knots :from-end t)))
;; 	     (k2 (position u main-knots))
;; 	     (new-knots-1 (list (concatenate 'vector
;; 					     (subseq main-knots 0 k1) (list u))
;; 				(second knots)))
;; 	     (new-knots-2 (list (concatenate 'vector
;; 					     (list u) (subseq main-knots k2))
;; 				(second knots))))
;; 	(list (make-bspline-surface (copy-list degrees)
;; 				    new-knots-1
;; 				    (subarray-2d net 0 0
;; 						 (- k1 (first degrees)) nil))
;; 	      (make-bspline-surface (copy-list degrees)
;; 				    new-knots-2
;; 				    (subarray-2d net (1- k2) 0 nil nil)))))))

;; (defun split-surface-v (bspline-surface v)
;;   (let ((surface bspline-surface)
;; 	(degrees (degrees bspline-surface)))
;;     (loop repeat (- (second degrees)
;; 		    (count v (second (knot-vectors bspline-surface)))) do
;; 	  (setf surface (insert-knot surface v :u-direction nil)))
;;     (with-accessors ((knots knot-vectors)
;; 		     (net control-net))
;; 	surface
;;       (let* ((main-knots (second knots))
;; 	     (k1 (1+ (position v main-knots :from-end t)))
;; 	     (k2 (position v main-knots))
;; 	     (new-knots-1 (list (first knots)
;; 				(concatenate 'vector
;; 					     (subseq main-knots 0 k1)
;; 					     (list v))))
;; 	     (new-knots-2 (list (first knots)
;; 				(concatenate 'vector
;; 					     (list v)
;; 					     (subseq main-knots k2)))))
;; 	(list (make-bspline-surface (copy-list degrees)
;; 				    new-knots-1
;; 				    (subarray-2d net 0 0
;; 						 nil (- k1 (first degrees))))
;; 	      (make-bspline-surface (copy-list degrees)
;; 				    new-knots-2
;; 				    (subarray-2d net 0 (1- k2) nil nil)))))))

;; (defun split-surface (bspline-surface uv &key (u-direction t))
;;   (if u-direction
;;       (split-surface-u bspline-surface uv)
;;       (split-surface-v bspline-surface uv)))

;; (defun subsurface-1 (bspline-surface min max u-dir)
;;   (cond ((and min max)
;; 	 (second (split-surface (first (split-surface bspline-surface max
;; 						      :u-direction u-dir))
;; 				min :u-direction u-dir)))
;; 	(min
;; 	 (second (split-surface bspline-surface min :u-direction u-dir)))
;; 	(max
;; 	 (first (split-surface bspline-surface max :u-direction u-dir)))
;; 	(t
;; 	 bspline-surface)))

;; (defun subsurface (bspline-surface min-u min-v max-u max-v)
;;   (subsurface-1 (subsurface-1 bspline-surface min-u max-u t) min-v max-v nil))
