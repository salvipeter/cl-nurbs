;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

(defun reparametrize-curve (bspline-curve min max)
  (let* ((old-min (lower-parameter bspline-curve))
	 (old-max (upper-parameter bspline-curve))
	 (knots (knot-vector bspline-curve))
	 (n (length knots))
	 (points (control-points bspline-curve))
	 (m (length points))
	 (new-knots (make-array n))
	 (new-points (make-array m)))
    (dotimes (i n)
      (setf (elt new-knots i)
	    (+ min (* (- max min) (/ (- (elt knots i) old-min)
				     (- old-max old-min))))))
    (dotimes (i m)
      (setf (elt new-points i) (copy-list (elt points i))))
    (make-bspline-curve (degree bspline-curve) new-knots new-points)))

(defun flip-uv (surface)
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let* ((n (array-dimension net 0))
	   (m (array-dimension net 1))
	   (new-net (make-array (list m n))))
      (dotimes (i n)
	(dotimes (j m)
	  (setf (aref new-net j i) (aref net i j))))
      (make-bspline-surface (reverse degrees) (reverse knots) new-net))))

(defun reverse-parametrization (surface &key u v)
  (with-accessors ((degrees degrees)
		   (knots knot-vectors)
		   (net control-net))
      surface
    (let* ((n (array-dimension net 0))
	   (m (array-dimension net 1))
	   (k (length (first knots)))
	   (l (length (second knots)))
	   (new-knots (list (if u (make-array k) (copy-seq (first knots)))
			    (if v (make-array l) (copy-seq (second knots)))))
	   (new-net (make-array (list n m))))
      (when u
	(let ((low (first (first knots)))
	      (high (elt (first knots) (1- k))))
	  (dotimes (i k)
	    (setf (elt (first new-knots) i)
		  (+ low (- high (elt (first knots) (- k i 1))))))))
      (when v
	(let ((low (first (second knots)))
	      (high (elt (second knots) (1- l))))
	  (dotimes (i l)
	    (setf (elt (second new-knots) i)
		  (+ low (- high (elt (second knots) (- l i 1))))))))
      (dotimes (i n)
	(dotimes (j m)
	  (setf (aref new-net i j)
		(aref net (if u (- n i 1) i) (if v (- m j 1) j)))))
      (make-bspline-surface (copy-list degrees) new-knots new-net))))

(defun timed-format (stream string &rest rest)
  (let ((now (multiple-value-list (get-decoded-time))))
  (apply #'format
	 stream (concatenate 'string "[~2,'0d:~2,'0d:~2,'0d] " string)
	 (third now) (second now) (first now) rest)))

(defun u-length (surface)
  (vector-length
   (vector-subtract
    (evaluate surface (lower-parameter surface))
    (evaluate surface (list (first (upper-parameter surface))
			    (second (lower-parameter surface)))))))

(defun v-length (surface)
  (vector-length
   (vector-subtract
    (evaluate surface (lower-parameter surface))
    (evaluate surface (list (first (lower-parameter surface))
			    (second (upper-parameter surface)))))))

;; (defun find-parameter (x min max fn iterations)
;;   (let* ((med (/ (+ min max) 2.0))
;; 	 (fn-med (apply fn med)))
;;     (cond ((or (= iterations 0) (= fn-med x)) med)
;; 	  ((< fn-med x) (find-parameter x med max fn (1- iterations)))
;; 	  (t (find-parameter x min med fn (1- iterations))))))

(defun write-deviation (surface reference filename resolution)
  (unless (listp resolution)
    (setf resolution (list resolution resolution)))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "# vtk DataFile Version 1.0~
               ~%B-Spline Surface~
               ~%ASCII~
               ~%DATASET POLYDATA~%~
               ~%POINTS ~d float~%" (* (first resolution) (second resolution)))
    (let* ((lower (lower-parameter surface))
	   (upper (upper-parameter surface))
	   (lower-ref (lower-parameter reference))
	   (upper-ref (upper-parameter reference))
	   (u-res (first resolution))
	   (v-res (second resolution))
	   distances)
      (dotimes (j v-res)
	(dotimes (i u-res)
	  (let* ((uv (list (interpolate (first lower)
					(/ i (1- u-res))
					(first upper))
			   (interpolate (second lower)
					(/ j (1- v-res))
					(second upper))))
		 (p (evaluate surface uv))
		 (r (evaluate reference
			      (list (interpolate (first lower-ref)
						 (/ i (1- u-res))
						 (first upper-ref))
				    (interpolate (second lower-ref)
						 (/ j (1- v-res))
						 (second upper-ref))))))
	    (format s "~{~f ~}~%" p)
	    (let ((deviation (vector-subtract p r))
		  (normal (surface-normal surface uv)))
	      (push (vector-scalar-product deviation normal) distances)))))
      (setf distances (nreverse distances))
      (format s "~%POLYGONS ~d ~d~%"
	      (* (1- u-res) (1- v-res)) (* (1- u-res) (1- v-res) 5))
      (dotimes (j (1- v-res))
	(dotimes (i (1- u-res))
	  (format s "4 ~d ~d ~d ~d~%"
		  (+ (* j u-res) i) (+ (* j u-res) i 1)
		  (+ (* j u-res) i u-res 1) (+ (* j u-res) i u-res))))
      (format s "~%POINT_DATA ~d~
                 ~%COLOR_SCALARS deviation 3~%"
	      (* (first resolution) (second resolution)))
      (let ((max (loop for i in distances maximize (abs i))))
	(loop for i in distances
	      for d = (abs (/ i max))
	      do (if (< i 0)
		     (format s "0.0 ~f ~f~%" (- 1.0 d) d)
		     (format s "~f ~f 0.0~%" d (- 1.0 d))))))))

(defun write-uv-surface (surface filename resolution)
  (unless (listp resolution)
    (setf resolution (list resolution resolution)))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "# vtk DataFile Version 1.0~
               ~%B-Spline Surface~
               ~%ASCII~
               ~%DATASET POLYDATA~%~
               ~%POINTS ~d float~%" (* (first resolution) (second resolution)))
    (let* ((lower (lower-parameter surface))
	   (upper (upper-parameter surface))
	   (u-res (first resolution))
	   (v-res (second resolution))
	   colors)
      (dotimes (j v-res)
	(dotimes (i u-res)
	  (let* ((uv (list (interpolate (first lower)
					(/ i (1- u-res))
					(first upper))
			   (interpolate (second lower)
					(/ j (1- v-res))
					(second upper))))
		 (p (evaluate surface uv)))
	    (format s "~{~f ~}~%" p)
	    (push (list (/ i (1- u-res)) (/ j (1- v-res))) colors))))
      (setf colors (nreverse colors))
      (format s "~%POLYGONS ~d ~d~%"
	      (* (1- u-res) (1- v-res)) (* (1- u-res) (1- v-res) 5))
      (dotimes (j (1- v-res))
	(dotimes (i (1- u-res))
	  (format s "4 ~d ~d ~d ~d~%"
		  (+ (* j u-res) i) (+ (* j u-res) i 1)
		  (+ (* j u-res) i u-res 1) (+ (* j u-res) i u-res))))
      (format s "~%POINT_DATA ~d~
                 ~%COLOR_SCALARS deviation 3~%"
	      (* (first resolution) (second resolution)))
      (loop for i in colors do (format s "~{~f ~}0.0~%" i)))))

(defun zap-to-curve (surface curve)
  "Zaps the u=0 isocurve of SURFACE to CURVE."
  (let ((curve (reparametrize-curve curve
				    (second (lower-parameter surface))
				    (second (upper-parameter surface))))
	(surface (copy-bspline-surface surface)))
    (loop with i = 0
	  for curve-short = (>= i (length (knot-vector curve)))
	  for surface-short = (>= i (length (second (knot-vectors surface))))
	  while (and (not curve-short) (not surface-short)) do
	  (if (or curve-short
		  (and (not surface-short)
		       (> (elt (knot-vector curve) i)
			  (elt (second (knot-vectors surface)) i))))
	      (setf curve
		    (insert-knot curve
				 (elt (second (knot-vectors surface)) i)))
	      (when (or surface-short
			(< (elt (knot-vector curve) i)
			   (elt (second (knot-vectors surface)) i)))
		(setf surface
		      (insert-knot surface
				   (elt (knot-vector curve) i)
				   :u-direction nil))))
	  (setf i (1+ i)))
    (dotimes (i (length (control-points curve)))
      (setf (aref (control-net surface) 0 i) (elt (control-points curve) i)))
    surface))

(defun points-bbox (points)
  (let ((flat (do ((i 0 (1+ i))
		   acc)
		  ((= i (array-dimension points 0)) acc)
		(do ((j 0 (1+ j)))
		    ((= j (array-dimension points 1)))
		  (push (aref points i j) acc)))))
    (let ((lst (apply #'mapcar #'list flat)))
      (list (mapcar #'(lambda (x) (apply #'min x)) lst)
	    (mapcar #'(lambda (x) (apply #'max x)) lst)))))

(defun scale-surface (surface &optional (scaling 1000))
  (let* ((net (control-net surface))
	 (size (array-dimensions net))
	 (new-net (make-array size)))
    (dotimes (i (first size))
      (dotimes (j (second size))
	(setf (aref new-net i j) (v* (aref net i j) scaling))))
    (make-bspline-surface (degrees surface)
			  (knot-vectors surface)
			  new-net)))

(defun angle (v w)
  (acos (scalar-product (vnormalize v) (vnormalize w))))



;;; Call out to sfview
(asdf:oos 'asdf:load-op 'trivial-shell)
(defun sfview (surface-or-list)
  (write-rbn surface-or-list "/tmp/ntest-sfview.rbn")
  (trivial-shell:shell-command "sfview /tmp/ntest-sfview.rbn")
  t)
