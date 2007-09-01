;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Input/Output

(defgeneric write-rbn (obj filename &key if-exists)
  (:documentation "Writes OBJ to FILENAME in RBN format.
Uses the IF-EXISTS parameter as in WITH-OPEN-FILE."))

(defmethod write-rbn ((curve bspline-curve) filename
		      &key (if-exists :supersede))
  (with-open-file (s filename :direction :output :if-exists if-exists)
    (with-accessors ((degree degree)
		     (knots knot-vector)
		     (points control-points))
	curve
      (format s "(:bspline-curve~
		 ~% :degree ~d~
		 ~% :knot-vector ~a~
		 ~% :control-points ~a)~%"
	      degree (coerce knots 'list) (coerce points 'list)))))

(defmethod write-rbn ((surface bspline-surface) filename
		      &key (if-exists :supersede))
  (with-open-file (s filename :direction :output :if-exists if-exists)
    (with-accessors ((degrees degrees)
		     (knots knot-vectors)
		     (net control-net))
	surface
      (format s "(:bspline-surface~
                 ~% :degrees ~a~
                 ~% :knot-vectors (~a~%                ~a)~
                 ~% :control-net ~a)~%"
	      degrees
	      (coerce (first knots) 'list) (coerce (second knots) 'list)
	      (iter (for i from 0 below (array-dimension net 0))
		    (collect (iter (for j from 0 below (array-dimension net 1))
				   (collect (aref net i j)))))))))

(defgeneric create-object (type parameters)
  (:documentation "Creates an object of type TYPE with the given PARAMETERS."))

(defmethod create-object ((type (eql :bspline-curve)) parameters)
  (make-bspline-curve (getf parameters :degree)
		      (getf parameters :knot-vector)
		      (getf parameters :control-points)))

(defmethod create-object ((type (eql :bspline-surface)) parameters)
  (make-bspline-surface (getf parameters :degrees)
			(getf parameters :knot-vectors)
			(getf parameters :control-net)))

(defun read-rbn (filename)
  "Reads all objects in the RBN file designated by FILENAME and returns
a list of newly created objects."
  (with-open-file (s filename)
    (iter (for obj = (read s nil nil))
	  (while obj)
	  (collect (create-object (first obj) (rest obj))))))

(defun write-bss (surface filename)
  'TODO)

(define-constant +ps-width+ 595 "A4 paper width.")
(define-constant +ps-height+ 841 "A4 paper height.")

(defparameter *ps-margin* 20
  "Margin of the paper used by WRITE-PS.")

(defparameter *ps-comb-scale* 40
  "Scale factor of the curvature comb used by WRITE-PS.")

(defun write-ps-control-points (curve stream convert)
  (let ((points (control-points curve)))
    (format stream "newpath~% 1.0 0.0 0.0 setrgbcolor~%")
    (iter (for i from 0 below (length points))
	  (for movement first "moveto" then "lineto")
	  (for p = (funcall convert (elt points i)))
	  (format stream " ~f ~f ~a~
		          ~% -2 -2 rmoveto~% 4 0 rlineto~
			  ~% 0 4 rlineto~% -4 0 rlineto~
			  ~% 0 -4 rlineto~% 2 2 rmoveto~%"
		  (first p) (second p) movement))
    (format stream "stroke~%")))

(defun write-ps-curve (curve stream convert resolution)
  (format stream "newpath~% 0.0 0.0 0.0 setrgbcolor~%")
  (iter (with upper = (bsc-upper-parameter curve))
	(with lower = (bsc-lower-parameter curve))
	(with step  = (/ (- upper lower) (1- resolution)))
	(for i from 0 below resolution)
	(for movement first "moveto" then "lineto")
	(for p = (funcall convert (bsc-evaluate curve (+ lower (* step i)))))
	(format stream " ~f ~f ~a~%" (first p) (second p) movement))
  (format stream "stroke~%"))

(defun write-ps-curvature-comb (curve stream convert resolution
				&optional (comb-scale *ps-comb-scale*))
  (format stream "newpath~% 0.0 0.0 0.0 setrgbcolor~%")
  (iter (with upper = (bsc-upper-parameter curve))
	(with lower = (bsc-lower-parameter curve))
	(with step  = (/ (- upper lower) (1- resolution)))
	(for i from 0 below resolution)
	(for u = (+ lower (* step i)))
	(for begin = (funcall convert (bsc-evaluate curve u)))
	(for end = (v+ begin (v* (bsc-2d-normal curve u)
				 (abs (bsc-curvature curve u)) comb-scale)))
	(format stream " ~f ~f moveto~% ~f ~f lineto~% ~f ~f moveto~%"
		(first begin) (second begin) (first end)
		(second end) (first begin) (second begin)))
  (format stream "stroke~%"))

(defun write-ps-target (curve stream convert parameters target
			&optional (comb-scale *ps-comb-scale*))
  (format stream "newpath~% 0.0 1.0 0.0 setrgbcolor~%")
  (iter (for i from 0 below (length parameters))
	(for begin = (funcall convert (bsc-evaluate curve (elt parameters i))))
	(for u = (elt parameters i))
	(for normal = (bsc-2d-normal curve u))
	(for curvature = (bsc-curvature curve u))
	(for targeted = (elt target i))
	(for end = (v+ begin (v* normal targeted comb-scale
				 (if (> (* curvature targeted) 0) 1 -1))))
	(for command first "moveto" then "lineto")
	(format stream " ~f ~f ~a~%" (first end) (second end) command))
  (format stream "stroke~%"))

(defun write-ps-polyline (stream convert points color)
  (format stream "newpath~% ~a setrgbcolor~%" color)
  (iter (for i from 0 below (length points))
	(for p = (funcall convert (elt points i)))
	(for movement first "moveto" then "lineto")
	(format stream " ~f ~f ~a~%" (first p) (second p) movement))
  (format stream "stroke~%"))

(defun write-ps (curve filename resolution &key (control-points t)
		 (bspline t) (curvature-comb t) target-curvature start-curvature
		 end-curvature (iteration 100) left-right distance blended
		 scaling margin (comb-scale *ps-comb-scale*))
  (let* ((bbox (bsc-bounding-box curve))
	 (width (- (caadr bbox) (caar bbox)))
	 (height (- (cadadr bbox) (cadar bbox)))
	 (lower-left (car bbox))
	 (margin (or margin (list *ps-margin* *ps-margin*)))
	 (scaling (or scaling (min (/ (- +ps-width+ (* (first margin) 2))
				      width)
				   (/ (- +ps-height+ (* (second margin) 2))
				      height))))
	 (margin (list (first margin) (- +ps-height+ (second margin)
					 (* height scaling))))
	 parameters target left right)
    (with-open-file (s filename :direction :output :if-exists :supersede)
      (flet ((convert (p) (v+ (v* (v- p lower-left) scaling) margin)))
	(format s "%!PS~%")
	(when control-points
	  (write-ps-control-points curve s #'convert))
	(when bspline
	  (write-ps-curve curve s #'convert resolution))
	(when curvature-comb
	  (write-ps-curvature-comb curve s #'convert resolution comb-scale))
	(when (or target-curvature left-right blended)
	  (let ((lower (bsc-lower-parameter curve))
		(upper (bsc-upper-parameter curve)))
	    (setf parameters (arc-length-sampling curve lower upper resolution)
		  target (target-curvature curve parameters iteration
					   :start-value start-curvature
					   :end-value end-curvature))))
	(when target-curvature
	  (write-ps-target curve s #'convert parameters target comb-scale))
	(when (or left-right blended)
	  (unless distance (error "Distance must be supplied."))
	  (setf left (integrate curve parameters target
				distance :from-right nil)
		right (integrate curve parameters target
				 distance :from-right t)))
	(when left-right
	  (write-ps-polyline s #'convert left "0.0 0.0 1.0")
	  (write-ps-polyline s #'convert right "0.0 0.0 1.0"))
	(when blended
	  (write-ps-polyline s #'convert (blend-points left (nreverse right))
			     "1.0 0.0 1.0"))
	(format s "showpage~%")))))

(defgeneric write-vtk (obj filename resolution)
  (:documentation "Writes OBJ at RESOLUTION to FILENAME in VTK format."))

(defmethod write-vtk ((curve bspline-curve) filename resolution)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let* ((points (control-points curve))
	   (n (length points)))
      (format s "# vtk DataFile Version 1.0~
		 ~%B-spline Curve~
		 ~%ASCII~
		 ~%DATASET POLYDATA~%~
		 ~%POINTS ~d float~%" (+ n resolution))
      (iter (for i from 0 below n)
	    (format s "~{~f ~}~%" (2d->3d (elt points i))))
      (iter (with upper = (bsc-upper-parameter curve))
	    (with lower = (bsc-lower-parameter curve))
	    (with step  = (/ (- upper lower) (1- resolution)))
	    (for i from 0 below resolution)
	    (for point = (bsc-evaluate curve (+ lower (* step i))))
	    (format s "~{~f ~}~%" (2d->3d point)))
      (let ((nlines (+ n resolution (- 2))))
	(format s "~%LINES ~d ~d~%" nlines (* 3 nlines)))
      (iter (for i from 0 below (1- n))
	    (format s "2 ~d ~d~%" i (1+ i)))
      (iter (for i from n below (+ n (1- resolution)))
	    (format s "2 ~d ~d~%" i (1+ i))))))

(defmethod write-vtk ((surface bspline-surface) filename resolution)
  (unless (listp resolution)
    (setf resolution (list resolution resolution)))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "# vtk DataFile Version 1.0~
               ~%B-Spline Surface~
               ~%ASCII~
               ~%DATASET POLYDATA~%~
               ~%POINTS ~d float~%" (* (first resolution) (second resolution)))
    (let* ((lower (bss-lower-parameter surface))
	   (upper (bss-upper-parameter surface)))
      (dotimes (j (second resolution))
	(dotimes (i (first resolution))
	  (let ((u (interpolate (first lower)
				(/ i (1- (first resolution)))
				(first upper)))
		(v (interpolate (second lower)
				(/ j (1- (second resolution)))
				(second upper))))
	    (format s "~{~f ~}~%" (bss-evaluate surface (list u v)))))))
    (let ((u-res (first resolution))
	  (v-res (second resolution)))
      (format s "~%POLYGONS ~d ~d~%"
	      (* (1- u-res) (1- v-res)) (* (1- u-res) (1- v-res) 5))
      (dotimes (j (1- v-res))
	(dotimes (i (1- u-res))
	  (format s "4 ~d ~d ~d ~d~%"
		  (+ (* j u-res) i) (+ (* j u-res) i 1)
		  (+ (* j u-res) i u-res 1) (+ (* j u-res) i u-res)))))))

(defun write-pts (surface filename resolution)
  (unless (listp resolution)
    (setf resolution (list resolution resolution)))
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (format s "~{~d~^ ~}~%" resolution)
    (let* ((lower (bss-lower-parameter surface))
	   (upper (bss-upper-parameter surface))
	   (step (mapcar #'(lambda (x y z) (/ (- x y) (1- z)))
			 upper lower resolution)))
      (dotimes (j (second resolution))
	(dotimes (i (first resolution))
	  (let ((u (+ (first lower) (* (first step) i)))
		(v (+ (second lower) (* (second step) j))))
	    (format s "~{~f~^ ~}~%" (bss-evaluate surface (list u v)))))))))

(defun write-points2-pts (points filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let ((dim (array-dimensions points)))
      (format s "~{~d~^ ~}~%" dim)
      (dotimes (j (second dim))
	(dotimes (i (first dim))
	  (format s "~{~f~^ ~}~%" (aref points i j)))))))

(defun write-points-vtk (points filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let ((n (length points)))
      (format s "# vtk DataFile Version 1.0~
		 ~%B-spline Curve~
		 ~%ASCII~
		 ~%DATASET POLYDATA~%
		 ~%POINTS ~d float~%" n)
      (dotimes (i n)
	(format s "~{~f ~}~%" (2d->3d (elt points i))))
      (format s "~%LINES ~d ~d~%" (1- n) (* 3 (1- n)))
      (dotimes (i (1- n))
	(format s "2 ~d ~d~%" i (1+ i))))))

(defun write-points2-vtk (points filename)
  (with-open-file (s filename :direction :output :if-exists :supersede)
    (let ((n1 (array-dimension points 0))
	  (n2 (array-dimension points 1)))
      (format s "# vtk DataFile Version 1.0~
		 ~%B-spline Surface~
		 ~%ASCII~
		 ~%DATASET POLYDATA~%
		 ~%POINTS ~d float~%" (* n1 n2))
      (dotimes (i n1)
	(dotimes (j n2)
	  (format s "~{~f ~}~%" (aref points i j))))
      (format s "~%POLYGONS ~d ~d~%" (* (1- n1) (1- n2)) (* (1- n1) (1- n2) 5))
      (dotimes (i (1- n1))
	(dotimes (j (1- n2))
	  (format s "4 ~d ~d ~d ~d~%"
		  (+ (* i n2) j) (+ (* i n2) j 1)
		  (+ (* i n2) j n2 1) (+ (* i n2) j n2)))))))

