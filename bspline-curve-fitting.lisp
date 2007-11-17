;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-nurbs)

;; Fitting

#+fff
(defun sequence->double-array (seq)
  "FOREIGN-ALLOCates an array of doubles and copies the elements of SEQ.
The returned value should be FOREIGN-FREEd after use."
  (let* ((n (length seq))
	 (array (foreign-alloc :double :count n)))
    (dotimes (i n)
      (setf (mem-aref array :double i) (coerce (elt seq i) 'double-float)))
    array))

#+fff
(defmacro with-double-arrays (binds &body body)
  "Binds variables to array representations of the given sequences.
BINDS must be of the form ((ARRAY1 SEQUENCE1) (ARRAY2 SEQUENCE2) ...).
The arrays are freed in the end."
  `(let ,(mapcar #'(lambda (x)
		     `(,(first x) (sequence->double-array ,(second x))))
		 binds)
    (unwind-protect
	 (progn ,@body)
      ,@(mapcar #'(lambda (x) `(foreign-free ,(first x))) binds))))

#+fff
(defun double-array->vector (array n)
  "Returns a vector of N elements, containing the numbers in ARRAY."
  (let ((result (make-array n)))
    (dotimes (i n)
      (setf (elt result i) (mem-aref array :double i)))
    result))

#+fff
(defun dimensionate (seq d)
  "Groups the elements of SEQ in a vector of D-long lists."
  (let* ((n (/ (length seq) d))
	 (result (make-array n)))
    (dotimes (i n)
      (setf (elt result i) (coerce (subseq seq (* i d) (* (1+ i) d)) 'list)))
    result))

#+fff
(defun bspline-curve-from-gcf (gcf)
  "Extracts the fitted B-spline curve from a GCF object."
  (let* ((nr-knots (gcf-get-nr-knots gcf))
	 (dimension (gcf-dim gcf))
	 (np (* (gcf-get-nr-ctrl-points gcf) dimension)))
    (make-bspline-curve (gcf-get-degree gcf)
			(with-foreign-object (knots :double nr-knots)
			  (gcf-get-knot-vector gcf knots)
			  (double-array->vector knots nr-knots))
			(dimensionate (with-foreign-object (points :double np)
					(gcf-get-control-points gcf points)
					(double-array->vector points np))
				      dimension))))

(defun uniform-parameter-points (points &optional (start 0.0) (end 1.0))
  "Returns a list of the form (P1 X1 Y1... P2 X2 Y2... P3...),
where P1...PN are equidistant points of the [START, END] interval."
  (let ((n (length points))
	(len (- end start)))
    (iter (for i from 0 below n)
	  (for ppts = (append ppts (list (+ (/ (* len i) (1- n)) start))
			      (copy-list (elt points i))))
	  (finally (return ppts)))))

#+fff
(defun bsc-fit (curve points tolerance)
  "Approximates POINTS with a curve, while retaining the features of CURVE.
POINTS should be a sequence of points."
  (with-accessors ((degree degree)
		   (knots knot-vector)
                   (control-points control-points))
      curve
    (let ((resolution (length points))
          (parameter-points (uniform-parameter-points points))
          (gcf (gcf-create (bsc-dimension curve)))
          result)
      (with-double-arrays ((ppoint-array parameter-points)
                           (knot-array knots)
                           (start-point (elt control-points 0))
                           (end-point (elt control-points
                                           (1- (length control-points)))))
        (gcf-add-point-group gcf (coerce tolerance 'double-float)
                             resolution ppoint-array)
        (gcf-set-degree gcf degree)
        (gcf-set-closed gcf nil)
        (gcf-set-smoothness-functional gcf :smf-crv)
        (gcf-set-optimize-parameters gcf nil)
        (gcf-set-knot-vector gcf (length knots) knot-array)
        (gcf-set-start-point gcf start-point)
        (gcf-set-end-point gcf end-point)
        (unwind-protect
             (setf result (if (eql (gcf-fit gcf) :success)
                              (bspline-curve-from-gcf gcf)
                              nil))
          (gcf-destroy gcf)))
      result)))

#-fff
(defun bsc-fit (curve points tolerance)
  "Approximates POINTS with a curve, while retaining the features of CURVE.
POINTS should be a sequence of points."
  'TODO)
