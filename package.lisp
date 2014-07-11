;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-user)

(defpackage :cl-nurbs
  #+fff (:use :common-lisp :iterate :cffi :fff)
  #-fff (:use :common-lisp :iterate)
  (:export :vlength :vlength2 :v+ :v- :v* :vnormalize :scalar-product
	   :cross-product :point-distance :affine-combine :interpolate
           ;; B-Spline Curves
	   :bspline-curve :make-bspline-curve :copy-bspline-curve
	   :degree :knot-vector :control-points
	   :bsc-dimension
	   :bsc-lower-parameter :bsc-upper-parameter
	   :bsc-bounding-box :bsc-bounding-box-axis
	   :bsc-evaluate :bsc-evaluate-on-parameters
	   :bsc-2d-normal :bsc-out-direction :bsc-out-direction-on-parameters
	   :bsc-curvature :bsc-torsion
	   :bsc-insert-knot
	   :bsc-split-curve :bsc-subcurve
	   :bsc-estimate-arc-length :bsc-iterative-arc-length
	   :bsc-knot-removal-reinsertion
	   :bsc-faired-polygon
	   :bsc-fair
	   :bsc-fit
	   :bsc-continuous-point :bsc-continuous-curve
	   :bsc-extrude
	   ;; B-Spline Surfaces
	   :bspline-surface :make-bspline-surface :copy-bspline-surface
	   :degrees :knot-vectors :control-net
	   :bss-lower-parameter :bss-upper-parameter
	   :bss-bounding-box
	   :bss-evaluate
	   :bss-surface-normal
	   :bss-principal-curvatures :bss-gaussian-curvature :bss-mean-curvature
	   :bss-insert-knot
	   :bss-split-surface
	   :bss-subsurface-one-direction
	   :bss-subsurface
	   :bss-knot-removal-reinsertion
	   :bss-construction-curve :bss-get-surface-curve
	   :bss-faired-mesh
	   :bss-fit
	   :bss-continuous-surface
	   ;; I/O
	   :write-rbn :read-rbn
	   :write-bss
	   :write-ps
	   :write-pts :read-pts
	   :write-points2-pts
	   :write-vtk
	   :write-points-vtk
	   :write-points2-vtk))
