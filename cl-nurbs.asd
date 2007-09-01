;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-user)

(defpackage :cl-nurbs-asd
  (:use :cl :asdf))

(in-package :cl-nurbs-asd)

(asdf:defsystem :cl-nurbs
  :author "Peter Salvi"
  :licence "Public Domain"
  :description "NURBS Library"
;;; FFF is only needed for interfacing the
;;; Raindrop Geomagic FFF (proprietary) library
;;   :depends-on (:iterate :downhill-simplex)
  :depends-on (:iterate :fff :downhill-simplex)
  :components ((:file "package")
	       (:file "utils")
	       (:file "affine")
	       (:file "bspline-curve")
	       (:file "bspline-curve-fairing")
	       (:file "bspline-curve-fitting")
	       (:file "bspline-curve-continuity")
	       (:file "bspline-surface")
	       (:file "bspline-surface-fairing")
	       (:file "bspline-surface-fitting")
	       (:file "bspline-surface-continuity")
	       (:file "io")))
