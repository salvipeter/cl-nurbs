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
	       (:file "utils" :depends-on ("package"))
	       (:file "affine" :depends-on ("package"))
	       (:file "bspline-curve" :depends-on ("utils" "affine"))
	       (:file "bspline-curve-fairing" :depends-on ("bspline-curve"))
	       (:file "bspline-curve-fitting" :depends-on ("bspline-curve"))
	       (:file "bspline-curve-continuity" :depends-on ("bspline-curve"))
	       (:file "bspline-surface" :depends-on ("bspline-curve"))
	       (:file "bspline-surface-fairing"
		      :depends-on ("bspline-surface"))
	       (:file "bspline-surface-fitting"
		      :depends-on ("bspline-surface"))
	       (:file "bspline-surface-continuity"
		      :depends-on ("bspline-surface"))
	       (:file "io" :depends-on ("bspline-curve-fairing"
					"bspline-surface"))))
