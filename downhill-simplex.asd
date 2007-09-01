;; -*- mode: lisp; syntax: common-lisp -*-

(in-package :cl-user)

(defpackage :downhill-simplex-asd
  (:use :cl :asdf))

(in-package :downhill-simplex-asd)

(asdf:defsystem :downhill-simplex
  :author "Peter Salvi"
  :licence "Public Domain"
  :description "Downhill Simplex Method as in Numerical Recipes"
  :components ((:file "downhill-simplex")))
