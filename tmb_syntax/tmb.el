;;; tmb.el --- Major mode for creating statistical models with TMB

;; Copyright (C) 2015 Arni Magnusson

;; Author:   Arni Magnusson
;; Keywords: languages

(defconst tmb-mode-version "1.0" "TMB Mode version number.")

;; This file is not part of GNU Emacs.

;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:
;;
;; Major mode for editing Template Model Builder (TMB) code, derived from
;; `c++-mode'. Provides syntax highlighting. The syntax groups for highlighting
;; are:
;;
;; Face                          Example
;; tmb-data-face                 DATA_VECTOR
;; tmb-parameter-face            PARAMETER
;; tmb-report-face               ADREPORT
;; font-lock-builtin-face        error
;; font-lock-comment-face        //
;; font-lock-constant-face       isNumeric
;; font-lock-function-name-face  [Type] myfunction
;; font-lock-keyword-face        log
;; font-lock-type-face           Type
;;
;; Installation:
;;
;; 1. Copy this file (tmb.el) to a directory in `load-path', or edit .emacs to
;;    add the directory to `load-path':
;;      (add-to-list 'load-path "mypath/tmb")
;; 2. Byte-compile this file to tmb.elc for faster startup:
;;      M-x byte-compile-file
;; 3. Edit .emacs so that `tmb-mode' is loaded at startup:
;;      (require 'tmb)
;;
;; Customization:
;;
;; If you want to set syntax colors, here is an example that does that:
;; (defun my-tmb-hook ()
;;   (set-face-attribute 'tmb-data-face      nil :foreground "dodgerblue")
;;   (set-face-attribute 'tmb-parameter-face nil :foreground "dodgerblue")
;;   (set-face-attribute 'tmb-report-face    nil :foreground "dodgerblue")
;;   (set-face-attribute 'font-lock-variable-name-face nil
;;                       :foreground 'unspecified))
;; (add-hook 'tmb-mode-hook 'my-tmb-hook)
;;
;; Usage:
;;
;; See documentation of function `tmb-mode'.
;;
;; Known issues:
;;
;; Cursor motion and deletions swallow entire underscore_separated_object_name,
;; instead of pausing at each underscore.

;;; History:
;;
;; 01 Sep 2015  1.0  Created main function `tmb-mode', derived from `c++-mode'.

;;; Code:

;; 1  Preamble

(require 'cc-mode)
(defgroup tmb nil "Major mode for editing Template Model Builder code."
  :tag "TMB" :group 'languages)
;; Switch to `tmb-mode' if C++ buffer includes TMB header
(setq magic-mode-regexp-match-limit 40000)
(defun tmb-include-p () "Search for #include <TMB.hpp> in C++ file."
       (if (string-equal (file-name-extension (buffer-name)) "cpp")
           (re-search-forward "#include ?<TMB.hpp>"
                              magic-mode-regexp-match-limit t)))
(setq magic-mode-alist '((tmb-include-p . tmb-mode)))

;; 2  User variables

(defface tmb-data-face '((t :inherit font-lock-type-face))
  "Font Lock face to highlight TMB data macros." :group 'tmb)
(defvar tmb-data-face 'tmb-data-face "Face name for TMB data macros.")
(defface tmb-parameter-face '((t :inherit font-lock-type-face))
  "Font Lock face to highlight TMB parameter macros." :group 'tmb)
(defvar tmb-parameter-face 'tmb-parameter-face
  "Face name for TMB parameter macros.")
(defface tmb-report-face '((t :inherit font-lock-keyword-face))
  "Font Lock face to highlight TMB report macros." :group 'tmb)
(defvar tmb-report-face 'tmb-report-face "Face name for TMB report macros.")

;; 3  Internal variables

(defvar tmb-font-lock-keywords
  (eval-when-compile
    (let ((TYPE '("scalartype" "Type"))
          (DATA
           '("DATA_INTEGER" "DATA_IVECTOR" "DATA_IARRAY"
             "DATA_SCALAR"  "DATA_VECTOR"  "DATA_ARRAY"
             "DATA_FACTOR"  "DATA_STRUCT"
             "DATA_MATRIX"  "DATA_SPARSE_MATRIX"
             "DATA_VECTOR_INDICATOR" "DATA_ARRAY_INDICATOR"
             "TMB_ATOMIC_VECTOR_FUNCTION"))
          (PARAMETERS
           '("PARAMETER" "PARAMETER_VECTOR"
             "PARAMETER_ARRAY" "PARAMETER_MATRIX"))
          (REPORT '("ADREPORT" "REPORT"))
          (FUNCTIONS
           '(;; C
             "free" "malloc"
             ;; C++
             "clear" "cout" "erase"
             ;; AD
             "ForTwo" "Forward" "myReverse" "Reverse" "RevTwo"
             "Dependent" "Independent" "Value" "Variable"
             "ForSparseJac" "Hessian" "Jacobian" "RevSparseHes"
             "generalized_symbol" "optimize" "tape_symbol" "traceforward0sweep"
             "Domain" "Range"
             "CondExpEq" "CondExpGe" "CondExpGt" "CondExpLt"
             "REGISTER_ATOMIC" "SEPARABLE"
             ;; OpenMP
             "omp_get_max_threads" "omp_get_thread_num" "start_parallel"
             ;; Rcpp
             "allocVector" "asSEXP" "defineVar" "install" "mkChar" "set"
             "getAttrib" "setAttrib"
             "R_FlushConsole" "Rprintf" "REprintf"
             "LENGTH"
             "SET_STRING_ELT" "SET_VECTOR_ELT" "STRING_ELT" "VECTOR_ELT"
             "PROTECT" "UNPROTECT"
             "R_do_slot_assign" "R_ExternalPtrAddr" "R_ExternalPtrTag"
             "R_MakeExternalPtr" "R_RunExitFinalizers"
             ;; TMB housekeeping
             "conservativeResize" "increase" "pushParname"
             "evalUserTemplate"
             "KEEP_COL" "KEEP_ROW"
             "CallCFinalizer" "R_RegisterCFinalizer" "RegisterCFinalizer"
             "HessianSparsityPattern" "optimizeTape"
             "count_parallel_regions" "set_parallel_region"
             ;; Basic math
             "abs" "exp" "log" "pow" "sqrt"
             "ceil" "floor" "max" "min" "trunc"
             "mod"
             "acos" "asin" "atan" "cos" "cosh" "sin" "sinh" "tan" "tanh"
             "besselK" "invlogit"
             ;; Is what
             "isEnvironment" "isMatrix" "isNewList" "isReal"
             "R_IsNA" "RObjectTestExpectedType"
             ;; Cast
             "asDouble" "asMatrix" "asVector" "coerceVector"
             "CHAR" "INTEGER" "Integer" "REAL" "real"
             "mat2vec" "vec2mat"
             ;; Vectors and matrices, basics
             "array" "c" "matrix" "tuple" "value" "vec"
             "begin" "end" "head" "index" "reverse" "segment" "tail"
             "dim" "length" "resize" "push" "setdim" "size" "squeeze"
             "ncols" "nrows" "outerSize" "size1" "size2"
             "col" "cols" "colwise" "row" "rows" "rowwise"
             "asDiagonal" "block"  "diagonal"
             "init" "initZeroArray" "setIdentity" "setZero"
             "maxCoeff" "prod" "sum"
             "perm"
             "fill" "fillShape" "getShape"
             ;; Vectors and matrices, linear algebra
             "determinant" "inverse" "norm" "solve" "trace" "transpose"
             "rotate" "scale"
             "expm" "matmul" "matmult"
             "invpd" "logdet" "matinv" "matinvpd"
             "cov" "UNSTRUCTURED_CORR"
             "eigenvalues" "eigenvectors"
             "jacobian"
             ;; Time series
             "AR1"
             ;; Lists
             "getListElement"
             ;; Factors
             "NLEVELS"
             ;; Distributions
             "dnorm" "pnorm" "qnorm"
             "dnorm1" "pnorm1" "qnorm1"
             "pnorm_approx" "qnorm_approx"
             "dt"
             "df"
             "MVNORM"
             "dbinom"
             "dzinbinom" "dzinbinom2"
             "dmultinom"
             "dpois" "ppois"
             "dzipois"
             "dnbinom" "dnbinom2"
             "dexp" "pexp" "qexp"
             "dbeta"
             "dgamma" "dlgamma" "lgamma" "lgammafn" "pgamma" "psigamma" "qgamma"
             "dweibull" "qweibull"
             "SCALE" "VECSCALE"))
          (CONSTANTS '("isNumeric" "isNumericScalar"))
          (IMPORTANT '("error" "TMB_CATCH" "TMB_ERROR_BAD_ALLOC" "TMB_TRY")))
      (list
       (cons (regexp-opt TYPE 'words) font-lock-type-face)
       (cons (regexp-opt DATA 'words) tmb-data-face)
       (cons (regexp-opt PARAMETERS 'words) tmb-parameter-face)
       (cons (regexp-opt REPORT 'words) tmb-report-face)
       (cons (regexp-opt FUNCTIONS 'words) font-lock-keyword-face)
       (cons (regexp-opt CONSTANTS 'words) font-lock-constant-face)
       (cons (regexp-opt IMPORTANT 'words) font-lock-builtin-face)))))
(nconc tmb-font-lock-keywords c++-font-lock-keywords)

;; 4  User functions

(defun tmb-help () "Show help page for `tmb-mode'." (interactive)
       (describe-function 'tmb-mode)(switch-to-buffer "*Help*")
       (delete-other-windows)(message nil))
(defun tmb-mode-version () "Show TMB Mode version number." (interactive)
       (message (concat "TMB Mode version " tmb-mode-version)))

;; 5  Main function

;;;###autoload
(define-derived-mode tmb-mode c++-mode "TMB"
  "Major mode for creating statistical models with TMB.\n
The `tmb-help' command shows this page.\n
\\{tmb-mode-map}"
  (set (make-local-variable 'font-lock-defaults)
       '(tmb-font-lock-keywords nil nil)))
(modify-syntax-entry ?_ "w" tmb-mode-syntax-table)

(provide 'tmb)
