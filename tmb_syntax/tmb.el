;;; tmb.el --- Major mode for creating statistical models with TMB

;; Copyright (C) 2015 Arni Magnusson

;; Author:   Arni Magnusson
;; Keywords: languages

(defconst tmb-mode-version "1.2" "TMB Mode version number.")

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
;; `c++-mode'. Provides syntax highlighting, IDE compilation, file manipulation,
;; and smaller tools. The syntax groups for highlighting are:
;;
;; Face                          Example
;; tmb-data-face                 DATA_VECTOR
;; tmb-parameter-face            PARAMETER
;; tmb-report-face               ADREPORT
;; font-lock-builtin-face        error
;; font-lock-comment-face        //
;; font-lock-constant-face       dim
;; font-lock-function-name-face  [Type] myfunction
;; font-lock-keyword-face        log
;; font-lock-type-face           Type
;;
;; Keybindings defined in `tmb-mode' are overridden by user `c-mode-common-hook'
;; and `c++-mode-hook', which can then be overridden by user `tmb-mode-hook'.
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
;; If you want to set the default R command, syntax colors, or keybindings, here
;; is an example that does that:
;; (defun my-tmb-hook ()
;;   (setq tmb-r-command "R --slave <")
;;   (set-face-attribute 'tmb-data-face      nil :foreground "dodgerblue")
;;   (set-face-attribute 'tmb-parameter-face nil :foreground "dodgerblue")
;;   (set-face-attribute 'tmb-report-face    nil :foreground "dodgerblue")
;;   (set-face-attribute 'font-lock-variable-name-face nil
;;                       :foreground 'unspecified)
;;   (local-set-key [f9]  'tmb-run-r  )
;;   (local-set-key [f10] 'tmb-open-r))
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
;; 04 Sep 2015  1.2  Added user functions `tmb-clean', `tmb-for',
;;                   `tmb-kill-process', `tmb-open', `tmb-open-r',
;;                   `tmb-run-any', `tmb-run-make', `tmb-run-r', and
;;                   `tmb-toggle-section'. Added user variables
;;                   `tmb-make-command' and `tmb-r-command'. Disabled
;;                   `abbrev-mode'.
;; 03 Sep 2015  1.1  Shortened list of recognized FUNCTIONS for maintainability.
;; 01 Sep 2015  1.0  Created main function `tmb-mode', derived from `c++-mode'.

;;; Code:

;; 1  Preamble

(require 'cc-mode)
(defgroup tmb nil
  "Major mode for editing Template Model Builder code."
  :tag "TMB" :group 'languages)
;; Switch to `tmb-mode' if C++ buffer includes TMB header
(setq magic-mode-regexp-match-limit 40000)
(defun tmb-include-p ()
  "Search for #include <TMB.hpp> in C++ file."
  (if (string-equal (file-name-extension (buffer-name)) "cpp")
      (re-search-forward "#include ?<TMB.hpp>"
                         magic-mode-regexp-match-limit t)))
(setq magic-mode-alist '((tmb-include-p . tmb-mode)))

;; 2  User variables

(defcustom tmb-make-command "make"
  "Shell command to run makefile using `tmb-run-make'."
  :tag "Make" :type 'string :group 'tmb)
(defcustom tmb-r-command "R --quiet --vanilla <"
  "Shell command to run R script using `tmb-run-r'."
  :tag "R" :type 'string :group 'tmb)
(defface tmb-data-face '((t :inherit font-lock-type-face))
  "Font Lock face to highlight TMB data macros." :group 'tmb)
(defvar tmb-data-face 'tmb-data-face
  "Face name for TMB data macros.")
(defface tmb-parameter-face '((t :inherit font-lock-type-face))
  "Font Lock face to highlight TMB parameter macros." :group 'tmb)
(defvar tmb-parameter-face 'tmb-parameter-face
  "Face name for TMB parameter macros.")
(defface tmb-report-face '((t :inherit font-lock-keyword-face))
  "Font Lock face to highlight TMB report macros." :group 'tmb)
(defvar tmb-report-face 'tmb-report-face
  "Face name for TMB report macros.")

;; 3  Internal variables

(defvar tmb-font-lock-keywords
  (eval-when-compile
    (let ((TYPE '("Integer" "Type"))
          (DATA
           '("DATA_INTEGER" "DATA_IVECTOR" "DATA_IARRAY"
             "DATA_SCALAR"  "DATA_VECTOR"  "DATA_ARRAY"
             "DATA_FACTOR"  "DATA_STRUCT"
             "DATA_MATRIX"  "DATA_SPARSE_MATRIX"
             "DATA_VECTOR_INDICATOR" "DATA_ARRAY_INDICATOR"))
          (PARAMETERS
           '("PARAMETER" "PARAMETER_VECTOR"
             "PARAMETER_ARRAY" "PARAMETER_MATRIX"))
          (REPORT '("ADREPORT" "REPORT"))
          (FUNCTIONS
           '(;; I/O
             "cout" "endl"
             ;; Basic math
             "abs" "exp" "log" "pow" "sqrt"
             "ceil" "floor" "mod" "trunc"
             "acos" "asin" "atan" "cos" "cosh" "sin" "sinh" "tan" "tanh"
             "invlogit" "lgamma" "logit"
             ;; Arrays, basics
             "array" "matrix"
             "head" "segment" "size" "tail"
             "col" "cols" "colwise" "row" "rows" "rowwise"
             "block" "diagonal" "transpose"
             "max" "min" "prod" "sum"
             "maxCoeff" "minCoeff"
             "setIdentity" "setZero"
             ;; Arrays, linear algebra
             "determinant" "expm" "inverse" "norm" "trace"
             ;; Distributions
             "dnorm" "pnorm" "qnorm"
             "dnorm1" "pnorm1" "qnorm1"
             "pnorm_approx" "qnorm_approx"
             "AR1" "MVNORM"
             "dt" "df"
             "dbinom" "dmultinom"
             "dpois" "ppois" "dzipois"
             "dnbinom" "dnbinom2" "dzinbinom" "dzinbinom2"
             "dbeta"
             "dexp" "pexp" "qexp"
             "dgamma" "dlgamma" "pgamma" "qgamma"
             "dweibull" "pweibull" "qweibull"))
          (CONSTANTS '("dim"))
          (WARNINGS '("error")))
      (list
       (cons (regexp-opt TYPE 'words) font-lock-type-face)
       (cons (regexp-opt DATA 'words) 'tmb-data-face)
       (cons (regexp-opt PARAMETERS 'words) 'tmb-parameter-face)
       (cons (regexp-opt REPORT 'words) 'tmb-report-face)
       (cons (regexp-opt FUNCTIONS 'words) font-lock-keyword-face)
       (cons (regexp-opt CONSTANTS 'words) font-lock-constant-face)
       (cons (regexp-opt WARNINGS 'words) font-lock-warning-face)))))
(nconc tmb-font-lock-keywords c++-font-lock-keywords)

;; 4  User functions

(defun tmb-clean ()
  "Remove C++ binary files (*.o *.so *.dll)." (interactive)
  (let* ((prog (file-name-sans-extension (buffer-name)))
         (pattern (concat prog "\\.o\\|" prog "\\.so\\|" prog "\\.dll"))
         (files (directory-files "." nil pattern t)))
    (dolist (x files)(delete-file x)))(message "Removed binary files."))
(defun tmb-for ()
  "Insert for(int i=0; i<; i++)." (interactive)
  (insert "for(int i=0; i<; i++)")(search-backward ";"))
(defun tmb-help ()
  "Show help page for `tmb-mode'." (interactive)
  (describe-function 'tmb-mode)(switch-to-buffer "*Help*")(delete-other-windows)
  (message nil))
(defun tmb-kill-process ()
  "Stop the current process, usually TMB compilation or model run."
  (interactive)(kill-process (car (process-list))))
(defun tmb-mode-version ()
  "Show TMB Mode version number." (interactive)
  (message (concat "TMB Mode version " tmb-mode-version)))
(defun tmb-open (ext)
  "Open file with extension EXT in other window." (interactive "sExtension: ")
  (let ((file (concat (file-name-sans-extension (buffer-name)) "." ext)))
    (if (not (file-regular-p file))(error "File %s not found" file)
      (find-file-other-window file))))
(defun tmb-open-r ()
  "Open R script with same filename prefix as current buffer." (interactive)
  (tmb-open "R"))
(defun tmb-run-any (script)
  "Run any R script, querying user for SCRIPT filename.\n
If the R script has the same filename prefix as the current buffer, then use
`tmb-run-r' instead.\n
Filename history is accessible using M-p and up arrow."
  (interactive "fRun R script: ")(save-buffer)
  (compile (concat tmb-r-command " " script))
  (with-current-buffer "*compilation*" (setq show-trailing-whitespace nil)))
(defun tmb-run-make ()
  "Run makefile in current directory, using `tmb-make-command'."
  (interactive)(save-buffer)(compile tmb-make-command)
  (with-current-buffer "*compilation*" (setq show-trailing-whitespace nil)))
(defun tmb-run-r ()
  "Run R script with same filename prefix as current buffer.\n
If the R script has a different filename, then use `tmb-run-any' instead."
  (interactive)
  (tmb-run-any (concat (file-name-sans-extension (buffer-name)) ".R")))
(defun tmb-toggle-section ()
  "Toggle whether the current section is indicated in the mode line."
  (interactive)(which-function-mode (if which-function-mode 0 1))
  (message "Section indicator %s" (if which-function-mode "ON" "OFF")))

;; 5  Main function

;;;###autoload
(define-derived-mode tmb-mode c++-mode "TMB"
  "Major mode for creating statistical models with TMB.\n
The `tmb-help' command shows this page.\n
\\{tmb-mode-map}"
  (set (make-local-variable 'font-lock-defaults)
       '(tmb-font-lock-keywords nil nil))
  (abbrev-mode 0))
(modify-syntax-entry ?_ "w" tmb-mode-syntax-table)
(define-key tmb-mode-map [?\C-c ?\C-/] 'tmb-help          )
(define-key tmb-mode-map [?\C-c ?\C-a] 'tmb-run-any       )
(define-key tmb-mode-map [?\C-c ?\C-d] 'tmb-clean         )
(define-key tmb-mode-map [?\C-c ?\C-f] 'tmb-for           )
(define-key tmb-mode-map [?\C-c ?\C-k] 'tmb-kill-process  )
(define-key tmb-mode-map [?\C-c ?\C-m] 'tmb-run-make      )
(define-key tmb-mode-map [?\C-c ?\C-o] 'tmb-open          )
(define-key tmb-mode-map [?\C-c ?\C-p] 'tmb-open-r        )
(define-key tmb-mode-map [?\C-c ?\C-r] 'tmb-run-r         )
(define-key tmb-mode-map [?\C-c ?\C-s] 'tmb-toggle-section)

(provide 'tmb)
