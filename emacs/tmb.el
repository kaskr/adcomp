;;; tmb.el --- Major mode for creating statistical models with TMB

;; Copyright (C) 2015 Arni Magnusson

;; Author:   Arni Magnusson
;; Keywords: languages

(defconst tmb-mode-version "2.2" "TMB Mode version number.")

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
;; basic templates, and smaller tools. The syntax groups for highlighting are:
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
;; 4. In Windows, the PATH environment variable should include the directory
;;      containing the R executables (path/to/bin/x64).
;; 5. Some features of this package require that Emacs Speaks Statistics (ESS)
;;      is also installed.
;;
;; Customization:
;;
;; This mode lets Emacs scan C++ files to check if they #include <TMB.hpp> in
;; the first `magic-mode-regexp-match-limit' characters of the file. To make
;; Emacs search deeper than the first few lines of C++ files, you can increase
;; the limit in your .emacs:
;;
;; (setq magic-mode-regexp-match-limit 40000)
;;
;; If you want to set the default R command, syntax colors, or keybindings, here
;; is an example that does that:
;;
;; (defun my-tmb-hook ()
;;   (setq tmb-r-command "R --slave <")
;;   (set-face-attribute 'tmb-data-face      nil :foreground "dodgerblue")
;;   (set-face-attribute 'tmb-parameter-face nil :foreground "dodgerblue")
;;   (set-face-attribute 'tmb-report-face    nil :foreground "dodgerblue")
;;   (set-face-attribute 'font-lock-variable-name-face nil
;;                       :foreground 'unspecified)
;;   (local-set-key [f9]           'tmb-run )
;;   (local-set-key [f10]          'tmb-open)
;;   (local-set-key [down-mouse-3] 'imenu  ))
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
;; 17 Sep 2015  2.2  Added GUI menu and toolbar. Added user function
;;                   `tmb-new-buffer'. Added internal variables `tmb-menu',
;;                   `tmb-mode-map', and `tmb-tool-bar-map'. Improved
;;                   `tmb-template-mini'.
;; 10 Sep 2015  2.1  Added internal function `tmb-windows-os-p'. Improved
;;                   `tmb-run-debug' and `tmb-template-mini'.
;; 07 Sep 2015  2.0  Added user functions `tmb-run-debug', `tmb-scroll-down',
;;                   `tmb-scroll-up', `tmb-show-compilation', and `tmb-show-r'.
;;                   Renamed `tmb-open' to `tmb-open-any' and `tmb-run-r' to
;;                   `tmb-run'. Improved `tmb-open', `tmb-open-any', `tmb-run',
;;                   `tmb-run-any', and `tmb-template-mini'.
;; 05 Sep 2015  1.3  Added user function `tmb-template-mini'. Renamed
;;                   `tmb-toggle-section' to `tmb-toggle-function'.
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

(require 'cc-mode) ; c++-font-lock-keywords
(require 'compile) ; compilation-scroll-output
;; (require 'ess-site) ; slow, load only when needed
(add-to-list 'magic-mode-alist '(tmb-include-p . tmb-mode))
(declare-function ess-eval-linewise "ess-inf")
(declare-function ess-process-live-p "ess-inf")
(declare-function ess-show-buffer "ess-inf")
(defgroup tmb nil
  "Major mode for editing Template Model Builder code."
  :tag "TMB" :group 'languages)

;; 2  User variables

(defcustom tmb-make-command "make"
  "Shell command to run makefile using `tmb-run-make'."
  :tag "Make" :type 'string :group 'tmb)
(defcustom tmb-r-command "R --quiet --vanilla <"
  "Shell command to run R script using `tmb-run'."
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
(defvar tmb-menu
  '("TMB"
    ["Run"              tmb-run          ]
    ["Make"             tmb-run-make     ]
    "--"
    ["Stop"             tmb-kill-process ]
    ["Clean"            tmb-clean        ]
    ["Debug"            tmb-run-debug    ]
    "--"
    ["View Script"      tmb-open         ]
    ["Mini Template"    tmb-template-mini]
    "--"
    ["Help"             tmb-help         ]
    ["TMB Mode Version" tmb-mode-version ]))
(defvar tmb-mode-map
  ;; Don't use C-c C-                        x
  ;; Special   C-c C-        h
  ;; Custom    C-c C- a cd f    klm opqrs
  ;; Available C-c C-  b  e g ij   n     tuvw yz
  (let ((map (make-sparse-keymap)))
    (easy-menu-define nil map nil tmb-menu)
    (define-key map [f12]               'tmb-template-mini   )
    (define-key map [?\C-c C-backspace] 'tmb-clean           )
    (define-key map [M-up]              'tmb-scroll-up       )
    (define-key map [M-down]            'tmb-scroll-down     )
    (define-key map [?\C-c ?\C-.]       'tmb-mode-version    )
    (define-key map [?\C-c ?\C-/]       'tmb-help            )
    (define-key map [?\C-c ?\C-a]       'tmb-run-any         )
    (define-key map [?\C-c ?\C-c]       'tmb-run             )
    (define-key map [?\C-c ?\C-d]       'tmb-run-debug       )
    (define-key map [?\C-c ?\C-f]       'tmb-for             )
    (define-key map [?\C-c ?\C-k]       'tmb-kill-process    )
    (define-key map [?\C-c ?\C-l]       'tmb-show-compilation)
    (define-key map [?\C-c ?\C-m]       'tmb-run-make        )
    (define-key map [?\C-c ?\C-o]       'tmb-open-any        )
    (define-key map [?\C-c ?\C-p]       'tmb-open            )
    (define-key map [?\C-c ?\C-q]       'tmb-kill-process    )
    (define-key map [?\C-c ?\C-r]       'tmb-show-r          )
    (define-key map [?\C-c ?\C-s]       'tmb-toggle-function )
    (define-key map [?\C-\M-v]          'ignore              )
    map))
(defvar tmb-tool-bar-map
  (let ((tool-bar-map (make-sparse-keymap)) ; ; undo-form from menu-bar.el
        (undo-form '(and (not buffer-read-only)(not (eq t buffer-undo-list))
                         (if (eq last-command 'undo)(listp pending-undo-list)
                           (consp buffer-undo-list)))))
    (tool-bar-add-item "new" 'tmb-new-buffer 'tmb-new-buffer :help "New")
    (tool-bar-add-item "open" 'find-file 'find-file :help "Open")
    (tool-bar-add-item "save" 'save-buffer 'save-buffer :help "Save"
                       :enable '(buffer-modified-p))
    (tool-bar-add-item "cut" 'kill-region 'kill-region :help "Cut")
    (tool-bar-add-item "copy" 'copy-region-as-kill 'copy-region-as-kill
                       :help "Copy")
    (tool-bar-add-item "paste" 'yank 'yank :help "Paste")
    (tool-bar-add-item "undo" 'undo 'undo :help "Undo" :enable undo-form)
    (tool-bar-add-item "close" 'kill-this-buffer 'kill-this-buffer
                       :help "Close")
    (tool-bar-add-item "jump-to" 'tmb-run 'tmb-run :help "Run")
    tool-bar-map))

;; 4  User functions

(defun tmb-clean ()
  "Remove C++ binary files (*.o *.so *.dll)." (interactive)
  (let* ((prog (file-name-sans-extension (buffer-name)))
         (pattern (concat prog "\\.o\\|" prog "\\.so\\|" prog "\\.dll"))
         (files (directory-files "." nil pattern t)))
    (dolist (x files)(delete-file x)))(message "Removed binary files"))
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
  (message "TMB Mode version %s" tmb-mode-version))
(defun tmb-new-buffer ()
  "Create new buffer." (interactive)
  (switch-to-buffer (generate-new-buffer "Untitled"))
  (eval (list (default-value 'major-mode))))
(defun tmb-open ()
  "Open R script with same filename prefix as current buffer." (interactive)
  (tmb-open-any "R"))
(defun tmb-open-any (ext)
  "Open file with extension EXT in other window." (interactive "sExtension: ")
  (let ((file (concat (file-name-sans-extension (buffer-name)) "." ext)))
    (if (not (file-regular-p file))(error "File %s not found" file)
      (find-file-other-window file)))(other-window 1))
(defun tmb-run ()
  "Run R script with same filename prefix as current buffer.\n
If the R script has a different filename, then use `tmb-run-any' instead.\n
Navigate compilation errors with \\<tmb-mode-map>\\[tmb-scroll-down] and \
\\[tmb-scroll-up]."
  (interactive)(save-buffer)
  (tmb-run-any (concat (file-name-sans-extension (buffer-name)) ".R")))
(defun tmb-run-any (script)
  "Run any R script, querying user for SCRIPT filename.\n
If the R script has the same filename prefix as the current buffer, then use
`tmb-run' instead.\n
Filename history is accessible in the minibuffer prompt \
(\\<minibuffer-local-map>\\[previous-history-element],\
 \\[next-history-element]).\n
Navigate compilation errors with \\<tmb-mode-map>\\[tmb-scroll-down] and \
\\[tmb-scroll-up]."
  (interactive "fRun R script: ")(save-buffer)
  (compile (concat tmb-r-command " " script))
  (with-current-buffer "*compilation*" (setq show-trailing-whitespace nil)))
(defun tmb-run-debug ()
  "Debug model with GDB.\n
The R session stays alive if it was running when this function was called."
  (interactive)(save-buffer)(message "Invoking debug session...")
  (require 'ess-site)
  ;; Platform-specific: -O1 and DLLFLAGS in Windows, -O0 otherwise
  (let* ((ess-dialect "R")
         (inferior-R-args "--quiet --vanilla")
         (ess-ask-for-ess-directory nil)
         (prefix (file-name-sans-extension (buffer-name)))
         (cmd-1 "require(TMB)")
         (cmd-2 (concat "; compile(\"" prefix ".cpp\",\"-g -O"
                        (if (tmb-windows-os-p) "1\",DLLFLAGS=\"\")" "0\")")))
         (cmd-3 (concat "; gdbsource(\"" prefix ".R\",TRUE)"))
         (cmd-4 (if (ess-process-live-p) "" "; q()"))
         (cmd (concat cmd-1 cmd-2 cmd-3 cmd-4)))
    (ess-eval-linewise cmd)(message "Invoking debug session...done")))
(defun tmb-run-make ()
  "Run makefile in current directory, using `tmb-make-command'."
  (interactive)(save-buffer)(compile tmb-make-command)
  (with-current-buffer "*compilation*" (setq show-trailing-whitespace nil)))
(defun tmb-scroll-down (n)
  "Scroll other window down N lines, or visit next error message.\n
The behavior of this command depends on whether the compilation buffer is
visible."
  (interactive "p")
  (if (null (get-buffer-window "*compilation*"))(scroll-other-window n)
    (next-error n)))
(defun tmb-scroll-up (n)
  "Scroll other window up N lines, or visit previous error message.\n
The behavior of this command depends on whether the compilation buffer is
visible."
  (interactive "p")
  (if (null (get-buffer-window "*compilation*"))(scroll-other-window (- n))
    (previous-error n)))
(defun tmb-show-compilation ()
  "Show compilation buffer." (interactive)
  (if (null (get-buffer "*compilation*"))
      (error "*compilation* buffer not found")(display-buffer "*compilation*")))
(defun tmb-show-r ()
  "Show R interactive buffer." (interactive)
  (require 'ess-site)
  (if (null (get-buffer "*R*"))(error "*R* interactive buffer not found")
    (ess-show-buffer "*R*")))
(defun tmb-template-mini ()
  "Create minimal TMB files (mini.cpp, mini.R) in current directory."
  (interactive)
  ;; Platform-specific: -O1 in Windows, -O0 otherwise
  (if (file-exists-p "mini.cpp")
      (error "Error: file mini.cpp already exists in current directory"))
  (if (file-exists-p "mini.R")
      (error "Error: file mini.R already exists in current directory"))
  (if (get-buffer "mini.cpp")(error "Error: buffer mini.cpp already exists"))
  (if (get-buffer "mini.R")(error "Error: buffer mini.R already exists"))
  (delete-other-windows)(find-file "mini.cpp")(insert "\
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  PARAMETER(mu);
  PARAMETER(logSigma);

  Type f;
  f = -sum(dnorm(x, mu, exp(logSigma), true));

  return f;
}
")
  (goto-char (point-min))(save-buffer "mini.cpp")
  (find-file-other-window "mini.R")(insert "\
data <- list(x=rivers)
parameters <- list(mu=0, logSigma=0)

require(TMB)
compile(\"mini.cpp\", \"-O" (if (tmb-windows-os-p) "1" "0") " -Wall\")
dyn.load(dynlib(\"mini\"))

################################################################################

model <- MakeADFun(data, parameters)
fit <- nlminb(model$par, model$fn, model$gr)
rep <- sdreport(model)

rep
")
  (goto-char (point-min))(save-buffer "mini.R")(other-window 1)(tmb-mode)
  (message (concat "Ready to run R script ("
                   (substitute-command-keys "\\<tmb-mode-map>\\[tmb-run]")
                   ") or edit code.")))
(defun tmb-toggle-function ()
  "Toggle whether the current function is indicated in the mode line."
  (interactive)(which-function-mode (if which-function-mode 0 1))
  (message "Function indicator %s" (if which-function-mode "ON" "OFF")))

;; 5  Internal functions

(defun tmb-include-p ()
  "Check if C++ file has #include <TMB.hpp>."
  (if (string-equal (file-name-extension (buffer-name)) "cpp")
      (re-search-forward "#include ?<TMB.hpp>"
                         magic-mode-regexp-match-limit t)))
(defun tmb-windows-os-p ()
  "Check if TMB is running in a Windows operating system."
  (if (string-match "windows" (prin1-to-string system-type)) t nil))

;; 6  Main function

;;;###autoload
(define-derived-mode tmb-mode c++-mode "TMB"
  "Major mode for creating statistical models with TMB.\n
The `tmb-help' command shows this page.\n
Start a new model from `tmb-template-mini'. Navigate between functions using
`imenu' and indicate the current function in the mode line with
`tmb-toggle-function'. Prototype for-loops with `tmb-for'.\n
Use `tmb-open' to open the corresponding R script and `tmb-open-any' to open
other model-related files. Show interactive *compilation* and *R* buffers with
`tmb-show-compilation' and `tmb-show-r'.\n
Build and run the model using `tmb-run', `tmb-run-any', or `tmb-run-make'.
Stop the compilation or model run with `tmb-kill-process'. The C++ binary files
(*.o, *.so, *.dll) can be removed using `tmb-clean'. Invoke a GDB debug session
with `tmb-run-debug'.\n
While staying in the TMB window, navigate the secondary window with
\\<tmb-mode-map>\
\\[beginning-of-buffer-other-window], \\[scroll-other-window-down], \
\\[tmb-scroll-up] (scroll home, page up, line up), and
\\[end-of-buffer-other-window], \\[scroll-other-window], \
\\[tmb-scroll-down] (scroll end, page down, line down).
This is particularly efficient for navigating error messages listed
in the compilation buffer.\n
\\{tmb-mode-map}"
  (abbrev-mode 0)
  (modify-syntax-entry ?_ "w" tmb-mode-syntax-table)
  (set (make-local-variable 'font-lock-defaults)
       '(tmb-font-lock-keywords nil nil))
  (set (make-local-variable 'tool-bar-map) tmb-tool-bar-map)
  (setq compilation-scroll-output 'first-error)
)

(provide 'tmb)
