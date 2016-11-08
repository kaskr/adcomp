;;; tmb.el --- Major mode for creating statistical models with TMB

;; Copyright (C) 2015 Arni Magnusson

;; Author:   Arni Magnusson
;; Keywords: languages
;; URL:      http://www.hafro.is/~arnima/tmb.html

(defconst tmb-mode-version "3.1" "TMB Mode version number.")

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
;; 5. Emacs Speaks Statistics (ESS) is a required package.
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
;; If you want to set the compilation args, syntax colors, or keybindings, here
;; is an example that does that:
;;
;; (defun my-tmb-hook ()
;;   (setq tmb-compile-args ",'-fno-gnu-unique -O0 -Wall'")
;;   (setq tmb-debug-args ",'-fno-gnu-unique -g -O0'")
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
;; 10 Nov 2015  3.1  Added user function `tmb-toggle-window', internal function
;;                   `tmb-split-window', and user variable `tmb-window-right'.
;; 01 Oct 2015  3.0  Added user functions `tmb-compile' and `tmb-multi-window'.
;;                   Added user variables `tmb-compile-args' and
;;                   `tmb-debug-args'. Renamed `tmb-run-debug' to `tmb-debug',
;;                   `tmb-run-make' to `tmb-make', and `tmb-r-command' to
;;                   `tmb-compile-command'. Removed `tmb-tool-bar-map'.
;; 28 Sep 2015  2.3  Improved `tmb-toggle-nan-debug'.
;; 22 Sep 2015  2.2  Added user function `tmb-toggle-nan-debug' and internal
;;                   functions `tmb-nan-off' and `tmb-nan-on'. Added internal
;;                   variables `tmb-menu', `tmb-mode-map', and
;;                   `tmb-tool-bar-map'. Added GUI menu and toolbar. Renamed
;;                   `tmb-toggle-function' to `tmb-toggle-show-function'.
;;                   Improved `tmb-template-mini'.
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
(require 'ess-site) ; ess-*
(add-to-list 'magic-mode-alist '(tmb-include-p . tmb-mode))
(defgroup tmb nil
  "Major mode for editing Template Model Builder code."
  :tag "TMB" :group 'languages)
(defun tmb-windows-os-p () ; define here for user variables
  "Check if TMB is running in a Windows operating system."
  (if (string-match "windows" (prin1-to-string system-type)) t nil))

;; 2  User variables

(defcustom tmb-compile-args
  ;; Platform-specific: -O1 in Windows, -O0 otherwise
  (concat ",'-O" (if (tmb-windows-os-p) "1" "0") " -Wall'")
  "Arguments for compile() function in `tmb-compile'  and `tmb-template-mini'."
  :tag "Compile args" :type 'string)
(defcustom tmb-compile-command "R --quiet --vanilla"
  "Shell command to compile model using `tmb-compile'."
  :tag "Compile" :type 'string)
(defcustom tmb-debug-args
  ;; Platform-specific: -O1 and DLLFLAGS in Windows, -O0 otherwise
  (concat ",'-g -O" (if (tmb-windows-os-p) "1',DLLFLAGS=''" "0'"))
  "Arguments for compile() function in `tmb-debug'."
  :tag "Debug args" :type 'string)
(defcustom tmb-make-command "make"
  "Shell command to run makefile using `tmb-make'." :tag "Make" :type 'string)
(defcustom tmb-window-right t
  "Non-nil places secondary window on the right, nil places it below.\n
The secondary window shows compilation and model runs, among other things."
  :tag "Window right" :type 'boolean)
(defface tmb-data-face '((t :inherit font-lock-type-face))
  "Font Lock face to highlight TMB data macros." :tag "Data")
(defvar tmb-data-face 'tmb-data-face
  "Face name for TMB data macros.")
(defface tmb-parameter-face '((t :inherit font-lock-type-face))
  "Font Lock face to highlight TMB parameter macros." :tag "Parameter")
(defvar tmb-parameter-face 'tmb-parameter-face
  "Face name for TMB parameter macros.")
(defface tmb-report-face '((t :inherit font-lock-keyword-face))
  "Font Lock face to highlight TMB report macros." :tag "Report")
(defvar tmb-report-face 'tmb-report-face
  "Face name for TMB report macros.")

;; 3  Internal variables

(defvar tmb-font-lock-keywords
  (eval-when-compile
    (let ((TYPE '("Integer" "Type"))
          (DATA
           '("DATA_INTEGER" "DATA_IVECTOR" "DATA_IARRAY"
             "DATA_SCALAR"  "DATA_VECTOR"  "DATA_ARRAY"
             "DATA_FACTOR"  "DATA_STRING"  "DATA_STRUCT"
             "DATA_MATRIX"  "DATA_SPARSE_MATRIX"
             "DATA_VECTOR_INDICATOR" "DATA_ARRAY_INDICATOR"))
          (BLOCK
           '("SIMULATE"
             "PARALLEL_REGION"))
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
             "dweibull" "pweibull" "qweibull"
             ;; Debug
             "feenableexcept"))
          (CONSTANTS
           '("dim" "FE_DIVBYZERO" "FE_INVALID" "FE_OVERFLOW" "FE_UNDERFLOW"))
          (WARNINGS '("error")))
      (list
       (cons (regexp-opt TYPE 'words) font-lock-type-face)
       (cons (regexp-opt DATA 'words) 'tmb-data-face)
       (cons (regexp-opt BLOCK 'words) 'tmb-report-face)
       (cons (regexp-opt PARAMETERS 'words) 'tmb-parameter-face)
       (cons (regexp-opt REPORT 'words) 'tmb-report-face)
       (cons (regexp-opt FUNCTIONS 'words) font-lock-keyword-face)
       (cons (regexp-opt CONSTANTS 'words) font-lock-constant-face)
       (cons (regexp-opt WARNINGS 'words) font-lock-warning-face)))))
(nconc tmb-font-lock-keywords c++-font-lock-keywords)
(defvar tmb-menu
  '("TMB"
    ["View Script"         tmb-open            ]
    ["View Compilation"    tmb-show-compilation]
    ["View R Session"      tmb-show-r          ]
    "--"
    ["Compile"             tmb-compile         ]
    ["Run"                 tmb-run             ]
    ["Make"                tmb-make            ]
    "--"
    ["Stop"                tmb-kill-process    ]
    ["Clean"               tmb-clean           ]
    "--"
    ["Debug"               tmb-debug           ]
    ["Toggle NaN Debug"    tmb-toggle-nan-debug]
    "--"
    ["Mini Template"       tmb-template-mini   ]
    ["Multi-Window Layout" tmb-multi-window    ]
    ["Toggle Window"       tmb-toggle-window   ]
    "--"
    ["Help"                tmb-help            ]
    ["TMB Mode Version"    tmb-mode-version    ]))
(defvar tmb-mode-map
  ;; Don't use C-c C-                        x
  ;; Special   C-c C-       gh
  ;; Custom    C-c C- abcd f    klmnopqrst vw
  ;; Available C-c C-     e   ij          u   yz
  (let ((map (make-sparse-keymap)))
    (easy-menu-define nil map nil tmb-menu)
    (define-key map [f12]               'tmb-template-mini       )
    (define-key map [?\C-c C-backspace] 'tmb-clean               )
    (define-key map [M-up]              'tmb-scroll-up           )
    (define-key map [M-down]            'tmb-scroll-down         )
    (define-key map [?\C-c ?\C-.]       'tmb-mode-version        )
    (define-key map [?\C-c ?\C-/]       'tmb-help                )
    (define-key map [?\C-c ?\C-a]       'tmb-run-any             )
    (define-key map [?\C-c ?\C-b]       'tmb-run                 )
    (define-key map [?\C-c ?\C-c]       'tmb-compile             )
    (define-key map [?\C-c ?\C-d]       'tmb-debug               )
    (define-key map [?\C-c ?\C-f]       'tmb-for                 )
    (define-key map [?\C-c ?\C-k]       'tmb-kill-process        )
    (define-key map [?\C-c ?\C-l]       'tmb-show-compilation    )
    (define-key map [?\C-c ?\C-m]       'tmb-make                )
    (define-key map [?\C-c ?\C-n]       'tmb-toggle-nan-debug    )
    (define-key map [?\C-c ?\C-o]       'tmb-open-any            )
    (define-key map [?\C-c ?\C-p]       'tmb-open                )
    (define-key map [?\C-c ?\C-q]       'tmb-kill-process        )
    (define-key map [?\C-c ?\C-r]       'tmb-show-r              )
    (define-key map [?\C-c ?\C-s]       'tmb-toggle-show-function)
    (define-key map [?\C-c ?\C-t]       'tmb-toggle-window       )
    (define-key map [?\C-c ?\C-v]       'tmb-run                 )
    (define-key map [?\C-c ?\C-w]       'tmb-multi-window        )
    (define-key map [?\C-\M-v]          'ignore                  )
    map))

;; 4  User functions

(defun tmb-clean ()
  "Remove C++ binary files (*.o *.so *.dll)." (interactive)
  (let* ((prog (file-name-sans-extension (buffer-name)))
         (pattern (concat prog "\\.o\\|" prog "\\.so\\|" prog "\\.dll"))
         (files (directory-files "." nil pattern t)))
    (dolist (x files)(delete-file x)))(message "Removed binary files"))
(defun tmb-compile ()
  "Compile model, using `tmb-compile-command' and `tmb-compile-args'."
  (interactive)(save-buffer)(tmb-split-window)
  (compile (concat tmb-compile-command " -e \"require(TMB); compile('"
                   (buffer-name) "'" tmb-compile-args ")\""))
  (with-current-buffer "*compilation*" (setq show-trailing-whitespace nil)))
(defun tmb-debug ()
  "Debug model with GDB, using `tmb-debug-args'.\n
The R session stays alive if it was running when this function was called."
  (interactive)(save-buffer)(tmb-split-window)
  (let* ((ess-dialect "R")
         (inferior-R-args "--quiet --vanilla")
         (ess-ask-for-ess-directory nil)
         (prefix (file-name-sans-extension (buffer-name)))
         (cmd-1 "require(TMB)")
         (cmd-2 (concat "; compile('" prefix ".cpp'" tmb-debug-args ")"))
         (cmd-3 (concat "; gdbsource('" prefix ".R',TRUE)"))
         (cmd-4 (if (ess-process-live-p) "" "; q()"))
         (cmd (concat cmd-1 cmd-2 cmd-3 cmd-4)))
    (ess-eval-linewise cmd)))
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
(defun tmb-multi-window ()
  "Arrange windows with C++ left, R right, and *R* below." (interactive)
  (let ((ess-dialect "R")
        (inferior-R-args "--quiet --vanilla")
        (ess-ask-for-ess-directory nil)
        (r-script (concat (file-name-sans-extension (buffer-name)) ".R")))
    (tmb-open)(ess-eval-linewise "")(delete-other-windows)
    (split-window-vertically -16)(set-window-buffer (next-window) "*R*")
    (split-window-horizontally)(set-window-buffer (next-window) r-script)))
(defun tmb-make ()
  "Run makefile in current directory, using `tmb-make-command'."
  (interactive)(save-buffer)(tmb-split-window)(compile tmb-make-command)
  (with-current-buffer "*compilation*" (setq show-trailing-whitespace nil)))
(defun tmb-mode-version ()
  "Show TMB Mode version number." (interactive)
  (message "TMB Mode version %s" tmb-mode-version))
(defun tmb-open ()
  "Open R script with same filename prefix as current buffer." (interactive)
  (tmb-open-any "R"))
(defun tmb-open-any (ext)
  "Open file with extension EXT in other window." (interactive "sExtension: ")
  (let ((file (concat (file-name-sans-extension (buffer-name)) "." ext)))
    (if (not (file-regular-p file))(error "File %s not found" file)
      (tmb-split-window)(save-selected-window (find-file-other-window file)))))
(defun tmb-run ()
  "Run R script with same filename prefix as current buffer.\n
If the R script has a different filename, then use `tmb-run-any' instead.\n
The script is sourced in an existing R session, or a new R session is started."
  (interactive)
  (tmb-run-any (concat (file-name-sans-extension (buffer-name)) ".R")))
(defun tmb-run-any (script)
  "Run any R script, querying user for SCRIPT filename.\n
If the R script has the same filename prefix as the current buffer, then use
`tmb-run' instead.\n
Filename history is accessible in the minibuffer prompt \
(\\<minibuffer-local-map>\\[previous-history-element],\
 \\[next-history-element]).\n
The script is sourced in an existing R session, or a new session is started."
  (interactive "fRun R script: ")(save-buffer)
  (let* ((ess-dialect "R")
         (inferior-R-args "--quiet --vanilla")
         (ess-ask-for-ess-directory nil))
    (tmb-split-window)(ess-load-file script)))
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
      (error "*compilation* buffer not found")
    (tmb-split-window)(display-buffer "*compilation*")))
(defun tmb-show-r ()
  "Show R interactive buffer." (interactive)
  (if (null (get-buffer "*R*"))(error "*R* interactive buffer not found")
    (tmb-split-window)(ess-show-buffer "*R*")))
(defun tmb-template-mini ()
  "Create minimal TMB files (mini.cpp, mini.R) in current directory.\n
The user variable `tmb-compile-args' is passed to the compile() function."
  (interactive)(delete-other-windows)(find-file "mini.cpp")
  (delete-region (point-min)(point-max))(insert "\
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
  (goto-char (point-min))(write-file "mini.cpp" t)
  (save-selected-window
    (tmb-split-window)(find-file-other-window "mini.R")
    (delete-region (point-min)(point-max))(insert "\
data <- list(x=rivers)
parameters <- list(mu=0, logSigma=0)

require(TMB)
compile('mini.cpp'" tmb-compile-args ")
dyn.load(dynlib('mini'))

################################################################################

model <- MakeADFun(data, parameters)
fit <- nlminb(model$par, model$fn, model$gr)
rep <- sdreport(model)

print(rep)
")
    (goto-char (point-min))(write-file "mini.R" t))
  (message (concat "Ready to compile R script ("
                   (substitute-command-keys "\\<tmb-mode-map>\\[tmb-compile]")
                   ") or edit code.")))
(defun tmb-toggle-nan-debug ()
  "Toggle floating point exceptions, for debugging." (interactive)
  (save-excursion
    (goto-char (point-min))
    (if (re-search-forward "#include *<fenv.h>" nil t)
        (tmb-nan-off)(tmb-nan-on))))
(defun tmb-toggle-show-function ()
  "Toggle whether to show the current function name in the mode line."
  (interactive)(which-function-mode (if which-function-mode 0 1))
  (message "Function indicator %s" (if which-function-mode "ON" "OFF")))
(defun tmb-toggle-window ()
  "Toggle whether secondary window is on the right or below." (interactive)
  (delete-other-windows)(setq tmb-window-right (not tmb-window-right))
  (message "Secondary window %s" (if tmb-window-right "RIGHT" "BELOW")))

;; 5  Internal functions

(defun tmb-include-p ()
  "Check if C++ file has #include <TMB.hpp>."
  (if (string-equal (file-name-extension (buffer-name)) "cpp")
      (re-search-forward "#include *<TMB.hpp>"
                         magic-mode-regexp-match-limit t)))
(defun tmb-nan-off ()
  "Disable floating point exceptions."
  (flush-lines "#include *<fenv.h>" (point-min)(point-max))
  (flush-lines "feenableexcept\(.*\);" (point-min)(point-max))
  (message "Floating point exceptions disabled"))
(defun tmb-nan-on ()
  "Enable floating point exceptions."
  (goto-char (point-min))(re-search-forward "#include *<TMB.hpp>")
  (insert "\n#include <fenv.h>")
  (search-forward "objective_function<Type>::operator()")(search-forward "{")
  (insert "\n  feenableexcept"
          "(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO /* | FE_UNDERFLOW */ );")
  (message "Floating point exceptions enabled"))
(defun tmb-split-window ()
  "Split window if it is the only window, otherwise do nothing.\n
The orientation of the split depends on the value of `tmb-window-right'."
  (if (one-window-p)(if tmb-window-right
                        (split-window-horizontally)(split-window-vertically))))

;; 6  Main function

;;;###autoload
(define-derived-mode tmb-mode c++-mode "TMB"
  "Major mode for creating statistical models with TMB.\n
The `tmb-help' command shows this page.\n
Start a new model from `tmb-template-mini'. Navigate between functions using
`imenu' and show the current function name in the mode line with
`tmb-toggle-show-function'. Prototype for-loops with `tmb-for'.\n
Use `tmb-open' to open the corresponding R script and `tmb-open-any' to open
other model-related files. Show interactive *compilation* and *R* buffers with
`tmb-show-compilation' and `tmb-show-r'. Use `tmb-toggle-window' to set
`tmb-window-right' to your viewing preference. The `tmb-multi-window' layout is
an alternative to the default left-right layout.\n
Build and run the model using `tmb-compile', `tmb-run', and `tmb-make'.
Stop the compilation or model run with `tmb-kill-process'.\n
The C++ binary files (*.o, *.so, *.dll) can be removed using `tmb-clean'.
Invoke a GDB debug session with `tmb-debug', where `tmb-toggle-nan-debug'
can be helpful.\n
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
  (setq compilation-scroll-output 'first-error))

(provide 'tmb)
