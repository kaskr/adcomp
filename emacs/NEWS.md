# TMB Mode 3.6 (2026-10-09)

* Maintenance release, adapting to recent changes in Emacs and ESS.




# TMB Mode 3.51 (2019-04-02)

* Added keywords "DATA_IMATRIX" and "vec".




# TMB Mode 3.5 (2018-01-28)

* Changed `tmb-run` and `tmb-run-any` so they source(echo=TRUE).




# TMB Mode 3.4 (2018-01-14)

* Added GUI toolbar for common tasks.

* Bound f11 to `tmb-open`.




# TMB Mode 3.3 (2017-10-14)

* Adjusted compiler flags to fix the combination precompile + debugging.

* Adapted mini example to follow style of `tmb_examples` and `setupRStudio`.




# TMB Mode 3.2 (2016-11-16)

* Added user variable `tmb-block-face`.

* Added keywords "DATA_STRING", "SIMULATE", "PARALLEL_REGION", "rnorm", "rpois",
  "rnbinom", "rnbinom2", "rgamma", and "simulate".

* Improved `tmb-template-mini` so it asks for model name.




# TMB Mode 3.1 (2015-11-10)

* Added user function `tmb-toggle-window`.

* Added internal function `tmb-split-window`.

* Added user variable `tmb-window-right`.




# TMB Mode 3.0 (2015-10-01)

* Added user functions `tmb-compile` and `tmb-multi-window`.

* Added user variables `tmb-compile-args` and `tmb-debug-args`.

* Renamed `tmb-run-debug` to `tmb-debug`, `tmb-run-make` to `tmb-make`, and
  `tmb-r-command` to `tmb-compile-command`.

* Removed `tmb-tool-bar-map`.




# TMB Mode 2.3 (2015-09-28)

* Improved `tmb-toggle-nan-debug`.




# TMB Mode 2.2 (2015-09-22)

* Added user function `tmb-toggle-nan-debug` and internal functions
  `tmb-nan-off` and `tmb-nan-on`.

* Added internal variables `tmb-menu`, `tmb-mode-map`, and `tmb-tool-bar-map`.

* Added GUI menu and toolbar.

* Renamed `tmb-toggle-function` to `tmb-toggle-show-function`.

* Improved `tmb-template-mini`.




# TMB Mode 2.1 (2015-09-10)

* Added internal function `tmb-windows-os-p`.

* Improved `tmb-run-debug` and `tmb-template-mini`.




# TMB Mode 2.0 (2015-09-07)

* Added user functions `tmb-run-debug`, `tmb-scroll-down`, `tmb-scroll-up`,
  `tmb-show-compilation`, and `tmb-show-r`.

* Renamed `tmb-open` to `tmb-open-any` and `tmb-run-r` to `tmb-run`.

* Improved `tmb-open`, `tmb-open-any`, `tmb-run`, `tmb-run-any`, and
  `tmb-template-mini`.




# TMB Mode 1.3 (2015-09-05)

* Added user function `tmb-template-mini`.

* Renamed `tmb-toggle-section` to `tmb-toggle-function`.




# TMB Mode 1.2 (2015-09-04)

* Added user functions `tmb-clean`, `tmb-for`, `tmb-kill-process`, `tmb-open`,
  `tmb-open-r`, `tmb-run-any`, `tmb-run-make`, `tmb-run-r`, and
  `tmb-toggle-section`.

* Added user variables `tmb-make-command` and `tmb-r-command`.

* Disabled `abbrev-mode`.




# TMB Mode 1.1 (2015-09-03)

* Shortened list of recognized FUNCTIONS for maintainability.




# TMB Mode 1.0 (2015-09-01)

* Created main function `tmb-mode`, derived from `c++-mode`.
