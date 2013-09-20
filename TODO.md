TODO
====
- Clean up R-code: newton, options, MakeADFun: class(output).
- Roxygen documentation
- MakeADFun: Call dyn.load if DLL not loaded (?)
- parList: use last.par.best as default (?)
- parallel_start.hpp does not yet have parallel Hessian member.
- Metis symbolic analysis with 3d example.
- Move methods, such as "parList", out of the object.
- Rinterface should remember to set DLL="..."
- sdreport() on R-side to get sd of ADREPORT().
- solveSubset to get sd's of all random effects.
- Make importence sampler work in high dim - need a GMRFsample.
- Compile and testing workflow:
  - Eliminate need to restart R.
  - Give better message if PARAMETER(name) evaluates to NULL.
  - Improve compile times by either precompiling or allow incomplete compilations.
  - Array bounds checking option