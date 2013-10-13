README
======

This folder contain the ADMB examples.

- Folder names should be the same as RcppAD example names.
- Each ADMB example should return sdreport output which will be automatically compared with the corresponding RcppAD example.
- For the same reason, each RcppAD example should include "sdreport(obj)".
- The order of random effects in RcppAD and ADMB should be the same.
- Each ADMB example should be timed with and without sdreport so that the sdreport time can be identified.

