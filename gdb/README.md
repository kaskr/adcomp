## Some useful gdb functions for TMB debugging

### Install

Copy file to `~/.gdbinit`.

### Example

- `cd tmb_examples` and run spde example `make spde`. We want to check that matrix 'Q' is as expected.
- Start gdb by `R -d gdb` and set a breakpoint `b spde.cpp:55`.
- Now `run` will start R, and we can `source("spde.R")`.

The debugger stops with something like

```
Breakpoint 1, objective_function<double>::operator() (
    this=this@entry=0x7fffffff30c0) at spde.cpp:57
57	}
```

We want to check the content of Q, but printing it is not so useful:

```
(gdb) print Q
$3 = {<Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >> = {<Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> >> = {<Eigen::EigenBase<Eigen::SparseMatrix<double, 0, int> >> = {<No data fields>}, m_isRValue = false}, <No data fields>}, 
  m_outerSize = 1721, m_innerSize = 1721, m_outerIndex = 0x55555f4b97d0, 
  m_innerNonZeros = 0x0, m_data = {m_values = 0x55555cdeda30, m_indices = 0x55555e6813f0, 
    m_size = 33759, m_allocatedSize = 55102}}
```

To get it into R we can instead do:

```
set logging on
getspmat(Q)
set logging off
tidyoutput
```

This generates a file `gdb.txt` that can be sourced from R.
