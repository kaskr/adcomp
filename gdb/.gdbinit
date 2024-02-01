## These helper functions are useful when debugging TMB internal functions

define printconfig
  set pagination off
  set print elements 0
  set max-value-size unlimited
  set print repeats unlimited
end

define tidyoutput
  shell sed -i 's/{/(/g' gdb.txt
  shell sed -i 's/}/)/g' gdb.txt
end
document tidyoutput
  Tidy the gdb.txt file so it can be read by R using 'dget'. Usage:
    tidyoutput
end

define getvec
  printconfig
  set $rows=$arg0.m_storage.m_rows
  echo =c
  eval "output *(double (*)[%i])$arg0.m_storage.m_data", $rows
  echo \n
end
document getvec
  Write R representation of Eigen vector to file 'gdb.txt'. Usage:
    set logging on
    getvec(M)
    set logging off
end

define getmat
  printconfig
  set $rows=$arg0.m_storage.m_rows
  set $cols=$arg0.m_storage.m_cols
  echo =matrix(c
  eval "output *(double (*)[%i])$arg0.m_storage.m_data", $rows*$cols
  printf ", %i)", $rows
  echo \n
end
document getmat
  Write R representation of Eigen vector to file 'gdb.txt'. Usage:
    set logging on
    getmat(M)
    set logging off
end

define getspmat
  printconfig
  set $rows=$arg0.m_innerSize
  set $cols=$arg0.m_outerSize
  set $size=$arg0.m_data.m_size
  echo require(Matrix)\n
  echo new("dgCMatrix",
  echo     p=as.integer(c(c
  eval     "output *(int (*)[%i])$arg0.m_outerIndex", $cols
  echo     ,
  output   $size
  echo     ))
  echo ,
  echo     i=as.integer(c
  eval     "output *(int (*)[%i])$arg0.m_data.m_indices", $size
  echo     )
  echo ,
  echo     x=as.double(c
  eval     "output *(double (*)[%i])$arg0.m_data.m_values", $size
  echo     )
  echo ,
  echo     Dim=as.integer(c(
  output       $rows
  echo ,
  output       $cols
  echo     ))
  echo    )  
  echo \n
end
document getspmat
  Write R representation of Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> > to file 'gdb.txt'.
  Usage:
    set logging on
    getspmat(M)
    set logging off
    tidyoutput
end
