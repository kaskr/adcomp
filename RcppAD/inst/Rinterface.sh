echo "library(RcppAD)"
echo "dyn.load(\"$1\")" | sed s/\.cpp/\.so/g
echo "MakeADFun("
echo " data=list("
grep "DATA[_(].*)" $1 | sed s/.*\(/\ \ /g | sed s/\).*/=,/g
echo "  )"
echo " ),"
echo " parameters=list("
grep "PARAMETER[_(].*)" $1 | sed s/.*\(/\ \ /g | sed s/\).*/=,/g
echo " )"
echo ")"
