/* One dimensional normal density */
template<class Type>
Type dnorm(Type x, Type mean, Type sd, int give_log=0)
{
  Type logres;
  logres=-log(Type(sqrt(2*M_PI))*sd)-Type(.5)*pow((x-mean)/sd,2);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}
