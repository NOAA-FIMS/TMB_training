// Simple linear regression
#include <TMB.hpp>

template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(lnSigma);
  
  Type nll = 0;
  vector<Type> mu = X*beta;
  vector<Type> sigma = exp(lnSigma);

  for(int i=0; i<y.size(); i++){
    nll -= dnorm(y(i), mu(i), sigma(i));
  }
  
  REPORT(sigma);
  ADREPORT(sigma);
  
  return nll;
}
