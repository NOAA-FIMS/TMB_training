// Simple linear regression
#include <TMB.hpp>
template <class Type>
  
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta);
  PARAMETER(lnSigma);
  
  Type nll = 0;
  Type sigma = exp(lnSigma);
  int n = y.size();
  vector<Type> mu = X * beta;
  for(int i=0; i<n; i++){
    nll -= dnorm(y(i), mu(i), sigma, true);
  }
  SIMULATE{
    y = rnorm(mu, sigma);
    REPORT(y);
  }
  
  Type sig2 = pow(sigma,2);
  REPORT(sig2);
  ADREPORT(sig2);
  
  return nll;
}
