// Normal linear mixed model specified through sparse design matrices.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);         
  DATA_SPARSE_MATRIX(B);  
  DATA_SPARSE_MATRIX(X); 
  PARAMETER_VECTOR(u);    
  PARAMETER_VECTOR(beta); 
  PARAMETER(lnSDu);      
  PARAMETER(lnSDy);     

  Type nll = 0;
  Type sdu = exp(lnSDu);
  nll -= dnorm(u, Type(0), sdu, true).sum();

  vector<Type> mu = X * beta + B * u;
  Type sdy = exp(lnSDy);
  nll -= dnorm(y, mu, sdy, true).sum();

   return nll;
}

