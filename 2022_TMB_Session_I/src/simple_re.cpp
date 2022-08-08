// Normal linear mixed model specified through sparse design matrices.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);         
  DATA_SPARSE_MATRIX(B);  
  DATA_SPARSE_MATRIX(X); 
  DATA_VECTOR_INDICATOR(keep, y); 
  PARAMETER_VECTOR(u);    
  PARAMETER_VECTOR(beta); 
  PARAMETER(lnSDu);      
  PARAMETER(lnSDy);   

  Type nll1 = 0;
  Type nll2 = 0;
  Type nll = 0;
  Type sdu = exp(lnSDu);
  nll1 -= sum(dnorm(u, Type(0), sdu, true));

  vector<Type> mu = X * beta + B * u;
  Type sdy = exp(lnSDy);
  for(int i=0; i<y.size(); i++){
    nll2 -= keep(i) * dnorm(y(i), mu(i), sdy, true);
  }

  nll = nll1 + nll2;
  REPORT(nll1);
  REPORT(nll2);
  REPORT(nll);

   return nll;
}

