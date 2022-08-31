// Spatial poisson GLMM 
#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(y);
  DATA_MATRIX(D); //distance matrix
  
  PARAMETER(beta0);
  PARAMETER(ln_phi);
  PARAMETER(logit_power);
  PARAMETER(ln_range); //spatial range parameter
  PARAMETER(ln_sigma2); //spatial sd
  PARAMETER_VECTOR(omega);

  Type phi = exp(ln_phi);
  Type power = invlogit(logit_power) + 1;
  Type range = exp(ln_range);
  Type sigma2 = exp(ln_sigma2);
  
  //Spatial Random Effects
  //Covariance Matrix
  matrix<Type> C(D);
  for(int i=0; i<C.rows(); i++){
    for(int j=0; j<C.cols(); j++){
      C(i,j) = sigma2 * matern(D(i,j), range, Type(1.0));
    }
  }
    
  Type nll_omega = MVNORM(C)(omega);
  SIMULATE{
    MVNORM(C).simulate(omega);
    REPORT(omega);
  }
  
  vector<Type> mu = exp(beta0 + omega);
  Type nll_y = -sum(dtweedie(y, mu, phi, power, true));
  SIMULATE{
    y = rtweedie(mu, phi, power);
    REPORT(y);
  }
  
  Type Total_Weight = mu.sum();
  
  Type nll = nll_omega + nll_y;
  
  REPORT(phi);
  REPORT(power);
  REPORT(range);
  REPORT(sigma2);
  REPORT(nll_omega);
  REPORT(nll_y);
  REPORT(nll);
  REPORT(C);
 
  ADREPORT(Total_Weight);
  
  return nll;
}
