#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";

template<class Type>
Type dbetabinom(Type x, Type n, Type p, Type phi, int do_log)
{
  Type ll = lgamma(n + 1.0) - lgamma(x + 1.0) - lgamma(n - x + 1.0) + 
    lgamma(x + p*phi) + lgamma(n - x +(1-p)*phi) - lgamma(n + phi) +
    lgamma(phi) - lgamma(p*phi) - lgamma((1-p)*phi);
  if(do_log) return(ll);
  else return(exp(ll));  
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y); //number mature
  DATA_VECTOR(N); //number of mature + not mature
  DATA_MATRIX(X); //includes age
  DATA_INTEGER(max_age); //to generate a maturity ogive
  DATA_INTEGER(age_col); //which column of X has age
  PARAMETER_VECTOR(beta); //length = NCOL(X)
  PARAMETER(log_phi);
  Type phi = exp(log_phi);

  int n_obs = Y.size();
  vector<Type> nll(n_obs);
  nll.setZero();
     
  vector<Type> logit_mat(n_obs), mat(n_obs);
  logit_mat = X * beta;
  mat = 1/(1 + exp(-logit_mat));
  for(int i = 0; i < n_obs; i++) {
    nll(i) = -dbetabinom(Y(i), N(i), mat(i), phi, 1); //negative log-likelihood
  }

  //vector<Type> logit_mat_at_age(max_age);//, mat_at_age(max_age);
  //for(int i = 0; i < max_age; i++) logit_mat_at_age(i) = beta(0)  + beta(1)*Type(i+1);
  //vector<Type> mat_at_age = 1/(1+ exp(-logit_mat_at_age));
  vector<Type> a50(X.cols()-1);
  int index = 0;
  for(int i = 0; i < a50.size(); i++){
    if(i != age_col-1){
      a50(index) = -beta(i)/beta(age_col-1);
      index++;
    }
  }
  //Type a50 = -beta(0)/beta(1);
  REPORT(mat);
  ADREPORT(phi);
  ADREPORT(a50);
  REPORT(nll);
  //ADREPORT(logit_mat_at_age);
  //ADREPORT(mat_at_age);
  return sum(nll);
}

