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
Type rbetabinom(Type n, Type mu, Type phi)
{
  Type x = rbeta(mu*phi,(Type(1.0)-mu)*phi);
  Type y = rbinom(n, x);
  return(y);  
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y); //number mature
  DATA_VECTOR(N); //number of mature + not mature
  DATA_VECTOR(age_obs); //age for each observation
  DATA_IVECTOR(re_ind); //indicator for random effect for each observation
  DATA_INTEGER(max_age); //to generate a maturity ogive
  PARAMETER_VECTOR(beta); //intercept, slope
  PARAMETER(log_phi);
  PARAMETER_VECTOR(AR_pars); //AR1 process (2 pars)
  PARAMETER_VECTOR(re); //annual ar1 process
  Type phi = exp(log_phi);

  int n_obs = Y.size();
  vector<Type> nll_re(re.size());
  vector<Type> nll(n_obs);
  nll.setZero();
  nll_re.setZero();

  Type sig_re = exp(AR_pars(0));
  Type rho_re = -1 + 2/(1+exp(-AR_pars(1)));
  nll_re(0) -= dnorm(re(0), Type(0), sig_re*exp(-Type(0.5) * log(1 - pow(rho_re,Type(2)))), 1);
  for(int y = 1; y < re.size(); y++) nll_re(y) -= dnorm(re(y), rho_re * re(y-1), sig_re, 1);
  SIMULATE{
    re(0) = rnorm(Type(0), sig_re*exp(-Type(0.5) * log(Type(1) - pow(rho_re,Type(2)))));
    for(int y = 1; y < re.size(); y++) re(y) = rnorm(rho_re * re(y-1), sig_re);
    REPORT(re);
  }

  vector<Type> logit_mat(n_obs), mat(n_obs);
  for(int i = 0; i < n_obs; i++) {
    logit_mat(i) = beta(0)  + beta(1)*age_obs(i) + re(re_ind(i)); 
    mat(i) = 1/(1 + exp(-logit_mat(i)));
    nll(i) = -dbetabinom(Y(i), N(i), mat(i), phi, 1); //negative log-likelihood
    //nll(i) = -dbinom(Y(i), N(i), mat(i), 1);   //negative log-likelihood
  }
  SIMULATE {
    for(int i = 0; i < n_obs; i++) Y(i) = rbetabinom(N(i), mat(i), phi);
    REPORT(Y);
  }

  vector<Type> logit_mat_at_age(max_age);//, mat_at_age(max_age);
  for(int i = 0; i < max_age; i++) logit_mat_at_age(i) = beta(0)  + beta(1)*Type(i+1);
  vector<Type> mat_at_age = 1/(1+ exp(-logit_mat_at_age));
  vector<Type> a50 = -(beta(0) + re)/beta(1);
  REPORT(mat);
  ADREPORT(phi);
  ADREPORT(a50);
  REPORT(nll);
  REPORT(nll_re);
  ADREPORT(logit_mat_at_age);
  ADREPORT(mat_at_age);
  return sum(nll) + sum(nll_re);
}

