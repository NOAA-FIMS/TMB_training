#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";
template<class Type>
Type dbetabinom(Type x, Type n, Type mu, Type phi, int do_log)
{
  Type ll = lgamma(n + 1.0) - lgamma(x + 1.0) - lgamma(n - x + 1.0) + 
	  lgamma(x + mu*phi) + lgamma(n - x +(1-mu)*phi) - lgamma(n + phi) +
	  lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi);
  if(do_log == 1) return(ll);
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
  DATA_IVECTOR(age_obs); //age for each observation
  DATA_IVECTOR(year_obs); //age for each observation
  DATA_IVECTOR(cy_obs); //cohort year for each observation (year of birth)
  DATA_IVECTOR(Ecov_maxages_k); //last year of Ecov to use (nobs)
  int max_ecov_ages = 0;
  for(int i = 0; i < Ecov_maxages_k.size(); i++) if(Ecov_maxages_k(i) > max_ecov_ages) max_ecov_ages = Ecov_maxages_k(i);
  DATA_IVECTOR(Ecov_maxages_a50); //last year of Ecov to use (nobs)
  for(int i = 0; i < Ecov_maxages_a50.size(); i++) if(Ecov_maxages_a50(i) > max_ecov_ages) max_ecov_ages = Ecov_maxages_a50(i);
  DATA_VECTOR(Ecov_obs);
  DATA_IVECTOR(use_Ecov_obs);
  DATA_VECTOR(Ecov_obs_sigma);
  DATA_INTEGER(fit_k);
  DATA_INTEGER(fit_a50);
  DATA_INTEGER(binomial);
  
  PARAMETER(beta_k);
  PARAMETER(beta_a50);
  PARAMETER_VECTOR(beta_Ecov_k);
  PARAMETER_VECTOR(beta_Ecov_a50);
  PARAMETER(Ecov_mu);
  PARAMETER_VECTOR(Ecov_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(Ecov_re);
  PARAMETER(beta_phi);
  PARAMETER_VECTOR(k_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(k_re);
  PARAMETER_VECTOR(a50_AR_pars); //AR1 process (2)
  PARAMETER_VECTOR(a50_re);

  int n_obs = Y.size();
  int n_years_Ecov = Ecov_obs.size();// + max_ecov_ages;
  //see(n_obs);
  
  Type zero = Type(0);
  Type one = Type(1);
  Type nll= zero; //negative log-likelihood
  
  vector<Type> Ecov_y(n_years_Ecov);
  Type Ecov_phi;
  Type Ecov_sig;
  Type nll_Ecov = zero;
  Ecov_phi = -one + Type(2)/(one + exp(-Ecov_AR_pars(0)));
  Ecov_sig = exp(Ecov_AR_pars(1));
  for(int y = 0; y < n_years_Ecov; y++) Ecov_y(y) = Ecov_mu + Ecov_re(y);
  nll_Ecov -= dnorm(Ecov_re(0), zero, Ecov_sig*exp(-Type(0.5) * log(one - pow(Ecov_phi,Type(2)))), 1);
  SIMULATE
  {
    Ecov_re(0) = rnorm(zero, Ecov_sig*exp(-Type(0.5) * log(one - pow(Ecov_phi,Type(2)))));
    Ecov_y(0) = Ecov_mu + Ecov_re(0);
  }
  for(int y = 1; y < n_years_Ecov; y++) nll -= dnorm(Ecov_re(y), Ecov_phi * Ecov_re(y-1), Ecov_sig, 1);
  SIMULATE
  {
    for(int y = 1; y < n_years_Ecov; y++) 
    {
      Ecov_re(y) = rnorm(Ecov_phi * Ecov_re(y-1), Ecov_sig);
      Ecov_y(y) = Ecov_mu + Ecov_re(y);
    }
    REPORT(Ecov_re);
  }
  nll += nll_Ecov;
  //see(nll_Ecov);
  Type nll_Ecov_obs = zero;
  for(int y = 0; y < Ecov_obs.size(); y++) if(use_Ecov_obs(y) == 1) nll_Ecov_obs -= dnorm(Ecov_obs(y), Ecov_y(y), Ecov_obs_sigma(y), 1);
  SIMULATE 
  {
    for(int y = 0; y < Ecov_obs.size(); y++) if(use_Ecov_obs(y) == 1) Ecov_obs(y) = rnorm(Ecov_y(y), Ecov_obs_sigma(y));  
    REPORT(Ecov_obs);
  }
  nll += nll_Ecov_obs;
  //see(nll_Ecov_obs);
  
  //k AR(1) process
  Type k_phi;
  Type k_sig;
  Type nll_k = zero;
  k_phi = -one + Type(2)/(one + exp(-k_AR_pars(0)));
  k_sig = exp(k_AR_pars(1));
  nll_k -= dnorm(k_re(0), zero, k_sig*exp(-Type(0.5) * log(one - pow(k_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll_k -= dnorm(k_re(y), k_phi * k_re(y-1), k_sig, 1);
  if(fit_k == 1) nll += nll_k;
  //see(nll_k);
  SIMULATE
  {
    k_re(0) = rnorm(zero, k_sig*exp(-Type(0.5) * log(one - pow(k_phi,Type(2)))));
    for(int y = 1; y < n_years_Ecov; y++) 
    {
      k_re(y) = rnorm(k_phi * k_re(y-1), k_sig);
    }
    REPORT(k_re);
  }
  //a50 AR(1) process
  Type a50_phi;
  Type a50_sig;
  Type nll_a50 = zero;
  a50_phi = -one + Type(2)/(one + exp(-a50_AR_pars(0)));
  a50_sig = exp(a50_AR_pars(1));
  nll_a50 -= dnorm(a50_re(0), zero, a50_sig*exp(-Type(0.5) * log(one - pow(a50_phi,Type(2)))), 1);
  for(int y = 1; y < n_years_Ecov; y++) nll_a50 -= dnorm(a50_re(y), a50_phi * a50_re(y-1), a50_sig, 1);
  if(fit_a50 == 1) nll += nll_a50;
  //see(nll_a50);
  SIMULATE
  {
    a50_re(0) = rnorm(zero, a50_sig*exp(-Type(0.5) * log(one - pow(a50_phi,Type(2)))));
    for(int y = 1; y < n_years_Ecov; y++) 
    {
      a50_re(y) = rnorm(a50_phi * a50_re(y-1), a50_sig);
    }
    REPORT(a50_re);
  }
  
  vector<Type> l_k(n_obs);
  vector<Type> l_a50(n_obs);
  vector<Type> logit_mat(n_obs);
  for(int i = 0; i < n_obs; i++) 
  {
    l_k(i) = beta_k; 
    for(int y = 0; y <= Ecov_maxages_k(i); y++) l_k(i) += Ecov_y(year_obs(i) - 1 - age_obs(i) + y) * beta_Ecov_k(y) + k_re(year_obs(i) -1 -age_obs(i) + y);
    l_a50(i) = beta_a50;
    for(int y = 0; y <= Ecov_maxages_a50(i); y++) l_a50(i) += Ecov_y(year_obs(i) - 1 - age_obs(i) + y) * beta_Ecov_a50(y) + a50_re(year_obs(i) -1 -age_obs(i) + y);
    logit_mat(i) = exp(l_k(i)) * (Type(age_obs(i)) - exp(l_a50(i)));
  }
  vector<Type> mat = one/(one + exp(-logit_mat));
  Type phi = exp(beta_phi);
  vector<Type> ll_mat(n_obs); 
  ll_mat.setZero();
  if(binomial == 1) 
  {
    ll_mat = dbinom(Y, N, mat, 1); //binomial?
    SIMULATE 
    {
      Y = rbinom(N, mat);
      REPORT(Y);
    }
  }
  else 
  {
    for(int i = 0; i < n_obs; i++) ll_mat(i) = dbetabinom(Y(i), N(i), mat(i), phi, 1);
    SIMULATE 
    {
      for(int i = 0; i < n_obs; i++) Y(i) = rbetabinom(N(i), mat(i), phi);
      REPORT(Y);
    }
  }
  //see(sum(ll_mat));
  nll -= sum(ll_mat);
  
  matrix<Type> log_a50(n_years_Ecov,beta_Ecov_a50.size()), log_k(n_years_Ecov,beta_Ecov_k.size());
  log_a50.setZero();
  for(int i = 0; i < beta_Ecov_a50.size(); i++) //age_obs
  {
    for(int y = i+1; y < n_years_Ecov; y++) //year_obs
    {
      log_a50(y,i) = beta_a50;
      for(int j = 0; j <= i; j++) log_a50(y,i) += beta_Ecov_a50(j) * Ecov_y(y - 1 - i + j);
    }
  }
  log_k.setZero();
  for(int i = 0; i < beta_Ecov_k.size(); i++) //age_obs
  {
    for(int y = i+1; y < n_years_Ecov; y++) //year_obs
    {
      log_k(y,i) = beta_k;
      for(int j = 0; j <= i; j++) log_k(y,i) += beta_Ecov_k(j) * Ecov_y(y - 1 - i + j);
    }
  }

  matrix<Type> logit_pmat(n_years_Ecov,3), pmat(n_years_Ecov,3);
  for(int y = 1; y < n_years_Ecov; y++) for(int a = 0; a < 3; a++) if(y-1>=a)
  {
    Type temp = beta_k;
    int maxage = a+1;
    if(a+1 > beta_Ecov_k.size()) maxage = beta_Ecov_k.size();
    //see(y);
    //see(a);
    for(int i = 0; i < maxage; i++) temp += beta_Ecov_k(i) * Ecov_y(y-1 - a + i) + k_re(y -1 - a + i);
    logit_pmat(y,a) = exp(temp);
    temp = beta_a50;
    maxage = a+1;
    if(a+1 > beta_Ecov_a50.size()) maxage = beta_Ecov_a50.size();
    for(int i = 0; i < maxage; i++) temp += beta_Ecov_a50(i) * Ecov_y(y-1 - a + i) + a50_re(y -1 - a + i);
    logit_pmat(y,a) *= Type(a+1) - exp(temp);
    pmat(y,a) = one/(one + exp(-logit_pmat(y,a)));
  }
  matrix<Type> log_a50_Ecov(n_years_Ecov,beta_Ecov_a50.size()), log_k_Ecov(n_years_Ecov,beta_Ecov_k.size()), logit_p_Ecov(n_years_Ecov,3);
  log_k_Ecov.fill(beta_k);
  log_a50_Ecov.fill(beta_a50);
  logit_p_Ecov.setZero();
  for(int y = 0; y < n_years_Ecov; y++)
  {
    for(int i = 0; i < beta_Ecov_k.size(); i++) for(int j = 0; j <= i; j++) log_k_Ecov(y,i) += beta_Ecov_k(j) * Ecov_y(y) + k_re(y);
    for(int i = 0; i < beta_Ecov_a50.size(); i++) for(int j = 0; j <= i; j++) log_a50_Ecov(y,i) += beta_Ecov_a50(j) * Ecov_y(y) + a50_re(y);
    Type temp;
    for(int i = 0; i < 3; i++) 
    {
      if(i+1 <= beta_Ecov_k.size()) temp = exp(log_k_Ecov(y,i));
      else temp = exp(log_k_Ecov(y,beta_Ecov_k.size()-1));
      logit_p_Ecov(y,i) = temp;
    }
    for(int i = 0; i < 3; i++)
    {
      if(i+1 <= beta_Ecov_a50.size()) temp = exp(log_a50_Ecov(y,i));
      else temp = exp(log_a50_Ecov(y,beta_Ecov_a50.size()-1));
      logit_p_Ecov(y,i) *= Type(i+1) - temp;
    }
  }
  //see(n_obs);
  REPORT(mat);
  ADREPORT(Ecov_y);
  ADREPORT(Ecov_phi);
  ADREPORT(Ecov_sig);
  ADREPORT(log_a50);
  ADREPORT(log_k);
  ADREPORT(logit_pmat);
  ADREPORT(log_k_Ecov);
  ADREPORT(log_a50_Ecov);
  ADREPORT(logit_p_Ecov);
  REPORT(log_k);
  REPORT(log_a50);
  REPORT(log_k_Ecov);
  REPORT(log_a50_Ecov);
  REPORT(logit_p_Ecov);
  REPORT(nll_Ecov);
  REPORT(nll_Ecov_obs);
  REPORT(ll_mat);
  REPORT(pmat);
  
  return nll;
}

