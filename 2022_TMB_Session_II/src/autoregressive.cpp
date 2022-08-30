//adapted from Jim Thorson, https://github.com/James-Thorson/2018_FSH556/blob/master/Week%205%20--%201D%20spatial%20models/Lecture/autoregressive_V1.cpp
#include <TMB.hpp>

// trace of a matrix
template<class Type>
Type trace( matrix<Type> mat ){
  Type Return = 0;
  for(int i=0; i<mat.col(0).size(); i++) Return += mat(i,i);
  return Return;
}

template<class Type>
Type logdet( matrix<Type> mat ){
  int dim = mat.col(0).size();
  matrix<Type> chol(dim,dim);
  chol = mat.llt().matrixL();   /* L0 L0' = Q0 */
  Type logdet_mat = 0;
  for(int i=0; i<dim; i++ ) logdet_mat += Type(2.0) * log(chol(i,i));
  return logdet_mat;
}

template<class Type>
Type dmvnorm( vector<Type> x, matrix<Type> Q, int give_log=0 ){
  int n_x = x.size();
  Type logres = 0;
  vector<Type> Tmp_x(n_x);
  Type logdet_Q = logdet( Q );
  Tmp_x =  Q * x.matrix();
  logres += ( Type(0.5)*logdet_Q );
  logres += Type(-0.5) * (x * Tmp_x).sum();  //
  if (give_log) return logres; else return exp(logres);
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Settings
  DATA_VECTOR( Options_vec );
  // Slot 0: method for calculating probability of random effects

  // Data
  DATA_VECTOR( y );

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_sigma2 );
  PARAMETER( logit_rho );

  // Random effects
  PARAMETER_VECTOR( u );

  // Objective funcction
  int n_i = y.size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();
  Type sigma2 = exp(ln_sigma2);
  Type rho = 1 / (1 + exp(-logit_rho));

  // Probability of random effects
  using namespace density;
  
  //// Manually set-up AD1 structure
  if( Options_vec(0)==0 ){
    jnll_comp(1) -= dnorm( u(0), Type(0.0), pow(sigma2,0.5), true );
    for(int i=1; i<n_i; i++) jnll_comp(1) -= dnorm( u(i), rho*u(i-1), pow(sigma2,0.5), true );
  }
  //// Calculate using precision matrix
  if( Options_vec(0)==1 ){
    matrix<Type> Q_ii(n_i,n_i);
    Q_ii.setZero();
    for(int i=0; i<n_i; i++) Q_ii(i,i) = (1+pow(rho,2))/sigma2;
    for(int i=1; i<n_i; i++){
      Q_ii(i-1,i) = -rho/sigma2;
      Q_ii(i,i-1) = -rho/sigma2;
    }
    REPORT( Q_ii )
    jnll_comp(1) -= dmvnorm( u, Q_ii, true );
  }
  //// Calculate using precision matrix and built-in GMRF
  if( Options_vec(0)==2 ){
    SparseMatrix<Type> Q_ii(n_i,n_i);
    //Q_ii.setZero();
    for(int i=0; i<n_i; i++) Q_ii.coeffRef(i,i) = (1+pow(rho,2));
    for(int i=1; i<n_i; i++){
      Q_ii.coeffRef(i-1,i) = -rho;
      Q_ii.coeffRef(i,i-1) = -rho;
    }
    REPORT( Q_ii );
    jnll_comp(1) += SCALE( GMRF( Q_ii ), sqrt(sigma2) )( u );
  }
  //// Calculate using covariance matrix
  if( Options_vec(0)==3 ){
    matrix<Type> Cov_ii(n_i,n_i);
    for(int i1=0; i1<n_i; i1++){
    for(int i2=i1; i2<n_i; i2++){
      Cov_ii(i1,i2) = sigma2 / (1-pow(rho,2)) * pow( rho, double(i2-i1) );
      if(i1!=i2) Cov_ii(i2,i1) = Cov_ii(i1,i2);
    }}
    REPORT( Cov_ii );
    jnll_comp(1) += MVNORM(Cov_ii)( u );
  }
  //// Calculate using built-in TMB functions
  if( Options_vec(0)==4 ){
    jnll_comp(1) += SCALE( AR1(rho), pow(sigma2 / (1-pow(rho,2)),0.5))( u );
  }

  // Probability of data conditional on random effects
  for( int i=0; i<n_i; i++){
    jnll_comp(0) -= dpois( y(i), exp(beta0 + u(i)), true );
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( sigma2 );
  REPORT( rho );

  return jnll;
}
