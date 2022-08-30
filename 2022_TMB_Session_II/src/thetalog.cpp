// Theta logistic population model from Pedersen et al 2012, Ecol. Modelling.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data section */
  DATA_VECTOR(Y);
  DATA_INTEGER(method);
  /* Parameter section */
  PARAMETER_VECTOR(u);
  PARAMETER(logr0);
  PARAMETER(logtheta);
  PARAMETER(logK);
  PARAMETER(logSigu);
  PARAMETER(logSigy);
  /* Procedure section */
  Type r0=exp(logr0);
  Type theta=exp(logtheta);
  Type K=exp(logK);
  Type sigu=exp(logSigu);
  Type sigy=exp(logSigy);
  int timeSteps=Y.size();
  Type ans=0;
  if(method == 1){
    for(int i=1;i<timeSteps;i++){
      Type m=u[i-1]+r0*(1.0-pow(exp(u[i-1])/K,theta));
      ans-=dnorm(u[i],m,sigu,true);
    }
  }
  if(method == 2){
    for(int i=1;i<timeSteps;i++){
      Type m=u[i-1]+r0*(1.0-pow(exp(u[i-1])/K,theta));
    }

  }
  for(int i=0;i<timeSteps;i++){
    ans-=dnorm(Y[i],u[i],sigy,true);
  }
  return ans;
}
