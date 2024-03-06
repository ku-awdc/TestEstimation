#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER( P );
  DATA_ARRAY( tally );
  DATA_VECTOR( se_alpha );
  DATA_VECTOR( se_beta );
  DATA_VECTOR( sp_alpha );
  DATA_VECTOR( sp_beta );

  PARAMETER_VECTOR(selg);
  PARAMETER_VECTOR(splg);
  PARAMETER_VECTOR(prevlg);
  
  /*
  PARAMETER_VECTOR(covse);
  PARAMETER_VECTOR(covsp);
  */
  
  vector<Type> se = invlogit(selg);
  vector<Type> sp = invlogit(splg);
  vector<Type> prev = invlogit(prevlg);
  const Type one = static_cast<Type>(1.0);

  vector<Type> covse(1L);
  vector<Type> covsp(1L);
  covse[0L] = 0.0;
  covsp[0L] = 0.0;

  vector<Type> se_prob(4L);
  vector<Type> sp_prob(4L);
  
  // Note: could potentially use a beta dist to impose constraints??
  se_prob[0L] = (1.0-se[0L])*(1.0-se[1L]) +covse[0L];
  se_prob[1L] = se[0L]*(1.0-se[1L]) -covse[0L];
  se_prob[2L] = (1.0-se[0L])*se[1L] -covse[0L];
  se_prob[3L] = se[0L]*se[1L] +covse[0L];

  sp_prob[0L] = sp[0L]*sp[1L] +covsp[0L];
  sp_prob[1L] = (1.0-sp[0L])*sp[1L] -covsp[0L];
  sp_prob[2L] = sp[0L]*(1.0-sp[1L]) -covsp[0L];
  sp_prob[3L] = (1.0-sp[0L])*(1.0-sp[1L]) +covsp[0L];

  Type neglogL = 0.0;
  
  for (int p=0L; p<P; ++p)
  {
    Type pr = static_cast<Type>(prev[p]);
    vector<Type> obs = tally.col(p);
    vector<Type> probs = se_prob*pr + sp_prob*(one-pr);
    neglogL -= dmultinom(obs, probs, true);
  }
  
  // NB: doesn't work for implementing constraints (if x<0 or x>1 then result is ignored):
  neglogL -= sum(dbeta(se, se_alpha, se_beta, true));
  neglogL -= sum(dbeta(sp, sp_alpha, sp_beta, true));
  neglogL -= sum(dbeta(prev, static_cast<Type>(1.0), static_cast<Type>(1.0), true));
  
  return neglogL;
}
