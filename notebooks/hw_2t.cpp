#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER( N );
  DATA_INTEGER( P );
  DATA_ARRAY( tally );
  DATA_VECTOR( se_alpha );
  DATA_VECTOR( se_beta );
  DATA_VECTOR( sp_alpha );
  DATA_VECTOR( sp_beta );

  PARAMETER_VECTOR(se);
  PARAMETER_VECTOR(sp);
  PARAMETER_VECTOR(covse);
  PARAMETER_VECTOR(covsp);
  PARAMETER_VECTOR(prev);

  vector<Type> se_prob(4L);
  vector<Type> sp_prob(4L);

  // TODO: use this rather than 1.0 ??
  Type one = static_cast<Type>(1.0);

  se_prob[0L] = (1.0-se[0L])*(1.0-se[1L]) +covse[0L];
  se_prob[1L] = se[0L]*(1.0-se[1L]) -covse[0L];
  se_prob[2L] = (1.0-se[0L])*se[1L] -covse[0L];
  se_prob[3L] = se[0L]*se[1L] +covse[0L];
  // if(!probs_ok(se_prob)) return(std::numeric_limits<T>::infinity());

  sp_prob[0L] = sp[0L]*sp[1L] +covsp[0L];
  sp_prob[1L] = (1.0-sp[0L])*sp[1L] -covsp[0L];
  sp_prob[2L] = sp[0L]*(1.0-sp[1L]) -covsp[0L];
  sp_prob[3L] = (1.0-sp[0L])*(1.0-sp[1L]) +covsp[0L];
  
  Type neglogL = 0.0;
  for(int p=0L; p<P; ++p)
  {
    vector<Type> probs(4L);
    vector<Type> obs(4L);
    Type pr = static_cast<Type>(prev[p]);
    for(int i=0L; i<4L; ++i)
    {
      probs[i] = se_prob[i]*pr + sp_prob[i]*(1-pr);
      obs[i] = tally[i,p];
    }
    // Works:
    neglogL -= dmultinom(obs, probs, true);
  }

  return neglogL;
}
