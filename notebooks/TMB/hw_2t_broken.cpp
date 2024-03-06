#include <TMB.hpp>

constexpr int C_combos = 4L;

template<class Type>
bool invalid_probs(const std::vector<Type>& probs){
  bool rv = false;
  for (auto & prob : probs)
  {
    rv = rv || (prob < 0.0) || (prob > 1.0);
  }
  return rv;
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_SCALAR( INF );
  DATA_INTEGER( P );
  DATA_ARRAY( tally );
  /*
  DATA_VECTOR( se_alpha );
  DATA_VECTOR( se_beta );
  DATA_VECTOR( sp_alpha );
  DATA_VECTOR( sp_beta );
  */

  PARAMETER_VECTOR(se);
  PARAMETER_VECTOR(sp);
  /*
  PARAMETER_VECTOR(covse);
  PARAMETER_VECTOR(covsp);
  */
  PARAMETER_VECTOR(prev);
  
  vector<Type> covse(1L);
  vector<Type> covsp(1L);
  covse[0L] = 0.0;
  covsp[0L] = 0.0;

  vector<Type> se_prob(4L);
  vector<Type> sp_prob(4L);

  // TODO: use this rather than 1.0 ??
  // Type one = static_cast<Type>(1.0);
  
  // This returns zero:
  // return std::numeric_limits<Type>::infinity();
  
  /* Note: can't use if with Type, but ? : is fine:
  return trv;
  */
  
  /*
  Type trv = 0.0;
  trv += (sp[1L]>0.0) ? INF : 0.0;
  return trv;
  */

  Type neglogL = 0.0;
  
  //negLogL += invalid_probs<Type>(sp) ? INF : 0.0;
  //negLogL += (sp[1L]>1.0) ? INF : 0.0;
  
  vector<Type> sec = se;
  vector<Type> spc = sp;
  vector<Type> prc = prev;
  if(invalid_probs<Type>(sec)) return INF;
  if(invalid_probs<Type>(spc)) return INF;
  if(invalid_probs<Type>(prc)) return INF;
  
  vector<Type> youden(2L);
  for (int i=0L; i<2L; ++i)
  {
    youden[i] = se[i]+sp[i]-static_cast<Type>(1.0);
  }
  if (invalid_probs<Type>(youden)) return INF;
  

  se_prob[0L] = (1.0-se[0L])*(1.0-se[1L]) +covse[0L];
  se_prob[1L] = se[0L]*(1.0-se[1L]) -covse[0L];
  se_prob[2L] = (1.0-se[0L])*se[1L] -covse[0L];
  se_prob[3L] = se[0L]*se[1L] +covse[0L];
  if (invalid_probs<Type>(se_prob)) return INF;

  sp_prob[0L] = sp[0L]*sp[1L] +covsp[0L];
  sp_prob[1L] = (1.0-sp[0L])*sp[1L] -covsp[0L];
  sp_prob[2L] = sp[0L]*(1.0-sp[1L]) -covsp[0L];
  sp_prob[3L] = (1.0-sp[0L])*(1.0-sp[1L]) +covsp[0L];
  if (invalid_probs<Type>(sp_prob)) return INF;

  for (int p=0L; p<P; ++p)
  {
    vector<Type> probs(4L);
    vector<Type> obs = tally.col(p);
    Type pr = static_cast<Type>(prev[p]);
    for(int i=0L; i<4L; ++i)
    {
      probs[i] = se_prob[i]*pr + sp_prob[i]*(1-pr);
    }
    neglogL -= dmultinom(obs, probs, true);
    //neglogL -= dbinom(obs[3L], obs[0L]+obs[3L], pr, true);
  }

  /*
  // Doesn't work:
  neglogL = (sp[1L]>1.0) ? INF : neglogL;
  neglogL = (sp[0L]>1.0) ? INF : neglogL;
  */

  return neglogL;
}
